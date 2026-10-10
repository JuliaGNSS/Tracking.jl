module PassengersTest

# Every record of a satellite goes through the estimator's `step_loop`, the
# passengers' (`signals[2:end]`) as well as the driver's, each built from its
# signal's state: Tracking never asks which estimator it drives. With
# TrackingLoops' `VectorPLLAndDLL(pilot => data)` that is what lets the vector
# loop range on the pilot and decode the data component.

using Test: @test, @testset
using Logging: NullLogger, with_logger
using Unitful: Hz
using StaticArrays: SVector
using GNSSSignals: GalileoE1B, GalileoE1C, get_carrier_phase_offset
using Tracking:
    TrackState,
    add_satellite!,
    append_correlator_output!,
    estimate_dopplers_and_filter_prompt!,
    get_carrier_doppler,
    get_doppler_estimator_state,
    get_sat_state
using TrackingLoops:
    AbstractDopplerEstimator,
    ConventionalAssistedPLLAndDLL,
    CorrelatorOutput,
    LoopRecord,
    VectorPLLAndDLL,
    get_default_correlator,
    satellite_report,
    update_accumulator
import TrackingLoops

const FS = 16.368e6Hz
const N = 65472 # 4 ms, one Galileo E1 code period

# An E1 record ending at `sample_index`, its prompt `p` on the VEML taps.
_record(p, sample_index) = CorrelatorOutput(
    update_accumulator(
        get_default_correlator(GalileoE1B()),
        ComplexF64(p * N) .* SVector(0.2, 0.7, 1.0, 0.7, 0.2),
    ),
    N,
    sample_index,
)

# An estimator that remembers the signal and summary of every record it is
# stepped with, then runs the conventional loop.
struct RecordingEstimator{E} <: AbstractDopplerEstimator
    inner::E
    seen::Vector{Tuple{Symbol,Int,Bool,Int}}
end
RecordingEstimator() =
    RecordingEstimator(ConventionalAssistedPLLAndDLL(), Tuple{Symbol,Int,Bool,Int}[])
TrackingLoops.init_estimator_state(e::RecordingEstimator, signal, carrier, code) =
    TrackingLoops.init_estimator_state(e.inner, signal, carrier, code)
function TrackingLoops.step_loop(
    e::RecordingEstimator,
    state,
    record::LoopRecord,
    words,
    landing::Int64,
)
    push!(
        e.seen,
        (
            TrackingLoops.get_signal_id(record.signal),
            record.sample_index,
            record.bit_synced,
            length(record.new_soft_bits),
        ),
    )
    TrackingLoops.step_loop(e.inner, state, record, words, landing)
end

# One chunk's records of satellite 11: E1B's (data, a sign per 4 ms symbol, in its
# carrier-phase frame) and E1C's (pilot), in the given order.
function feed!(track_state, k; data_first = true, carrier = 0.0)
    bit = isodd(hash(k) >> 3) ? 1.0 : -1.0
    rotation =
        cis(get_carrier_phase_offset(GalileoE1B()) - get_carrier_phase_offset(GalileoE1C()))
    pilot_prompt = cis(carrier * k)
    data =
        () -> append_correlator_output!(
            track_state,
            _record(rotation * bit * pilot_prompt, k * N),
            :e1,
            11,
            GalileoE1B,
        )
    pilot =
        () -> append_correlator_output!(
            track_state,
            _record(pilot_prompt, k * N),
            :e1,
            11,
            GalileoE1C,
        )
    data_first ? (data(); pilot()) : (pilot(); data())
    # No noise observations are fed, so the noise-referenced C/N₀ warns once;
    # nothing here reads it.
    with_logger(NullLogger()) do
        estimate_dopplers_and_filter_prompt!(track_state, (L1 = FS,))
    end
end

function e1_state(estimator)
    track_state = TrackState(;
        signals = (e1 = (GalileoE1C(), GalileoE1B()),),
        doppler_estimator = estimator,
    )
    add_satellite!(
        track_state;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    track_state
end

@testset "Every record of a satellite is stepped, each with its signal's state" begin
    estimator = RecordingEstimator()
    track_state = e1_state(estimator)
    for k = 1:20
        feed!(track_state, k)
    end
    signals = first.(estimator.seen)
    @test count(==(:GalileoE1C), signals) == 20
    @test count(==(:GalileoE1B), signals) == 20
    # The driver's records first, then the passenger's, each signal's in order.
    @test signals[1:2] == [:GalileoE1C, :GalileoE1B]
    for id in (:GalileoE1C, :GalileoE1B)
        @test issorted(getindex.(filter(s -> s[1] == id, estimator.seen), 2))
    end
    # E1B, one symbol per code period, syncs on its first record and appends one
    # soft bit per record, which its record hands on exactly once.
    data = filter(s -> s[1] == :GalileoE1B, estimator.seen)
    @test all(s -> s[3], data)
    @test sum(s -> s[4], data) == 20
end

@testset "A passenger's records leave the driver's loop as it is" begin
    # The same pilot records, with and without the data component's records next to
    # them: the Dopplers the loop commands are identical.
    dopplers(signals) = begin
        track_state = TrackState(;
            signals = (e1 = signals,),
            doppler_estimator = ConventionalAssistedPLLAndDLL(),
        )
        add_satellite!(
            track_state;
            group = :e1,
            prn = 11,
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        map(1:30) do k
            feed_pilot_only = length(signals) == 1
            if feed_pilot_only
                append_correlator_output!(
                    track_state,
                    _record(cis(0.05k), k * N),
                    :e1,
                    11,
                    GalileoE1C,
                )
                with_logger(NullLogger()) do
                    estimate_dopplers_and_filter_prompt!(track_state, (L1 = FS,))
                end
            else
                feed!(track_state, k; carrier = 0.05)
            end
            get_carrier_doppler(get_sat_state(track_state, :e1, 11))
        end
    end
    @test dopplers((GalileoE1C(),)) == dopplers((GalileoE1C(), GalileoE1B()))
end

@testset "A vector loop ranges on E1C and decodes E1B, in either record order" begin
    for data_first in (true, false)
        estimator = VectorPLLAndDLL(GalileoE1C() => GalileoE1B(); approximate_year = 2021)
        track_state = e1_state(estimator)
        for k = 1:60
            feed!(track_state, k; data_first)
        end
        report = satellite_report(estimator, GalileoE1C(), 11)
        @test report.tracked && report.bit_synced
        @test report.decoder.prn == 11
        @test report.decoder isa typeof(TrackingLoops.GNSSDecoderState(GalileoE1B(), 11))
        @test isnothing(satellite_report(estimator, GalileoE1B(), 11))
        @test estimator.navigation.registrations == 1
        @test get_doppler_estimator_state(get_sat_state(track_state, :e1, 11)).slot == 1
    end
end

end
