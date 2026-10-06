module CarrierLoopStagingTest

# The host side of TrackingLoops' carrier-loop staging: the fold hands every
# record the discriminators its bit buffer allows (a four-quadrant FLL on a
# wiped-off prompt, a four-quadrant PLL with the secondary-code sync's sign),
# drops the previous prompt where the FLL must not compare across a record, and
# reports the carrier phase's half-cycle state. The staging itself is tested in
# TrackingLoops.

using Test: @test, @testset, @inferred
using Unitful: Hz
using Random: Xoshiro
using Statistics: mean
using GNSSSignals:
    GPSL1CA,
    GPSL5Q,
    GalileoE1C,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency,
    get_code_length,
    get_secondary_code_length
import Tracking
using Tracking:
    TrackedSignal,
    TrackState,
    add_satellite!,
    get_carrier_doppler,
    get_carrier_phase,
    get_carrier_phase_polarity,
    get_sat_state,
    track!
import TrackingLoops
using TrackingLoops:
    BitBuffer,
    CorrelatorOutput,
    get_default_correlator,
    has_bit_or_secondary_code_been_found,
    sync_polarity

# Track a pilot, noise-free unless `cn0` is given. Returns the
# replica-minus-signal carrier phase after every chunk, the chunk the
# four-quadrant PLL took over at (the sync), and the state.
function track_pilot(
    signal,
    start_phase;
    prn = 1,
    duration_chunks,
    chunk_ms,
    cn0 = nothing,
    doppler_error = 5Hz,
    seed = 1,
)
    carrier_doppler = 200Hz
    sampling_frequency = max(2 * get_code_frequency(signal), 25e6Hz)
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(signal) +
        get_code_frequency(signal)
    chunk_time = chunk_ms / 1000Hz
    num_samples = round(Int, sampling_frequency * chunk_time)
    code_length = get_code_length(signal) * get_secondary_code_length(signal)
    rng = Xoshiro(seed)
    track_state = TrackState(; signal)
    add_satellite!(
        track_state;
        prn,
        code_phase = 0.0,
        carrier_doppler = carrier_doppler - doppler_error,
    )
    t = (0:(num_samples-1)) ./ sampling_frequency
    phase_offsets = Float64[]
    switched_at = nothing
    for i = 0:(duration_chunks-1)
        start_time = i * chunk_time
        code = gen_code(
            num_samples,
            signal,
            prn,
            sampling_frequency,
            code_frequency,
            mod(code_frequency * start_time, code_length),
        )
        carrier_phase = 2π * carrier_doppler * start_time + start_phase
        samples =
            cis.(2π .* carrier_doppler .* t .+ carrier_phase) .* code ./
            sqrt(mean(abs2, code))
        if !isnothing(cn0)
            fs = Float64(sampling_frequency / 1Hz)
            samples =
                10^(cn0 / 20) .* samples .+ randn(rng, ComplexF64, num_samples) .* sqrt(fs)
        end
        track!(ComplexF32.(samples), track_state, sampling_frequency)
        end_phase = 2π * carrier_doppler * (start_time + chunk_time) + start_phase
        push!(
            phase_offsets,
            rem2pi(get_carrier_phase(track_state) - end_phase, RoundNearest),
        )
        driver = first(get_sat_state(track_state, prn).signals)
        if isnothing(switched_at) && sync_polarity(signal, driver.bit_buffer, prn) != 0
            switched_at = i + 1
        end
    end
    phase_offsets, switched_at, track_state
end

# Half-cycle steps of the carrier phase from chunk `from` on.
function half_cycle_slips(phase_offsets, from)
    cumulative = phase_offsets[from]
    reference = round(cumulative / π)
    slips = 0
    for k = (from+1):length(phase_offsets)
        cumulative += rem2pi(phase_offsets[k] - phase_offsets[k-1], RoundNearest)
        step = round(cumulative / π)
        step != reference && (slips += 1; reference = step)
    end
    slips
end

@testset "$(nameof(typeof(signal))) switches to the four-quadrant PLL at the sync" for (
    signal,
    chunk_ms,
    duration_chunks,
) in (
    (GPSL5Q(), 1, 2000),
    (GalileoE1C(), 4, 500),
)
    # Both start phases, so that Costas pulls in on either sign of the prompt.
    # The switch may jump half a cycle; from shortly after it, the phase is
    # continuous.
    for start_phase in (0.0, π)
        phase_offsets, switched_at, track_state =
            track_pilot(signal, start_phase; duration_chunks, chunk_ms)
        @test has_bit_or_secondary_code_been_found(track_state)
        @test !isnothing(switched_at)
        @test abs(get_carrier_doppler(track_state) - 200Hz) < 0.1Hz
        settled = max(switched_at, round(Int, 300 / chunk_ms)) + round(Int, 200 / chunk_ms)
        @test half_cycle_slips(phase_offsets, settled) == 0
    end
end

@testset "Carrier phase polarity" begin
    # Unresolved until the driver pilot syncs.
    track_state = TrackState(; signal = GalileoE1C())
    add_satellite!(track_state; prn = 11, code_phase = 0.0, carrier_doppler = 0.0Hz)
    @test @inferred(get_carrier_phase_polarity(track_state)) === 0
    @test get_carrier_phase_polarity(get_sat_state(track_state, 11)) === 0
    @test Tracking._carrier_phase_polarity(nothing, nothing) === Int8(0)

    # Resolved after: corrected, the carrier phase is the signal's, whichever
    # sign Costas pulled in on.
    for (signal, chunk_ms, duration_chunks) in
        ((GPSL5Q(), 1, 2000), (GalileoE1C(), 4, 500)),
        start_phase in (0.0, π)

        phase_offsets, _, track_state =
            track_pilot(signal, start_phase; duration_chunks, chunk_ms)
        polarity = get_carrier_phase_polarity(track_state)
        @test polarity != 0
        corrected = phase_offsets[end] + (polarity < 0 ? π : 0.0)
        @test abs(rem2pi(corrected, RoundNearest)) < 0.2
    end
end

# A synced bit buffer, as a data signal's after its bit sync.
function _synced(b::BitBuffer)
    BitBuffer(
        b.code_block_buffer,
        b.code_block_buffer_length,
        true,
        0,
        Int8(1),
        b.prompt_accumulator,
        0,
        b.soft_bits,
        b.phase_acc,
    )
end

@testset "No FLL reading across a change of record length" begin
    fs = 5e6Hz
    signal = GPSL1CA()
    unsynced = TrackedSignal(signal)
    ts = TrackedSignal(
        unsynced;
        bit_buffer = _synced(unsynced.bit_buffer),
        last_fully_integrated_filtered_prompt = cis(0.1),
        last_fully_integrated_num_code_blocks = 1,
    )
    record(n) = CorrelatorOutput(get_default_correlator(signal), n, n)
    # One block after one block: the previous prompt holds.
    @test Tracking._fll_previous_prompt(ts, record(5000), fs) == cis(0.1)
    # Twenty after one, as for a data signal at its bit sync: none.
    @test iszero(Tracking._fll_previous_prompt(ts, record(100000), fs))
end

@testset "A pilot drops its previous prompt at the sync that wipes it off" begin
    # The four-quadrant FLL must not compare the first wiped-off prompt with one
    # correlated before the sync.
    pilot = TrackedSignal(
        TrackedSignal(GalileoE1C());
        last_fully_integrated_filtered_prompt = cis(0.2),
    )
    synced_pilot = TrackedSignal(pilot; bit_buffer = _synced(pilot.bit_buffer))
    reset = Tracking._reset_on_snap(pilot, synced_pilot)
    @test iszero(reset.last_fully_integrated_filtered_prompt)
    @test reset.integrated_samples == 0
    # A data signal's sync wipes nothing off: its previous prompt stays.
    data = TrackedSignal(
        TrackedSignal(GPSL1CA());
        last_fully_integrated_filtered_prompt = cis(0.2),
    )
    synced_data = TrackedSignal(data; bit_buffer = _synced(data.bit_buffer))
    @test Tracking._reset_on_snap(data, synced_data).last_fully_integrated_filtered_prompt ==
          cis(0.2)
end

end
