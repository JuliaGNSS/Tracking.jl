module SignalCombiningTest

# The host side of TrackingLoops' signal combining: the driver fold hands every
# passenger record to the estimator before the driver record it ends within,
# leaves the rest pending, drops them at the sync phase snap, and refers the
# passengers' code loop to the driver by their group delays. The combining
# itself is tested in TrackingLoops.

using Test: @test, @testset
using Unitful: Hz, ns, s, ustrip
using StaticArrays: SVector
using Random: Xoshiro
using Statistics: mean, std
using GNSSSignals:
    GalileoE1B,
    GalileoE1C,
    gen_code,
    get_center_frequency,
    get_code_frequency,
    get_code_length,
    get_secondary_code_length
import Tracking
using Tracking:
    TrackState,
    add_satellite!,
    append_correlator_output!,
    estimate_dopplers_and_filter_prompt!,
    get_code_phase,
    get_doppler_estimator_state,
    get_sat_state,
    set_group_delay!,
    track!
using TrackingLoops:
    BitBuffer,
    ConventionalAssistedPLLAndDLL,
    CorrelatorOutput,
    SignalCombiningSums,
    get_default_correlator,
    update_accumulator

# A VEML correlator (Galileo E1's default) with the taps around prompt `p`.
_veml(p; early = 0.7, late = 0.7) = update_accumulator(
    get_default_correlator(GalileoE1B()),
    p .* SVector(0.2, early, 1.0, late, 0.2),
)

_e1_state(signals = (GalileoE1C(), GalileoE1B()); combining = true) = TrackState(;
    signals = (e1 = signals,),
    doppler_estimator = ConventionalAssistedPLLAndDLL(; combine_signals = combining),
)

# Mark a signal synced up front. A signal with one symbol per code period
# (Galileo E1B) syncs on its first record, which would make that chunk the sync
# phase snap that drops pending passenger sums.
function _mark_synced!(track_state, group, prn, index)
    sats = Tracking.get_sat_states(track_state, group)
    signal = sats[prn].signals[index]
    b = signal.bit_buffer
    synced = BitBuffer(
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
    sats[prn] = Tracking.TrackedSat(
        sats[prn];
        signals = Base.setindex(
            sats[prn].signals,
            Tracking.TrackedSignal(signal; bit_buffer = synced),
            index,
        ),
    )
    track_state
end

_sums(track_state) =
    get_doppler_estimator_state(get_sat_state(track_state, :e1, 11)).signal_combining_sums

const FS = 16.368e6Hz
const N = 65472  # 4 ms
_record(p, sample_index) = CorrelatorOutput(_veml(p * N), N, sample_index)

@testset "Passenger records are combined into the driver record they end within" begin
    track_state = _e1_state()
    add_satellite!(
        track_state;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    _mark_synced!(track_state, :e1, 11, 2)
    estimate!() = estimate_dopplers_and_filter_prompt!(track_state, (L1 = FS,))

    # A chunk with a passenger record only: it stays pending. Group delays are
    # unknown, so it is combined into the carrier loops only (the FLL has no
    # previous prompt).
    append_correlator_output!(track_state, _record(cis(0.1), N), :e1, 11, GalileoE1B)
    estimate!()
    @test _sums(track_state).pll.weight ≈ 0.5 * 0.004s
    @test _sums(track_state).pll.sum / _sums(track_state).pll.weight ≈ 0.1 / 2π
    @test iszero(_sums(track_state).fll.weight)
    @test iszero(_sums(track_state).dll.weight)
    @test isempty(Tracking.get_correlator_outputs(track_state, :e1, 11, GalileoE1B))

    # The driver's next record consumes it with the passenger record ending by
    # then; one ending later waits for the following driver record.
    append_correlator_output!(track_state, _record(cis(0.1), N), :e1, 11, GalileoE1B)
    append_correlator_output!(track_state, _record(cis(0.0), 2N), :e1, 11, GalileoE1C)
    append_correlator_output!(track_state, _record(cis(0.3), 3N), :e1, 11, GalileoE1B)
    estimate!()
    @test _sums(track_state).pll.weight ≈ 0.5 * 0.004s
    @test _sums(track_state).pll.sum / _sums(track_state).pll.weight ≈ 0.3 / 2π
    # This one had a previous prompt of the passenger's own, of the same length.
    @test _sums(track_state).fll.weight > 0s^3

    # With both group delays known, the passenger joins the code loop too.
    set_group_delay!(track_state, :e1, 11, GalileoE1C, 0.0ns)
    set_group_delay!(track_state, :e1, 11, GalileoE1B, 0.0ns)
    append_correlator_output!(track_state, _record(cis(0.2), 4N), :e1, 11, GalileoE1B)
    append_correlator_output!(track_state, _record(cis(0.0), 4N), :e1, 11, GalileoE1C)
    append_correlator_output!(track_state, _record(cis(0.2), 5N), :e1, 11, GalileoE1B)
    estimate!()
    @test _sums(track_state).dll.weight ≈ 0.5 * 0.004s

    # Without combining nothing is pending.
    uncombined = _e1_state(; combining = false)
    add_satellite!(
        uncombined;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    append_correlator_output!(uncombined, _record(cis(0.1), N), :e1, 11, GalileoE1B)
    estimate_dopplers_and_filter_prompt!(uncombined, (L1 = FS,))
    @test _sums(uncombined) === SignalCombiningSums()
    @test isempty(Tracking.get_correlator_outputs(uncombined, :e1, 11, GalileoE1B))
end

@testset "The sync phase snap drops pending passenger sums" begin
    # Galileo E1B driving syncs on its first record; the E1C record that ends
    # after it would stay pending, but the snap drops it with the driver's
    # in-flight integration.
    function pending_after_first_fold(; driver_synced)
        track_state = _e1_state((GalileoE1B(), GalileoE1C()))
        add_satellite!(
            track_state;
            group = :e1,
            prn = 11,
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        driver_synced && _mark_synced!(track_state, :e1, 11, 1)
        append_correlator_output!(track_state, _record(cis(0.0), N), :e1, 11, GalileoE1B)
        append_correlator_output!(track_state, _record(cis(0.1), 2N), :e1, 11, GalileoE1C)
        estimate_dopplers_and_filter_prompt!(track_state, (L1 = FS,))
        _sums(track_state)
    end
    @test pending_after_first_fold(; driver_synced = false) === SignalCombiningSums()
    # Without a snap it stays pending.
    @test pending_after_first_fold(; driver_synced = true).pll.weight > 0s
end

# Galileo E1 with E1C (pilot) driving and E1B (data) combined, E1B delayed by
# `delay_chips` against E1C. Returns the driver's code-phase error per chunk.
function track_e1(; combining, delays, delay_chips = 0.0, cn0 = 45.0, duration, seed = 1)
    fs = 16.368e6
    n = round(Int, fs * 0.004)
    e1b, e1c = GalileoE1B(), GalileoE1C()
    fd = 1000.0
    prn = 11
    fc =
        ustrip(Hz, get_code_frequency(e1c)) *
        (1 + fd / ustrip(Hz, get_center_frequency(e1c)))
    rng = Xoshiro(seed)
    track_state = TrackState(;
        signals = (e1 = (e1c, e1b),),
        doppler_estimator = ConventionalAssistedPLLAndDLL(;
            code_loop_filter_bandwidth = 4.0Hz,
            combine_signals = combining,
        ),
    )
    add_satellite!(
        track_state;
        group = :e1,
        prn,
        code_phase = 0.0,
        carrier_doppler = (fd - 5) * Hz,
    )
    if !isnothing(delays)
        set_group_delay!(track_state, :e1, prn, GalileoE1C, delays[1])
        set_group_delay!(track_state, :e1, prn, GalileoE1B, delays[2])
    end
    t = (0:(n-1)) ./ fs
    errors = Float64[]
    code_length = get_code_length(e1c) * get_secondary_code_length(e1c)
    # Each component at half the power.
    amplitude = 10^(cn0 / 20) / sqrt(2)
    for i = 0:(round(Int, duration/0.004)-1)
        t0 = i * 0.004
        c = gen_code(n, e1c, prn, fs * Hz, fc * Hz, mod(fc * t0, code_length))
        b = gen_code(n, e1b, prn, fs * Hz, fc * Hz, mod(fc * t0 - delay_chips, 4092))
        bit = rand(rng, (-1.0, 1.0))
        samples =
            amplitude .* cis.(2π * fd .* (t .+ t0)) .*
            (c ./ sqrt(mean(abs2, c)) .+ bit .* b ./ sqrt(mean(abs2, b))) .+
            randn(rng, ComplexF64, n) .* sqrt(fs)
        track!(ComplexF32.(samples), track_state, fs * Hz)
        error = get_code_phase(track_state, :e1, prn) - mod(fc * (t0 + 0.004), 4092)
        push!(errors, mod(error + 2046, 4092) - 2046)
    end
    errors
end

@testset "Galileo E1B combined into E1C's loops end to end" begin
    # E1B arrives 0.1 chip after E1C. With both group delays the code loop stays
    # on the driver's code phase (a residual of the discriminators' CBOC
    # calibration remains); without them E1B is not combined into it; with a zero
    # difference instead of the true one it pulls the code phase half way.
    delay = 0.1 / 1.023e6 * 1e9 * ns
    bias(delays) = mean(
        track_e1(; combining = true, delays, delay_chips = 0.1, duration = 2.0)[(end-249):end],
    )
    @test abs(bias((0.0ns, delay))) < 0.01
    @test abs(bias(nothing)) < 0.005
    @test bias((0.0ns, 0.0ns)) ≈ -0.05 atol = 0.015

    # A second signal's code discriminator lowers the code-phase noise.
    noise(combining) = std(
        track_e1(; combining, delays = (0.0ns, 0.0ns), cn0 = 35.0, duration = 6.0)[500:end],
    )
    @test noise(true) < 0.7 * noise(false)
end

# Records are built outside the measurement, which covers the fold alone.
function estimate_allocations(combining)
    track_state = _e1_state(; combining)
    add_satellite!(
        track_state;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    set_group_delay!(track_state, :e1, 11, GalileoE1C, 0.0ns)
    set_group_delay!(track_state, :e1, 11, GalileoE1B, 0.0ns)
    driver_record = _record(cis(0.1), N)
    passenger_record = _record(cis(0.2), N)
    sampling_frequencies = (L1 = FS,)
    allocated = 0
    for _ = 1:10
        append_correlator_output!(track_state, driver_record, :e1, 11, GalileoE1C)
        append_correlator_output!(track_state, passenger_record, :e1, 11, GalileoE1B)
        allocated = @allocated estimate_dopplers_and_filter_prompt!(
            track_state,
            sampling_frequencies,
        )
    end
    allocated
end

@testset "Signal combining adds no allocation" begin
    @test estimate_allocations(true) <= estimate_allocations(false)
end

end
