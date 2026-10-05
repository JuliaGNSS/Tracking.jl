module CarrierLoopStagingTest

using Test: @test, @testset, @inferred
using Unitful: Hz, ms, s
using StaticArrays: SVector
using Random: Xoshiro
using Statistics: mean
using GNSSSignals:
    GPSL1CA,
    GPSL2CL,
    GPSL5I,
    GPSL5Q,
    GalileoE1C,
    GalileoE5aQP,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency,
    get_code_length,
    get_secondary_code_length
using TrackingLoopFilters: ThirdOrderAssistedBilinearLF, ThirdOrderBilinearLF
import Tracking
using Tracking:
    BitBuffer,
    FrequencyLockIndicator,
    EarlyPromptLateCorrelator,
    SatConventionalPLLAndDLL,
    SatVectorPLLAndDLL,
    TrackState,
    add_satellite!,
    frequency_lock_threshold,
    frequency_lock_window,
    get_carrier_doppler,
    get_carrier_phase,
    get_carrier_phase_polarity,
    get_doppler_estimator_state,
    get_sat_state,
    has_bit_or_secondary_code_been_found,
    track!

_prompt_correlator(p) = EarlyPromptLateCorrelator(SVector(p / 2, p, p / 2), 0.5)

include("frequency_lock_test_helpers.jl")  # `_locked`

@testset "Wiped-off prompts and the sync's sign" begin
    wiped_off = Tracking._is_wiped_off
    # Data signals never are.
    @test !wiped_off(GPSL1CA(), false)
    @test !wiped_off(GPSL1CA(), true)
    # Pilots with a secondary code once synced to it ...
    @test !wiped_off(GalileoE1C(), false)
    @test wiped_off(GalileoE1C(), true)
    @test !wiped_off(GPSL5Q(), false)
    @test wiped_off(GPSL5Q(), true)
    # ... and pilots without one from the start.
    @test wiped_off(GalileoE5aQP(), false)
    @test wiped_off(GPSL2CL(), false)

    # The sync polarity times secondary chip 0, which the pre-sync replica
    # carries on every block: +1 for GPS L5Q, -1 for every Galileo E1C PRN.
    synced(polarity) = Tracking.Accessors.setproperties(
        BitBuffer{UInt32}(),
        (; found = true, polarity = Int8(polarity)),
    )
    sync_sign = Tracking._sync_polarity
    @test @inferred(sync_sign(GPSL5Q(), synced(1), 1)) === Int8(1)
    @test sync_sign(GPSL5Q(), synced(-1), 1) === Int8(-1)
    @test sync_sign(GalileoE1C(), synced(1), 1) === Int8(-1)
    @test sync_sign(GalileoE1C(), synced(-1), 11) === Int8(1)
    # None where there is no such sync: before it, on data, and on the pilots
    # without a secondary code, which keep the Costas PLL.
    @test sync_sign(GalileoE1C(), BitBuffer{UInt32}(), 1) === Int8(0)
    @test sync_sign(GPSL5I(), synced(1), 1) === Int8(0)
    @test sync_sign(GalileoE5aQP(), synced(1), 1) === Int8(0)
    @test sync_sign(GPSL2CL(), BitBuffer{UInt32}(), 1) === Int8(0)
end

@testset "Frequency lock indicator" begin
    update = Tracking._update_frequency_lock
    signal = GPSL5Q()
    T = 1 / 1000Hz
    # Records until lock is declared, for a constant FLL reading.
    function records_to_lock(fll; integration_time = T, max_records = 10_000)
        indicator = FrequencyLockIndicator()
        for n = 1:max_records
            indicator = @inferred update(indicator, signal, fll, cis(0.0), integration_time)
            indicator.locked && return n
        end
        nothing
    end
    # Decided at the end of each 0.5 s window.
    @test records_to_lock(1.0Hz) in 499:501
    @test records_to_lock(-2.9Hz) in 499:501
    @test records_to_lock(1.0Hz; integration_time = 1 / 50Hz) in 24:26
    @test isnothing(records_to_lock(3.5Hz))
    @test isnothing(records_to_lock(-3.5Hz))

    # On long records the window spans at least four of them, and the threshold
    # shrinks to a quarter of the FLL's ±1/(4T) range, which 3 Hz would exceed:
    # GPS L2 CL's 1.5 s records read at most 0.17 Hz.
    @test frequency_lock_window(signal, T) == 0.5s
    @test frequency_lock_threshold(signal, T) == 3.0Hz
    @test frequency_lock_window(GPSL2CL(), 1.5s) == 6.0s
    @test frequency_lock_threshold(GPSL2CL(), 1.5s) ≈ 1 / 24 * 1Hz
    @test isnothing(records_to_lock(0.1Hz; integration_time = 1.5s, max_records = 100))
    @test records_to_lock(0.03Hz; integration_time = 1.5s) == 4

    # It is the window mean that counts: noise around zero averages out, a
    # residual error does not.
    function after_window(fll)
        indicator = FrequencyLockIndicator()
        for k = 1:500
            indicator = update(indicator, signal, fll(k) * 1Hz, cis(0.0), T)
        end
        indicator
    end
    @test after_window(k -> isodd(k) ? 20.0 : -20.0).locked
    @test after_window(k -> isodd(k) ? 24.0 : -16.0) === FrequencyLockIndicator()

    # A record without a previous prompt has no FLL reading, and lock latches.
    @test update(FrequencyLockIndicator(), signal, 0.0Hz, 0.0im, T) ===
          FrequencyLockIndicator()
    @test update(_locked(), signal, 10.0Hz, cis(0.0), T) === _locked()
end

# Run one record through the per-sat loop closure; returns the carrier update
# and the new state. The integration time is in 1/Hz, as the fold's
# `integrated_samples / sampling_frequency` is.
function close_loops(
    state,
    prompt,
    previous_prompt,
    wiped_off = false,
    polarity = 0;
    integration_time = 1 / 1000Hz,
)
    carrier_freq_update, _, state = Tracking._close_loops(
        state,
        GPSL5Q(),
        _prompt_correlator(prompt),
        previous_prompt,
        0.0Hz,
        5e6Hz,
        integration_time,
        18.0Hz,
        1.0Hz,
        wiped_off,
        Int8(polarity),
    )
    carrier_freq_update, state
end

@testset "Frequency-locked FLL-assisted filter is the pure PLL" begin
    state(carrier_loop_filter) = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter,
        frequency_lock = _locked(),
    )
    assisted, plain = state(ThirdOrderAssistedBilinearLF()), state(ThirdOrderBilinearLF())
    previous_prompt = cis(0.0)
    for k = 1:20
        prompt = cis(0.3 * sin(k))
        assisted_update, assisted = close_loops(assisted, prompt, previous_prompt)
        plain_update, plain = close_loops(plain, prompt, previous_prompt)
        @test assisted_update == plain_update
        previous_prompt = prompt
    end
    @test assisted.carrier_loop_filter.x1 == plain.carrier_loop_filter.x1
    @test assisted.carrier_loop_filter.x2 == plain.carrier_loop_filter.x2

    # Before lock the FLL branch reads the FLL discriminator.
    unlocked = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter = ThirdOrderAssistedBilinearLF(),
    )
    with_fll, _ = close_loops(unlocked, cis(0.3), cis(0.0))
    without_fll, _ = close_loops(unlocked, cis(0.3), 0.0im)
    @test with_fll != without_fll
end

@testset "Scalar loop: four-quadrant PLL from the sync, FLL dropped at frequency lock" begin
    # A closed loop on a synced pilot's prompt, locked half a cycle off, which
    # the sync's sign of -1 accounts for.
    integration_time = 1 / 1000Hz
    true_doppler = 5.0Hz
    state = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter = ThirdOrderAssistedBilinearLF(),
    )
    φ = 0.4  # signal minus replica phase, rad
    previous_prompt = 0.0im
    max_error = 0.0
    freq_update = 0.0Hz
    for k = 1:3000
        prompt = -cis(φ)
        freq_update, state = close_loops(state, prompt, previous_prompt, true, -1)
        state.frequency_lock.locked &&
            (max_error = max(max_error, abs(rem2pi(φ, RoundNearest))))
        previous_prompt = prompt
        φ += 2π * Float64((true_doppler - freq_update) * integration_time)
    end
    @test state.frequency_lock.locked
    @test max_error < 0.2
    @test abs(rem2pi(φ, RoundNearest)) < 0.01
    # Pure PLL from lock on: the FLL branch adds nothing to the frequency state.
    @test freq_update ≈ true_doppler atol = 0.01Hz

    # Four-quadrant with the sync's sign, Costas without one.
    pll(polarity, prompt) =
        Tracking.pll_disc(GPSL5Q(), _prompt_correlator(prompt); polarity = Int8(polarity))
    @test pll(0, -cis(2.5)) ≈ (2.5 - π) / 2π
    @test pll(-1, -cis(2.5)) ≈ 2.5 / 2π
    @test pll(1, -cis(0.1)) ≈ (0.1 - π) / 2π
end

@testset "Vector closure: FLL-assisted, four-quadrant from the sync" begin
    vt_state = SatVectorPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        vt_on = true,
    )
    # A 120° advance in 1 ms: 333 Hz four-quadrant, -167 Hz two-quadrant.
    _, wiped_off = close_loops(vt_state, cis(2π / 3), cis(0.0), true, -1)
    _, with_flips = close_loops(vt_state, cis(2π / 3), cis(0.0), false, 0)
    @test Tracking.mean_carrier_discr(wiped_off) ≈ (333 + 1 / 3) * 1Hz
    @test Tracking.mean_carrier_discr(with_flips) ≈ -(166 + 2 / 3) * 1Hz
    # No frequency lock indicator under vector closure: it stays FLL-assisted.
    for _ = 1:2000
        _, vt_state = close_loops(vt_state, cis(0.01), cis(0.0))
    end
    @test vt_state.frequency_lock === FrequencyLockIndicator()

    # With `vt_on` unset it stages like the scalar loop.
    fallback = SatVectorPLLAndDLL(vt_state; vt_on = false)
    for _ = 1:1500
        _, fallback = close_loops(fallback, cis(0.01), cis(0.0))
    end
    @test fallback.frequency_lock.locked
end

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
        if isnothing(switched_at) &&
           Tracking._sync_polarity(signal, driver.bit_buffer, prn) != 0
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

end
