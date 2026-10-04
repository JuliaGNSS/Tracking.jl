module SignalCombiningTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz, ns, s, ustrip
using StaticArrays: SVector
using Random: Xoshiro
using Statistics: mean, std
using GNSSSignals:
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GPSL5I,
    GPSL5Q,
    GalileoE1B,
    GalileoE1C,
    gen_code,
    get_center_frequency,
    get_code_frequency,
    get_code_length,
    get_secondary_code_length
using TrackingLoopFilters:
    SecondOrderBilinearLF, ThirdOrderAssistedBilinearLF, ThirdOrderBilinearLF, filter_loop
import Tracking
using Tracking:
    CorrelatorOutput,
    SignalCombiningSums,
    SatConventionalPLLAndDLL,
    SatVectorPLLAndDLL,
    SignalGroup,
    TrackState,
    TrackedSignal,
    ConventionalAssistedPLLAndDLL,
    VectorPLLAndDLL,
    add_satellite!,
    enable_vt!,
    mean_carrier_discr,
    mean_code_discr,
    reset_carrier_discr_accs!,
    reset_code_discr_accs!,
    append_correlator_output!,
    dll_disc,
    estimate_dopplers_and_filter_prompt!,
    get_code_phase,
    get_default_correlator,
    get_doppler_estimator_state,
    get_group_delay,
    get_sat_state,
    pll_disc,
    set_group_delay!,
    track!,
    update_accumulator

# A VEML correlator (Galileo E1's default) with the taps around prompt `p`.
_veml(p; early = 0.7, late = 0.7) = update_accumulator(
    get_default_correlator(GalileoE1B()),
    p .* SVector(0.2, early, 1.0, late, 0.2),
)

_e1_state(; combining = true) = TrackState(;
    signals = (e1 = (GalileoE1C(), GalileoE1B()),),
    doppler_estimator = ConventionalAssistedPLLAndDLL(; combine_signals = combining),
)

# Mark a passenger synced up front. A signal with one symbol per code period
# (Galileo E1B, GPS L1C-D) syncs on its first record, which would make that chunk
# the sync phase snap that drops pending passenger sums.
function _mark_synced!(track_state, group, prn, index)
    sats = Tracking.get_sat_states(track_state, group)
    signal = sats[prn].signals[index]
    synced = Tracking.BitBuffer(signal.bit_buffer.code_block_buffer, 0, true, 0.0im, 0)
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

@testset "Signal combining is an estimator setting, off by default" begin
    @test !Tracking.ConventionalPLLAndDLL().combine_signals
    @test !ConventionalAssistedPLLAndDLL().combine_signals
    @test !Tracking.VectorPLLAndDLL().combine_signals
    estimator = ConventionalAssistedPLLAndDLL(; combine_signals = true)
    @test estimator.combine_signals
    @test Tracking.VectorPLLAndDLL(; combine_signals = true).combine_signals
    # The keyword-update constructor keeps it unless told otherwise.
    @test Tracking.ConventionalPLLAndDLL(estimator; code_loop_filter_bandwidth = 2.0Hz).combine_signals
    @test !Tracking.ConventionalPLLAndDLL(estimator; combine_signals = false).combine_signals
    @test _e1_state().doppler_estimator.combine_signals
    @test !TrackState(; signal = GPSL1CA()).doppler_estimator.combine_signals
end

@testset "Weights and the weighted mean" begin
    T = 1 / 250Hz  # 4 ms
    @test Tracking._discriminator_weight(GalileoE1B(), T) ≈ 0.5 * 0.004s
    @test Tracking._fll_discriminator_weight(GalileoE1B(), T) ≈ 0.5 * (0.004s)^3
    @test Tracking._discriminator_weight(GPSL1CA(), T) ≈ 0.7079457843841379 * 0.004s
    # Without passengers the driver's own reading, bit for bit.
    @test Tracking._weighted_mean(0.1234567, 0.002s, Tracking.WeightedSum(0.0s, 0.0s)) ===
          0.1234567
    @test Tracking._weighted_mean(0.1, 0.002s, Tracking.WeightedSum(0.002s * 0.3, 0.002s)) ≈
          0.2
    @test Tracking._weighted_mean(0.1, 0.002s, Tracking.WeightedSum(0.006s * 0.3, 0.006s)) ≈
          0.25
    # The FLL's mean stays in Hz.
    @test Tracking._weighted_mean(
        1.0Hz,
        1.0s^3,
        Tracking.WeightedSum(3.0Hz * 1.0s^3, 1.0s^3),
    ) === 2.0Hz
end

@testset "A passenger record's contribution" begin
    T = 1 / 1000Hz
    fs = 20e6Hz
    code_doppler = 0.0Hz
    driver = TrackedSignal(GPSL5Q(); group_delay = 1.0ns)
    passenger = TrackedSignal(GPSL5I(); group_delay = 11.0ns)
    code_frequency = get_code_frequency(GPSL5Q())
    context = Tracking._passenger_context(
        passenger,
        (nothing, true),
        driver,
        Tracking.get_carrier_phase_offset(GPSL5Q()),
        code_frequency,
    )
    # 10 ns at 10.23 Mchip/s: the later passenger reads the driver's error minus
    # 0.1023 chip, which the offset adds back.
    @test context.differential_group_delay_chips ≈ 0.1023
    all_loops = (pll = true, fll = true, dll = true)
    # L5I sits a quarter cycle from L5Q, which the loops lock onto the real axis:
    # its PLL is read on the driver's frame.
    correlator = update_accumulator(
        get_default_correlator(GPSL5I()),
        cis(π / 2 + 0.1) .* SVector(0.5, 1.0, 0.4),
    )
    sums = Tracking._add_passenger_discriminators(
        SignalCombiningSums(),
        GPSL5I(),
        correlator,
        cis(π / 2),
        T,
        context,
        all_loops,
        code_doppler,
        fs,
    )
    weight = 0.5 * 0.001s
    @test sums.pll.weight ≈ weight
    @test sums.pll.sum / sums.pll.weight ≈ 0.1 / 2π
    @test sums.fll.weight ≈ weight * (0.001s)^2
    @test sums.fll.sum / sums.fll.weight ≈ 0.1 / (2π * 0.001) * Hz
    @test sums.dll.weight ≈ weight
    @test sums.dll.sum / sums.dll.weight ≈
          dll_disc(GPSL5I(), correlator, code_doppler, fs) + 0.1023

    # Only the loops it is combined into, and no FLL without a previous prompt or a DLL
    # without both group delays.
    none = Tracking._add_passenger_discriminators(
        SignalCombiningSums(),
        GPSL5I(),
        correlator,
        cis(π / 2),
        T,
        context,
        (pll = false, fll = false, dll = false),
        code_doppler,
        fs,
    )
    @test none === SignalCombiningSums()
    unknown = Tracking._passenger_context(
        TrackedSignal(GPSL5I()),
        (nothing, true),
        driver,
        0.0,
        code_frequency,
    )
    @test isnan(unknown.differential_group_delay_chips)
    partial = Tracking._add_passenger_discriminators(
        SignalCombiningSums(),
        GPSL5I(),
        correlator,
        0.0im,
        T,
        unknown,
        all_loops,
        code_doppler,
        fs,
    )
    @test iszero(partial.fll.weight) &&
          iszero(partial.dll.weight) &&
          partial.pll.weight > 0s
end

# Close one record's loops; returns the carrier and code update and the state.
function close(
    state,
    combining;
    wiped_off = false,
    polarity = 0,
    previous_prompt = cis(0.0),
)
    Tracking._close_loops(
        state,
        GalileoE1B(),
        _veml(cis(0.2); early = 0.8, late = 0.6),
        previous_prompt,
        0.0Hz,
        16.368e6Hz,
        1 / 250Hz,
        18.0Hz,
        1.0Hz,
        wiped_off,
        Int8(polarity),
        combining,
    )
end

include("frequency_lock_test_helpers.jl")  # `_locked`

@testset "Which loops the passengers are combined into" begin
    loops_to_combine = Tracking._loops_to_combine
    assisted = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter = ThirdOrderAssistedBilinearLF(),
    )
    @test loops_to_combine(assisted, true) == (pll = true, fll = true, dll = true)
    # Not into an FLL that is not formed: after frequency lock, or without an
    # FLL-assisted carrier filter.
    locked = SatConventionalPLLAndDLL(assisted; frequency_lock = _locked())
    @test !loops_to_combine(locked, true).fll
    plain = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter = ThirdOrderBilinearLF(),
    )
    @test !loops_to_combine(plain, true).fll
    # Without signal combining into none.
    @test loops_to_combine(assisted, false) == (pll = false, fll = false, dll = false)
    # Under vector closure only into the PLL.
    vt = SatVectorPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        code_discr_accs = ((0, 0.0),),
        carrier_discr_accs = ((0, 0.0Hz),),
        vt_on = true,
    )
    @test loops_to_combine(vt, true) == (pll = true, fll = false, dll = false)
end

@testset "The loop closure mixes the passengers in" begin
    state = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        carrier_loop_filter = ThirdOrderAssistedBilinearLF(),
    )
    w = 0.5 * 0.004s
    # No passenger weight: exactly the driver alone.
    @test close(state, SignalCombiningSums()) === Tracking._close_loops(
        state,
        GalileoE1B(),
        _veml(cis(0.2); early = 0.8, late = 0.6),
        cis(0.0),
        0.0Hz,
        16.368e6Hz,
        1 / 250Hz,
        18.0Hz,
        1.0Hz,
        false,
        Int8(0),
    )
    # The DLL reads the weighted mean of the driver's and the passengers'.
    own_dll =
        dll_disc(GalileoE1B(), _veml(cis(0.2); early = 0.8, late = 0.6), 0.0Hz, 16.368e6Hz)
    _, code_update, _ = close(
        state,
        SignalCombiningSums(
            Tracking.WeightedSum(0s, 0s),
            Tracking.WeightedSum(0Hz * 0s^3, 0s^3),
            Tracking.WeightedSum(3w * 0.05, 3w),
        ),
    )
    expected = first(
        filter_loop(SecondOrderBilinearLF(), (own_dll + 3 * 0.05) / 4, 1 / 250Hz, 1.0Hz),
    )
    @test code_update ≈ expected

    pll = SignalCombiningSums(
        Tracking.WeightedSum(w * 0.01, w),
        Tracking.WeightedSum(0Hz * 0s^3, 0s^3),
        Tracking.WeightedSum(0s, 0s),
    )
    fll = SignalCombiningSums(
        Tracking.WeightedSum(0s, 0s),
        Tracking.WeightedSum(w * (0.004s)^2 * 7.0Hz, w * (0.004s)^2),
        Tracking.WeightedSum(0s, 0s),
    )
    alone = first(close(state, SignalCombiningSums()))
    # The closure combines the sums it is handed (the fold decides which loops
    # get any, see "Which loops the passengers are combined into"); an FLL that
    # is no longer formed takes none.
    @test first(close(state, pll)) != alone
    @test first(close(state, fll)) != alone
    # Into a four-quadrant driver too while its reading lies within the
    # passengers' two-quadrant range, and not beyond it.
    @test first(close(state, pll; polarity = 1)) !=
          first(close(state, SignalCombiningSums(); polarity = 1))
    @test Tracking._gated_mean(0.2, 1.0s, Tracking.WeightedSum(0.1s, 1.0s), 0.25) ≈ 0.15
    @test Tracking._gated_mean(0.3, 1.0s, Tracking.WeightedSum(0.1s, 1.0s), 0.25) === 0.3
    locked = SatConventionalPLLAndDLL(state; frequency_lock = _locked())
    @test first(close(locked, fll)) == first(close(locked, SignalCombiningSums()))

    # Under vector closure the passengers are combined into the PLL only, and the
    # navigation filter gets the driver's own DLL and FLL readings.
    vt = SatVectorPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        code_discr_accs = ((0, 0.0),),
        carrier_discr_accs = ((0, 0.0Hz),),
        vt_on = true,
    )
    @test Tracking._locally_closed_loops(vt) == (pll = true, fll = false, dll = false)
    @test Tracking._locally_closed_loops(SatVectorPLLAndDLL(vt; vt_on = false)) ==
          (pll = true, fll = true, dll = true)
    mixed = SignalCombiningSums(
        Tracking.WeightedSum(w * 0.01, w),
        Tracking.WeightedSum(w * (0.004s)^2 * 7.0Hz, w * (0.004s)^2),
        Tracking.WeightedSum(w * 0.05, w),
    )
    combined_carrier, _, combined_vt = close(vt, mixed)
    alone_carrier, _, alone_vt = close(vt, SignalCombiningSums())
    @test combined_carrier != alone_carrier
    @test combined_vt.code_discr_accs == alone_vt.code_discr_accs
    @test combined_vt.carrier_discr_accs == alone_vt.carrier_discr_accs
    @test first(
        close(
            vt,
            SignalCombiningSums(
                Tracking.WeightedSum(0s, 0s),
                Tracking.WeightedSum(w * (0.004s)^2 * 7.0Hz, w * (0.004s)^2),
                Tracking.WeightedSum(w * 0.05, w),
            ),
        ),
    ) == alone_carrier
end

@testset "Several passengers add up" begin
    # GPS L1C-P drives; L1C-D (10 ms) and L1 C/A (1 ms) are its passengers.
    fs = 5e6Hz
    taps(c) =
        length(Tracking.get_accumulators(c)) == 3 ? SVector(0.5, 1.0, 0.5) :
        SVector(0.2, 0.7, 1.0, 0.7, 0.2)
    function record(signal, p, n, sample_index)
        c = get_default_correlator(signal)
        CorrelatorOutput(update_accumulator(c, p * n .* taps(c)), n, sample_index)
    end
    function pending_sums(records)
        track_state = TrackState(;
            signals = (l1 = (GPSL1C_P(), GPSL1C_D(), GPSL1CA()),),
            doppler_estimator = ConventionalAssistedPLLAndDLL(; combine_signals = true),
        )
        add_satellite!(
            track_state;
            group = :l1,
            prn = 1,
            code_phase = 0.0,
            carrier_doppler = 0.0Hz,
        )
        _mark_synced!(track_state, :l1, 1, 2)
        for (signal, output) in records
            append_correlator_output!(track_state, output, :l1, 1, signal)
        end
        # No driver record: every passenger record stays pending.
        estimate_dopplers_and_filter_prompt!(track_state, (L1 = fs,))
        get_doppler_estimator_state(get_sat_state(track_state, :l1, 1)).signal_combining_sums
    end
    d_records = [(GPSL1C_D, record(GPSL1C_D(), cis(0.1), 50000, 50000))]
    ca_records = [
        (GPSL1CA, record(GPSL1CA(), cis(0.2), 5000, 5000)),
        (GPSL1CA, record(GPSL1CA(), cis(0.25), 5000, 10000)),
    ]
    d_only = pending_sums(d_records)
    ca_only = pending_sums(ca_records)
    both = pending_sums([d_records; ca_records])
    @test d_only.pll.weight ≈ Tracking._discriminator_weight(GPSL1C_D(), 50000 / fs)
    @test ca_only.pll.weight ≈ 2 * Tracking._discriminator_weight(GPSL1CA(), 5000 / fs)
    # Each passenger adds its own records, whatever the others do.
    @test both.pll.weight ≈ d_only.pll.weight + ca_only.pll.weight
    @test both.pll.sum ≈ d_only.pll.sum + ca_only.pll.sum
    # Only L1 C/A's second record has a previous prompt for the FLL.
    @test iszero(d_only.fll.weight)
    @test both.fll.weight ≈ ca_only.fll.weight > 0s^3
    @test both.fll.sum ≈ ca_only.fll.sum
end

@testset "Passenger records are combined into the driver record they end within" begin
    fs = 16.368e6Hz
    n = 65472  # 4 ms
    track_state = _e1_state()
    add_satellite!(
        track_state;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    _mark_synced!(track_state, :e1, 11, 2)
    record(p, sample_index) = CorrelatorOutput(_veml(p * n), n, sample_index)
    sums() =
        get_doppler_estimator_state(get_sat_state(track_state, :e1, 11)).signal_combining_sums
    estimate!() = estimate_dopplers_and_filter_prompt!(track_state, (L1 = fs,))

    # A chunk with a passenger record only: it stays pending. Group delays are
    # unknown, so it is combined into the carrier loops only (the FLL has no
    # previous prompt).
    append_correlator_output!(track_state, record(cis(0.1), n), :e1, 11, GalileoE1B)
    estimate!()
    @test sums().pll.weight ≈ 0.5 * 0.004s
    @test sums().pll.sum / sums().pll.weight ≈ 0.1 / 2π
    @test iszero(sums().dll.weight)
    @test isempty(Tracking.get_correlator_outputs(track_state, :e1, 11, GalileoE1B))

    # The driver's next record consumes it with the passenger record ending by
    # then; one ending later waits for the following driver record.
    append_correlator_output!(track_state, record(cis(0.1), n), :e1, 11, GalileoE1B)
    append_correlator_output!(track_state, record(cis(0.0), 2n), :e1, 11, GalileoE1C)
    append_correlator_output!(track_state, record(cis(0.3), 3n), :e1, 11, GalileoE1B)
    estimate!()
    @test sums().pll.weight ≈ 0.5 * 0.004s
    @test sums().pll.sum / sums().pll.weight ≈ 0.3 / 2π
    @test sums().fll.weight > 0s^3  # this one had a previous prompt

    # Without combining nothing is pending.
    uncombined = _e1_state(; combining = false)
    add_satellite!(
        uncombined;
        group = :e1,
        prn = 11,
        code_phase = 0.0,
        carrier_doppler = 0.0Hz,
    )
    append_correlator_output!(uncombined, record(cis(0.1), n), :e1, 11, GalileoE1B)
    estimate_dopplers_and_filter_prompt!(uncombined, (L1 = fs,))
    @test get_doppler_estimator_state(get_sat_state(uncombined, :e1, 11)).signal_combining_sums ===
          SignalCombiningSums()
end

@testset "Vector tracking hands the navigation filter every signal's readings" begin
    fs = 16.368e6Hz
    n = 65472
    record(p, sample_index; early = 0.7, late = 0.7) =
        CorrelatorOutput(_veml(p * n; early, late), n, sample_index)
    function vt_state(; combining)
        track_state = TrackState(;
            signals = (e1 = (GalileoE1C(), GalileoE1B()),),
            doppler_estimator = VectorPLLAndDLL(; combine_signals = combining),
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
    function step!(track_state)
        append_correlator_output!(track_state, record(cis(0.0), n), :e1, 11, GalileoE1C)
        append_correlator_output!(
            track_state,
            record(cis(0.1), n ÷ 2; early = 0.8, late = 0.6),
            :e1,
            11,
            GalileoE1B,
        )
        append_correlator_output!(
            track_state,
            record(cis(0.1), n; early = 0.8, late = 0.6),
            :e1,
            11,
            GalileoE1B,
        )
        estimate_dopplers_and_filter_prompt!(track_state, (L1 = fs,))
    end
    state(track_state) = get_doppler_estimator_state(get_sat_state(track_state, :e1, 11))

    for combining in (false, true)
        track_state = vt_state(; combining)
        # Outside the vector loop nothing is accumulated.
        step!(track_state)
        @test state(track_state).code_discr_accs == ((0, 0.0), (0, 0.0))
        enable_vt!(track_state, :e1, (11,))
        step!(track_state)
        accs = state(track_state).code_discr_accs
        # One reading per record: the driver's, and the passenger's two.
        @test first.(accs) == (1, 2)
        @test first.(state(track_state).carrier_discr_accs) == (1, 2)
        passenger_dll =
            dll_disc(GalileoE1B(), _veml(cis(0.1); early = 0.8, late = 0.6), 0.0Hz, fs)
        @test abs(passenger_dll) > 0.01
        @test mean_code_discr(track_state, :e1, 11, GalileoE1B) ≈ passenger_dll
        @test mean_code_discr(get_sat_state(track_state, :e1, 11), 2) ≈ passenger_dll
        @test mean_carrier_discr(track_state, :e1, 11, 2) isa typeof(1.0Hz)
        @test_throws ArgumentError mean_code_discr(state(track_state))
        reset_code_discr_accs!(track_state)
        reset_carrier_discr_accs!(track_state)
        @test state(track_state).code_discr_accs == ((0, 0.0), (0, 0.0))
        @test state(track_state).carrier_discr_accs == ((0, 0.0Hz), (0, 0.0Hz))
    end

    # A record without a previous prompt has no FLL reading and is not counted:
    # here the driver's first and the passenger's first.
    track_state = vt_state(; combining = false)
    enable_vt!(track_state, :e1, (11,))
    step!(track_state)
    @test first.(state(track_state).code_discr_accs) == (1, 2)
    @test first.(state(track_state).carrier_discr_accs) == (0, 1)
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

@testset "The sync phase snap drops pending passenger sums" begin
    state = SatConventionalPLLAndDLL(;
        init_carrier_doppler = 0.0Hz,
        init_code_doppler = 0.0Hz,
        signal_combining_sums = SignalCombiningSums(
            Tracking.WeightedSum(0.001s, 0.002s),
            Tracking.WeightedSum(0Hz * 0s^3, 0s^3),
            Tracking.WeightedSum(0s, 0s),
        ),
    )
    @test Tracking._drop_pending_combining(state, true).signal_combining_sums ===
          SignalCombiningSums()
    # Without signal combining nothing is pending, and the state is kept as is.
    @test Tracking._drop_pending_combining(state, false) === state
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
    fs = 16.368e6Hz
    n = 65472
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
    driver_record = CorrelatorOutput(_veml(cis(0.1) * n), n, n)
    passenger_record = CorrelatorOutput(_veml(cis(0.2) * n), n, n)
    sampling_frequencies = (L1 = fs,)
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
