module MultiSignalCombiningTest

# Tests for multi-signal discriminator combining: a passenger signal's
# (`signals[2:end]`) record that ends on the same sample as a driver record is
# folded into that driver record's loop update as a weighted mean.
#
# The properties worth pinning, in order of how badly a regression would hurt:
#
#  1. **Bit-identity where nothing is combined.** A single-signal satellite, a
#     multi-signal one with `discriminator_combining = false`, and a driver record
#     no passenger record coincides with must all produce exactly the Dopplers
#     the driver alone produces. The whole feature is safe-by-default only
#     because of this.
#  2. **Coincident records only, and nothing carried.** A passenger record that
#     does not end on a driver record's sample is still applied to its own
#     signal, but its discriminators reach no loop — not this update, not a later
#     one.
#  3. **The weighted mean is a mean.** Two signals carrying the same information
#     must not double the loop gain — a plain sum would.
#  4. **The inter-signal correction is subtracted with the right sign**, and is
#     converted to chips at the code frequency in effect.
#  5. **An unset group delay gates the code loop only** — the carrier loops
#     combine from the first integration, with nothing supplied by the caller.

using Test: @allocated, @inferred, @test, @test_throws, @testset
using Unitful: Hz, NoUnits, s, uconvert
using GNSSSignals:
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GalileoE1C,
    GPSL5I,
    GPSL5Q,
    gen_code,
    get_carrier_phase_offset,
    get_code_frequency,
    get_code_center_frequency_ratio
using Tracking:
    Tracking,
    ConventionalAssistedPLLAndDLL,
    CorrelatorOutput,
    DiscriminatorAccumulator,
    EarlyPromptLateCorrelator,
    NoCN0Estimator,
    NumAnts,
    TrackState,
    TrackedSat,
    TrackedSignal,
    VeryEarlyPromptLateCorrelator,
    add_satellite!,
    append_correlator_output!,
    dll_disc_noise_gain,
    estimate_dopplers_and_filter_prompt!,
    get_carrier_doppler,
    get_code_doppler,
    get_group_delay,
    get_preferred_num_code_blocks_to_integrate,
    get_sat_state,
    set_group_delay!,
    set_preferred_num_code_blocks_to_integrate!,
    set_loop_filter_bandwidths!,
    track
using Tracking: RecordFoldContext, _normalized_mean, _relative_record_snr
using Tracking:
    VectorPLLAndDLL,
    enable_vt!,
    mean_carrier_discr,
    mean_code_discr,
    reset_carrier_discr_acc!,
    reset_code_discr_acc!

const FS = 5e6Hz

# One millisecond of a clean GPS L1 C/A signal at the given Doppler / code phase
# — the same fixture shape the other track tests use.
function l1ca_signal(prn, carrier_doppler, start_code_phase, num_samples, fs = FS)
    gpsl1 = GPSL1CA()
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + get_code_frequency(gpsl1)
    range_ = 0:(num_samples-1)
    cis.(2π .* carrier_doppler .* range_ ./ fs) .*
    gen_code(num_samples, gpsl1, prn, fs, code_frequency, start_code_phase)
end

# A satellite tracking `n` identical GPS L1 C/A signals. Identical signals make
# the combining arithmetic checkable in closed form: every signal sees the same
# samples, so every discriminator is the same number and any *mean* of them is
# that number, while a sum would be `n` times it.
function identical_l1ca_track_state(
    n,
    prn,
    carrier_doppler;
    discriminator_combining = true,
    kwargs...,
)
    gpsl1 = GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL(; kwargs...)
    tracked_signals = ntuple(
        _ -> TrackedSignal(
            gpsl1;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            # Two copies of one signal have no group delay difference between
            # them by construction, so pin 0.0 and the code loop combines. The
            # `nothing` default — and what it gates — is exercised separately
            # below.
            group_delay = 0.0s,
        ),
        n,
    )
    sat = TrackedSat(
        tracked_signals,
        prn,
        0.0,
        carrier_doppler;
        doppler_estimator = estimator,
    )
    TrackState(gpsl1, sat; doppler_estimator = estimator, discriminator_combining)
end

@testset "single-signal sats are bit-identical with combining on" begin
    # A satellite with no passengers has nothing to combine, so with combining
    # switched on it must track *exactly* as with it off — not merely to within a
    # tolerance. `_seed_and_fold` dispatches on the empty passenger tuple to
    # guarantee this at the type level.
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 20_000)

    ts_on = TrackState(;
        signal = GPSL1CA(),
        doppler_estimator = ConventionalAssistedPLLAndDLL(),
        discriminator_combining = true,
    )
    ts_on = add_satellite!(ts_on; prn, code_phase = 0.0, carrier_doppler)
    ts_off = TrackState(;
        signal = GPSL1CA(),
        doppler_estimator = ConventionalAssistedPLLAndDLL(),
    )
    ts_off = add_satellite!(ts_off; prn, code_phase = 0.0, carrier_doppler)

    on = track(copy(buf), ts_on, FS)
    off = track(copy(buf), ts_off, FS)

    @test get_carrier_doppler(get_sat_state(on, prn)) ===
          get_carrier_doppler(get_sat_state(off, prn))
    @test get_code_doppler(get_sat_state(on, prn)) ===
          get_code_doppler(get_sat_state(off, prn))
    @test Tracking.get_code_phase(get_sat_state(on, prn)) ===
          Tracking.get_code_phase(get_sat_state(off, prn))
end

@testset "combining two identical signals is a mean, not a sum" begin
    # Both signals correlate the same samples with the same replica, so both
    # discriminators equal the driver's and their records end on the same
    # sample. A weight-normalized mean of equal values is that value, so the
    # Dopplers must match the driver-only run closely; a summing combiner would
    # double the loop gain and visibly diverge.
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 40_000)

    combined = track(copy(buf), identical_l1ca_track_state(2, prn, carrier_doppler), FS)
    driver_only = track(
        copy(buf),
        identical_l1ca_track_state(
            2,
            prn,
            carrier_doppler;
            discriminator_combining = false,
        ),
        FS,
    )

    c_carrier = get_carrier_doppler(get_sat_state(combined, prn))
    d_carrier = get_carrier_doppler(get_sat_state(driver_only, prn))
    c_code = get_code_doppler(get_sat_state(combined, prn))
    d_code = get_code_doppler(get_sat_state(driver_only, prn))

    # The mean of two equal discriminators is that discriminator, so this holds
    # to floating-point round-off rather than to a loop-convergence tolerance.
    @test c_carrier ≈ d_carrier rtol = 1e-10
    @test c_code ≈ d_code rtol = 1e-10

    # Both signals really did integrate (i.e. the passenger fold ran and its
    # records were consumed), so the agreement above is not the trivial
    # "passenger produced nothing" case.
    sat = get_sat_state(combined, prn)
    @test length(Tracking.get_filtered_prompts(sat.signals[1])) > 1
    @test length(Tracking.get_filtered_prompts(sat.signals[2])) ==
          length(Tracking.get_filtered_prompts(sat.signals[1]))
    @test isempty(Tracking.get_correlator_outputs(sat.signals[2]))
end

@testset "a longer driver leaves the short passenger's records applied and consumed" begin
    # A 10 ms L1C-P driver with a 1 ms C/A passenger: at the default 1 ms chunk
    # most chunks hold passenger records and no driver record at all. With no
    # state carried across chunks those records reach no loop — but every one of
    # them must still be applied to its signal and cleared from its buffer, or it
    # would be folded twice.
    gpsl1, l1cp = GPSL1CA(), GPSL1C_P()
    prn, carrier_doppler = 1, 1000.0Hz
    ts = TrackState(;
        signals = (gps_l1 = (l1cp, gpsl1),),
        doppler_estimator = ConventionalAssistedPLLAndDLL(),
        discriminator_combining = true,
    )
    ts = add_satellite!(ts; group = :gps_l1, prn, code_phase = 0.0, carrier_doppler)

    # L1C-P is TMBOC(6,1,4/33), so its replica needs ≥ 12.276 MHz — this testset
    # runs at a higher rate than the rest of the file for that reason alone.
    fs = 15e6Hz
    ts = track(l1ca_signal(prn, carrier_doppler, 0.0, 75_000, fs), ts, fs)   # 5 ms
    sat = get_sat_state(ts, prn)

    # No driver (L1C-P) integration completed in 5 ms, but five C/A ones did.
    @test isempty(Tracking.get_correlator_outputs(sat.signals[1]))
    @test isempty(Tracking.get_correlator_outputs(sat.signals[2]))
    @test length(Tracking.get_filtered_prompts(sat.signals[2])) == 5
    @test isempty(Tracking.get_filtered_prompts(sat.signals[1]))
    # …and with no driver record the Dopplers held.
    @test get_carrier_doppler(sat) === carrier_doppler
    # Nothing is parked anywhere: combining keeps no per-satellite state.
    @test !hasfield(Tracking.SatConventionalPLLAndDLL, :pending_combining_sums)
end

# --------------------------------------------------------------------------
# One record's contribution
# --------------------------------------------------------------------------

# Run one passenger record through the passenger walk against a driver record
# ending at `driver_sample_index`, and return what it contributed to that driver
# record's loop update as a `DiscriminatorAccumulator`. Going through
# `_advance_one_passenger_to` rather than calling the weighting directly pins the
# wiring too: the coincidence test, the de-rotation and the group-delay referral
# the context carries.
function passenger_contribution(
    tracked_signal,
    correlator,
    integrated_samples,
    passenger_delay;
    code_doppler = 0.0Hz,
    # Default: driver and passenger co-phased, so the PLL de-rotation is the
    # identity and the accumulated phase equals the bare discriminator. Not
    # `0.0` — GPS L1 C/A's own carrier phase offset is −π/2, so a literal zero
    # here would inject a 90° rotation and mean "the driver is an L1C component".
    driver_carrier_phase_offset = get_carrier_phase_offset(
        Tracking.get_signal(tracked_signal),
    ),
    previous_prompt = complex(0.0, 0.0),
    # The driver is the datum unless a test says otherwise, so the passenger's
    # stored value is directly its referred bias.
    driver_delay = 0.0s,
    sample_index = integrated_samples,
    driver_sample_index = integrated_samples,
)
    signal = Tracking.get_signal(tracked_signal)
    # Two steps: `correlator_outputs` and the previous prompt are kwargs of the
    # copy-update constructor only. `NoCN0Estimator` keeps the noise-density
    # plumbing out of it — no noise estimator here, hence a `nothing` floor.
    ts = TrackedSignal(
        TrackedSignal(
            signal;
            correlator,
            group_delay = passenger_delay,
            cn0_estimator = NoCN0Estimator(),
        );
        last_fully_integrated_filtered_prompt = previous_prompt,
        correlator_outputs = [
            CorrelatorOutput(correlator, integrated_samples, sample_index),
        ],
    )
    ctx = RecordFoldContext(
        ts,
        1,
        2,
        FS,
        code_doppler,
        Float64(driver_carrier_phase_offset),
        driver_delay,
        (nothing, true),
    )
    new_ts, cursor, contribution, _ =
        Tracking._advance_one_passenger_to(ts, 1, ctx, driver_sample_index, nothing)
    # Coincident or not, the record is consumed from the walk's point of view
    # whenever it ended at or before the driver record.
    @test cursor == (sample_index <= driver_sample_index ? 2 : 1)
    contribution
end

@testset "the group delay is subtracted, in chips, with the right sign" begin
    gpsl1 = GPSL1CA()
    # An early-late-asymmetric correlator: a real code error, so the raw
    # discriminator is non-zero and the correction is visible against it.
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    base = TrackedSignal(gpsl1; correlator)
    raw = Tracking.dll_disc(gpsl1, correlator, 0.0Hz, FS)
    gain = dll_disc_noise_gain(gpsl1, correlator, 0.0Hz, FS)
    snr = _relative_record_snr(gpsl1, n)

    # A zero delay asserts "this signal's code phase *is* the driver's", so the
    # discriminator enters as measured, weighted by SNR over the noise gain.
    same = passenger_contribution(base, correlator, n, 0.0s)
    @test same.dll_weight ≈ snr / gain
    @test same.dll_sum / same.dll_weight ≈ raw

    # A *larger* group delay than the driver's means the passenger arrives later
    # and so sits at a **smaller** code phase: its raw discriminator under-reads
    # the shared error by `δ · f_code` chips, which the referral adds back.
    δt = 4.0e-9s
    # Computed independently of the implementation, from seconds × chips/second.
    δchips = uconvert(NoUnits, δt * get_code_frequency(gpsl1))
    later = passenger_contribution(base, correlator, n, δt)
    @test later.dll_sum / later.dll_weight ≈ raw + δchips
    # …and the opposite sign for a passenger that leads the driver.
    earlier = passenger_contribution(base, correlator, n, -δt)
    @test earlier.dll_sum / earlier.dll_weight ≈ raw - δchips
    # Only the difference is read: shifting both ends by one constant is a no-op.
    shifted =
        passenger_contribution(base, correlator, n, δt + 7.0e-9s; driver_delay = 7.0e-9s)
    @test shifted.dll_sum / shifted.dll_weight ≈ raw + δchips

    # Converted at the code frequency actually in effect, code Doppler included —
    # not at the nominal chip rate.
    code_doppler = 1.0Hz * 4000
    doppler_shifted = passenger_contribution(base, correlator, n, δt; code_doppler)
    raw_shifted = Tracking.dll_disc(gpsl1, correlator, code_doppler, FS)
    @test doppler_shifted.dll_sum / doppler_shifted.dll_weight ≈
          raw_shifted + uconvert(NoUnits, δt * (get_code_frequency(gpsl1) + code_doppler))
end

@testset "an unset group delay gates the code loop only" begin
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    base = TrackedSignal(gpsl1; correlator)
    raw = Tracking.dll_disc(gpsl1, correlator, 0.0Hz, FS)

    # No delay supplied: the passenger stays out of the code loop …
    unset = passenger_contribution(base, correlator, n, nothing)
    @test unset.dll_weight == 0.0
    # … but its carrier phase contribution is unaffected: carrier phase carries
    # no inter-signal group delay, so it combines from the first integration with
    # nothing supplied at all.
    @test unset.pll_weight > 0
    @test unset.pll_sum / unset.pll_weight ≈ Tracking.pll_disc(gpsl1, correlator)
    # The same when it is the *driver's* datum that is unknown.
    no_datum = passenger_contribution(base, correlator, n, 0.0s; driver_delay = nothing)
    @test no_datum.dll_weight == 0.0
    @test no_datum.pll_weight > 0

    # `0.0` is a different statement from `nothing`, and the only one of the two
    # that lets the code loop combine.
    zeroed = passenger_contribution(base, correlator, n, 0.0s)
    @test zeroed.dll_weight > 0
    @test zeroed.dll_sum / zeroed.dll_weight ≈ raw
end

@testset "only a record ending on the driver's sample contributes" begin
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    base = TrackedSignal(gpsl1; correlator)
    @test passenger_contribution(base, correlator, n, 0.0s).pll_weight > 0
    # One sample early: applied (the cursor advances), contributes nothing.
    @test passenger_contribution(base, correlator, n, 0.0s; sample_index = n - 1) ==
          DiscriminatorAccumulator()
    # One sample late: left for the next driver record, or for the drain.
    @test passenger_contribution(base, correlator, n, 0.0s; sample_index = n + 1) ==
          DiscriminatorAccumulator()
end

@testset "weights: SNR, discriminator noise gain, and the FLL's T²" begin
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    base = TrackedSignal(gpsl1; correlator)

    # Code-loop weight scales with the record length …
    short = passenger_contribution(base, correlator, 5000, 0.0s)
    long = passenger_contribution(base, correlator, 50_000, 0.0s)
    @test long.dll_weight ≈ 10 * short.dll_weight
    # … and so does the carrier phase loop's, which is the bare nominal SNR.
    @test short.pll_weight === _relative_record_snr(gpsl1, 5000)
    @test long.pll_weight ≈ 10 * short.pll_weight

    # A zeroed previous prompt means "no frequency measurement yet", not "zero
    # frequency error", so it must not pull the combined FLL towards 0 Hz.
    @test short.fll_weight == 0.0
    # With a previous prompt the FLL weight carries the extra `N²`: dividing the
    # inter-prompt phase difference by 2π·T scales its noise by 1/T.
    with_prev = passenger_contribution(
        base,
        correlator,
        5000,
        0.0s;
        previous_prompt = complex(0.9, 0.1),
    )
    with_prev_long = passenger_contribution(
        base,
        correlator,
        50_000,
        0.0s;
        previous_prompt = complex(0.9, 0.1),
    )
    @test with_prev.fll_weight ≈ _relative_record_snr(gpsl1, 5000) * 5000^2
    @test with_prev_long.fll_weight ≈ 1000 * with_prev.fll_weight
    @test with_prev.fll_sum / with_prev.fll_weight ≈
          Tracking.fll_disc(gpsl1, correlator, complex(0.9, 0.1), 5000 / FS)

    # A degenerate VEML tap layout carries no delay information: its `Inf` noise
    # gain gives the record weight *zero* in the code loop rather than a `NaN`
    # that would propagate into the driver's code loop filter.
    degenerate = VeryEarlyPromptLateCorrelator(
        [
            complex(0.60, 0.0),
            complex(0.90, 0.0),
            complex(1.00, 0.0),
            complex(0.80, 0.0),
            complex(0.50, 0.0),
        ],
        2.0,
        3.0,
    )
    excluded = passenger_contribution(
        TrackedSignal(GalileoE1C(); correlator = degenerate),
        degenerate,
        5000,
        0.0s,
    )
    @test excluded.dll_weight === 0.0
    @test _normalized_mean(excluded.dll_sum, excluded.dll_weight) === 0.0
    @test excluded.pll_weight > 0   # the carrier loops are unaffected
end

@testset "_normalized_mean divides, and reads zero weight as zero" begin
    @test _normalized_mean(4.0 * 2.0, 4.0) === 2.0
    # Equal values give that value regardless of the weights — the property that
    # keeps the loop gain fixed as contributors come and go.
    @test _normalized_mean(5.0 * 2.0 + 3.0 * 2.0, 5.0 + 3.0) ≈ 2.0
    # Zero weight answers a zero of the sum's own type, never `0/0`.
    @test _normalized_mean(0.0, 0.0) === 0.0
    @test _normalized_mean(0.0Hz, 0.0) === 0.0Hz

    # No passenger weight is the driver's own value *as is*, bit for bit.
    d = 0.1234567
    @test Tracking._combined_discriminator(d, 0.3, 0.0, 0.0) === d
    # With passenger weight, the weighted mean over all contributors.
    @test Tracking._combined_discriminator(1.0, 3.0, 1.0 * 2.0, 1.0) ≈ (3.0 + 2.0) / 4
    # A driver without an FLL measurement of its own (zero weight) takes the
    # passengers' mean.
    @test Tracking._combined_discriminator(0.0Hz, 0.0, 2.0 * 5.0Hz, 2.0) ≈ 5.0Hz
    # Zero total weight is a zero of the right type.
    @test Tracking._combined_discriminator(0.0Hz, 0.0, 0.0Hz, -0.0) === 0.0Hz
end

@testset "quadrature passengers are de-rotated before the PLL" begin
    # GPS L5 is a quadrature pair: with the pilot driving, the data component
    # sits on the imaginary axis. Without de-rotation `atan(imag/real)` reads
    # ±π/2 — a fabricated phase error that would drag the combined PLL off lock.
    l5i, l5q = GPSL5I(), GPSL5Q()
    driver_phase = get_carrier_phase_offset(l5q)
    on_axis = complex(1.0, 0.0)
    quadrature = on_axis * cis(get_carrier_phase_offset(l5i) - driver_phase)
    correlator =
        EarlyPromptLateCorrelator([complex(0.7, 0.0), quadrature, complex(0.7, 0.0)], 0.5)

    derotated = Tracking._derotate_correlator(
        correlator,
        Tracking._carrier_phase_derotation(driver_phase, l5i),
    )
    @test abs(Tracking.pll_disc(l5i, derotated)) < 1e-12
    @test abs(Tracking.pll_disc(l5i, correlator)) > 1.0     # the un-rotated trap
    # Every tap is rotated, and the magnitudes the code loop reads are untouched.
    @test Tracking.get_prompt(derotated) ≈ on_axis
    @test abs.(Tracking.get_accumulators(derotated)) ≈
          abs.(Tracking.get_accumulators(correlator))
    # The driver's own record is not rotated at all.
    @test Tracking._derotate_correlator(correlator, nothing) === correlator

    # …and the passenger walk really passes that rotation in. The delay is left
    # unset so the code loop abstains and only the carrier claim is under test.
    folded = passenger_contribution(
        TrackedSignal(l5i; correlator),
        correlator,
        5000,
        nothing;
        driver_carrier_phase_offset = driver_phase,
    )
    @test folded.pll_weight > 0
    @test abs(folded.pll_sum / folded.pll_weight) < 1e-12
end

@testset "only a pre-sync-invalidated passenger record is excluded" begin
    # Nominal weights say what a component's power share should be; they cannot
    # tell whether it is being received, so a record contributes on its nominal
    # weight regardless of how strong its prompt turned out to be.
    gpsl1 = GPSL1CA()
    strong = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    weak = EarlyPromptLateCorrelator(
        [complex(0.009, 0.0), complex(0.01, 0.0), complex(0.006, 0.0)],
        0.5,
    )
    base = TrackedSignal(gpsl1; correlator = strong)
    n = 5000
    a = passenger_contribution(base, strong, n, 0.0s)
    b = passenger_contribution(base, weak, n, 0.0s)
    @test a.pll_weight === b.pll_weight === _relative_record_snr(gpsl1, n)
    @test a.dll_weight ≈ b.dll_weight

    # The one exclusion: a record correlated with a pre-sync replica on a
    # secondary-coded signal, whose coherent sum partially cancels without the
    # secondary-code wipe-off, so its discriminators measure nothing.
    @test Tracking._record_invalidated_by_sync(GalileoE1C(), true)
    @test !Tracking._record_invalidated_by_sync(GalileoE1C(), false)
    # A signal with no secondary code correlates identically either side of the
    # sync instant, so it is unaffected.
    @test !Tracking._record_invalidated_by_sync(gpsl1, true)
end

@testset "the combining fold is inferable and adds no allocations" begin
    # The interleaved fold threads a tuple of signals, a tuple of cursors and a
    # tuple of per-signal contexts through a recursive walk. That has to fold at
    # inference time like the rest of the hot path, and must not put the
    # combining machinery on `track!`'s allocation budget.
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 25_000)
    measurements = (L1 = Tracking.BandMeasurement(buf, FS, 0.0Hz),)

    for nsig in (1, 2, 3)
        allocations = map((false, true)) do combining
            ts = identical_l1ca_track_state(
                nsig,
                prn,
                carrier_doppler;
                discriminator_combining = combining,
            )
            dc = Tracking.CPUDownconvertAndCorrelator()
            Tracking.track!(measurements, ts; downconvert_and_correlator = dc)  # warmup
            Tracking.track!(measurements, ts; downconvert_and_correlator = dc)
            allocated = @allocated Tracking.track!(
                measurements,
                ts;
                downconvert_and_correlator = dc,
            )
            # Inference of the estimator step itself, on a state holding real
            # unconsumed correlator outputs.
            Tracking.downconvert_and_correlate!(
                dc,
                measurements,
                ts;
                chunk_index = 0,
                chunk_duration = 1e-3s,
                stop_before_partial = true,
            )
            @inferred Tracking.estimate_dopplers_and_filter_prompt!(ts, measurements)
            allocated
        end
        # Combining is free: identical allocation to the driver-only path.
        @test allocations[1] == allocations[2]
    end

    # The same pin on a **heterogeneous** group, which is the realistic one and
    # the only one that can catch a regression in the tuple walk: three distinct
    # signal types, two distinct correlator types, and per-passenger contexts
    # whose group-delay states differ.
    het_fs = 15e6Hz
    het_buf = l1ca_signal(prn, carrier_doppler, 0.0, 45_000, het_fs)
    het_measurements = (L1 = Tracking.BandMeasurement(het_buf, het_fs, 0.0Hz),)
    het_allocations = map((false, true)) do combining
        estimator = ConventionalAssistedPLLAndDLL()
        # The allocation contract is about the fold's types, not its numerics:
        # this ordering is what puts three different signal types in one group at
        # a sampling rate all three can be generated at.
        sigs = (
            TrackedSignal(
                GPSL1CA();
                correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                group_delay = 0.0s,
            ),
            TrackedSignal(
                GPSL1C_D();
                correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                group_delay = 1.0e-9s,
            ),
            TrackedSignal(
                GPSL1C_P();
                correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                # Left unknown on purpose: the carrier-only passenger.
                group_delay = nothing,
            ),
        )
        sat = TrackedSat(sigs, prn, 0.0, carrier_doppler; doppler_estimator = estimator)
        ts = TrackState(
            GPSL1CA(),
            sat;
            doppler_estimator = estimator,
            discriminator_combining = combining,
        )
        dc = Tracking.CPUDownconvertAndCorrelator()
        Tracking.track!(het_measurements, ts; downconvert_and_correlator = dc)  # warmup
        Tracking.track!(het_measurements, ts; downconvert_and_correlator = dc)
        allocated = @allocated Tracking.track!(
            het_measurements,
            ts;
            downconvert_and_correlator = dc,
        )
        Tracking.downconvert_and_correlate!(
            dc,
            het_measurements,
            ts;
            chunk_index = 0,
            chunk_duration = 1e-3s,
            stop_before_partial = true,
        )
        @inferred Tracking.estimate_dopplers_and_filter_prompt!(ts, het_measurements)
        allocated
    end
    @test het_allocations[1] == het_allocations[2]
end

# A CN0 estimator that keeps the noise density it was last handed, so which floor
# reached which signal is observable from the outside.
struct DensityRecorder{D} <: Tracking.AbstractCN0Estimator
    density::D
end
Tracking.update(::DensityRecorder, prompt, ctx::Tracking.CN0UpdateContext) =
    DensityRecorder(ctx.noise_density)

@testset "every signal folds against its own noise floor" begin
    # The fold threads one `(density, ready)` pair per signal because the C/N₀
    # floor is per signal — handing one passenger another's would misreport every
    # record it folds.
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    signals = ntuple(
        _ -> TrackedSignal(
            TrackedSignal(gpsl1; correlator, cn0_estimator = DensityRecorder(nothing));
            correlator_outputs = [CorrelatorOutput(correlator, n, n)],
        ),
        2,
    )
    sat = TrackedSat(signals, 1, 0.0, 1000.0Hz)
    estimator = ConventionalAssistedPLLAndDLL()
    state = Tracking.init_estimator_state(estimator, sat)

    # Two floors five times apart, so a mix-up cannot hide in rounding.
    driver_density, passenger_density = 1.0e-7 / 1.0Hz, 5.0e-7 / 1.0Hz
    noise = ((driver_density, true), (passenger_density, true))
    for combine in (Val(true), Val(false))
        # Fresh buffers per run: the fold empties them.
        for s in signals
            empty!(s.correlator_outputs)
            push!(s.correlator_outputs, CorrelatorOutput(correlator, n, n))
        end
        new_driver, new_passengers, _, _, _ = Tracking._fold_satellite_chunk(
            sat.signals[1],
            Base.tail(sat.signals),
            sat,
            state,
            FS,
            noise,
            get_carrier_phase_offset(gpsl1),
            combine,
            nothing,
        )
        @test Tracking.get_cn0_estimator(new_driver).density ≈ driver_density
        @test Tracking.get_cn0_estimator(new_passengers[1]).density ≈ passenger_density
    end
end

# --------------------------------------------------------------------------
# The combined value actually reaches the loop filters
# --------------------------------------------------------------------------
#
# The testsets above pin the arithmetic of one record's contribution, or compare
# a combined run against a driver-only one and expect them to *agree* — which two
# copies of one signal make true whether the combiner runs or not. These feed
# the two signals records that disagree and exploit the fact that the loop
# filters and `aid_dopplers` are affine in the discriminator on a fresh filter
# state: if `d` is what the driver's record alone produces and `p` is what the
# passenger's record alone produces (installed on the driver, so the same code
# path evaluates it), then a weight-normalized mean must land exactly at the
# weighted point between them. No loop internals are reproduced here, which is
# what keeps the assertion honest.
#
# Records are appended through `append_correlator_output!` rather than correlated
# from samples, so each run is one record per signal at a known sample index and
# the arithmetic is exact rather than convergence-limited.

# Two records of one satellite that disagree about both the code error (E/L
# asymmetry, opposite ways round) and the carrier phase error (prompt argument,
# opposite signs). Same correlator layout and same sample count, so on two
# signals of equal power share they carry equal weight in all three loops and
# the mean is the midpoint.
const DISAGREEING_DRIVER_RECORD = EarlyPromptLateCorrelator(
    [complex(0.90, 0.00), complex(1.00, 0.10), complex(0.60, 0.00)],
    0.5,
)
const DISAGREEING_PASSENGER_RECORD = EarlyPromptLateCorrelator(
    [complex(0.55, 0.00), complex(1.00, -0.25), complex(0.95, 0.00)],
    0.5,
)
const RECORD_SAMPLES = 5000

# A two-signal satellite whose records arrive from an external producer.
# `NoCN0Estimator` keeps the noise-density plumbing out of it — this is a test
# about the loop update, and the C/N₀ path has its own.
function externally_fed_state(;
    combining = true,
    passenger_delay = 0.0s,
    signals = (GPSL1CA(), GPSL1CA()),
    cn0_estimator = NoCN0Estimator,
)
    estimator = ConventionalAssistedPLLAndDLL()
    tracked_signals = map(
        (signal, delay) -> TrackedSignal(
            signal;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            cn0_estimator = cn0_estimator(),
            group_delay = delay,
        ),
        signals,
        (0.0s, passenger_delay),
    )
    sat = TrackedSat(tracked_signals, 1, 0.0, 1000.0Hz; doppler_estimator = estimator)
    TrackState(
        first(signals),
        sat;
        doppler_estimator = estimator,
        discriminator_combining = combining,
    )
end

# One estimate call over one driver record and, optionally, one passenger record
# ending at `passenger_sample_index` (by default the driver's). Returns the
# satellite's Dopplers after it.
function one_update(
    driver_record;
    passenger_record = nothing,
    passenger_sample_index = RECORD_SAMPLES,
    kwargs...,
)
    ts = externally_fed_state(; kwargs...)
    n = RECORD_SAMPLES
    isnothing(passenger_record) || append_correlator_output!(
        ts,
        CorrelatorOutput(passenger_record, n, passenger_sample_index),
        1,
        1,
        2,
    )
    append_correlator_output!(ts, CorrelatorOutput(driver_record, n, n), 1, 1, 1)
    estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
    (get_carrier_doppler(ts, 1), get_code_doppler(ts, 1))
end

# The code loop is carrier-aided (`aid_dopplers`), so the raw code Doppler
# carries the PLL's update too. Subtracting that term leaves the quantity that is
# affine in the DLL discriminator *alone*, which is what the group-delay
# assertions below need.
_code_only(dopplers) =
    dopplers[2] - (dopplers[1] - 1000.0Hz) * get_code_center_frequency_ratio(GPSL1CA())

@testset "the combined discriminator is what the loop filters see" begin
    driver_alone = one_update(DISAGREEING_DRIVER_RECORD)
    # The passenger's record run through the *driver's* slot: the same loop code
    # evaluated at the passenger's measurement, which is the other end of the
    # interval the mean has to land inside.
    passenger_alone = one_update(DISAGREEING_PASSENGER_RECORD)
    combined = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )

    # The two records really do disagree, so the midpoint below is nowhere near
    # either end and no tolerance can hide a combiner that dropped one of them.
    @test abs(driver_alone[1] - passenger_alone[1]) > 1.0Hz
    @test abs(driver_alone[2] - passenger_alone[2]) > 0.1Hz

    # Equal weights, so the combined update is exactly the midpoint.
    @test combined[1] ≈ (driver_alone[1] + passenger_alone[1]) / 2 rtol = 1e-12
    @test combined[2] ≈ (driver_alone[2] + passenger_alone[2]) / 2 rtol = 1e-12

    # …and with combining off the passenger's record must not reach the loop at
    # all, bit-identically to the run where it was never appended.
    driver_only_run = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        combining = false,
    )
    @test driver_only_run === driver_alone
end

@testset "coincident records combine as a weighted mean, by power share" begin
    # GPS L1C splits its power 75/25 between pilot and data, and both components
    # share a 10 ms code period, so their records end on the same sample. With the
    # same correlator layout on both, every loop weights the pilot driver 3:1
    # against the data passenger, and the update lands a quarter of the way from
    # the driver's own towards the passenger's — not at the midpoint.
    l1c = (GPSL1C_P(), GPSL1C_D())
    driver_alone = one_update(DISAGREEING_DRIVER_RECORD; signals = l1c)
    passenger_alone = one_update(DISAGREEING_PASSENGER_RECORD; signals = l1c)
    combined = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        signals = l1c,
    )
    @test abs(driver_alone[1] - passenger_alone[1]) > 1.0Hz
    for loop = 1:2
        @test combined[loop] ≈
              driver_alone[loop] + (passenger_alone[loop] - driver_alone[loop]) / 4 rtol =
            1e-12
    end
end

# A CN0 estimator that counts the records it was handed.
struct CountingCN0 <: Tracking.AbstractCN0Estimator
    count::Int
end
CountingCN0() = CountingCN0(0)
Tracking.update(e::CountingCN0, prompt, ::Tracking.CN0UpdateContext) =
    CountingCN0(e.count + 1)

@testset "a record that coincides with no driver record is applied, never combined" begin
    # The passenger's record ends one sample before the driver's. It is applied
    # to its own signal in full — prompt, post-correlation filter, C/N₀
    # estimator, bit buffer, `last_fully_integrated_*` — but its discriminators
    # are dropped, so the loop update is bit-identical to the driver alone.
    driver_alone = one_update(DISAGREEING_DRIVER_RECORD)
    shifted = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_sample_index = RECORD_SAMPLES - 1,
    )
    @test shifted === driver_alone
    # One sample *late* is not combined either: it is drained after the last
    # driver record, and nothing carries it to the next call.
    late = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_sample_index = RECORD_SAMPLES + 1,
    )
    @test late === driver_alone

    # The passenger's own state advances exactly as it does when the same record
    # *is* combined: whether a record combines changes the loop update and
    # nothing else.
    n = RECORD_SAMPLES
    passenger_after(sample_index) = begin
        ts = externally_fed_state(; cn0_estimator = CountingCN0)
        append_correlator_output!(
            ts,
            CorrelatorOutput(DISAGREEING_PASSENGER_RECORD, n, sample_index),
            1,
            1,
            2,
        )
        append_correlator_output!(
            ts,
            CorrelatorOutput(DISAGREEING_DRIVER_RECORD, n, n),
            1,
            1,
            1,
        )
        before = get_sat_state(ts, 1).signals[2]
        bit_buffer_before = Tracking.get_bit_buffer(before)
        estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
        (bit_buffer_before, get_sat_state(ts, 1).signals[2])
    end
    bit_buffer_before, not_coincident = passenger_after(n - 1)
    _, coincident = passenger_after(n)
    @test isempty(Tracking.get_correlator_outputs(not_coincident))
    @test Tracking.get_filtered_prompts(not_coincident) ==
          Tracking.get_filtered_prompts(coincident)
    @test length(Tracking.get_filtered_prompts(not_coincident)) == 1
    @test Tracking.get_last_fully_integrated_filtered_prompt(not_coincident) ==
          Tracking.get_last_fully_integrated_filtered_prompt(coincident) !=
          complex(0.0, 0.0)
    @test Tracking.get_last_fully_integrated_correlator(not_coincident) ==
          DISAGREEING_PASSENGER_RECORD
    @test Tracking.get_cn0_estimator(not_coincident).count == 1
    @test Tracking.get_cn0_estimator(coincident).count == 1
    # The bit buffer's sync search saw the record: one more code block buffered,
    # the same one either way.
    bits(signal) = (
        Tracking.get_bit_buffer(signal).code_block_buffer,
        Tracking.get_bit_buffer(signal).code_block_buffer_length,
    )
    @test bits(not_coincident) == bits(coincident)
    @test bits(not_coincident)[2] == bit_buffer_before.code_block_buffer_length + 1
end

@testset "a passenger integrating longer than the driver aids only where records meet" begin
    # The equal-length assumption is documented, not enforced. A passenger that
    # integrates twice as long as the driver completes one record per two driver
    # records, and only the driver record ending on the same sample is combined.
    n = RECORD_SAMPLES
    run(passenger_record; combining = true) = begin
        ts = externally_fed_state(; combining)
        isnothing(passenger_record) || append_correlator_output!(
            ts,
            CorrelatorOutput(passenger_record, 2n, 2n),
            1,
            1,
            2,
        )
        append_correlator_output!(
            ts,
            CorrelatorOutput(DISAGREEING_DRIVER_RECORD, n, n),
            1,
            1,
            1,
        )
        append_correlator_output!(
            ts,
            CorrelatorOutput(DISAGREEING_DRIVER_RECORD, n, 2n),
            1,
            1,
            1,
        )
        estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
        (get_carrier_doppler(ts, 1), get_code_doppler(ts, 1))
    end
    driver_alone = run(nothing)
    @test run(DISAGREEING_PASSENGER_RECORD; combining = false) === driver_alone
    # The second driver record met the passenger's, so the chunk's final update
    # moved.
    @test run(DISAGREEING_PASSENGER_RECORD) !== driver_alone

    # Nothing refuses the configuration, either at construction or through the
    # setter — the setter is about one signal's integration length and says
    # nothing about what combining makes of it.
    prn, carrier_doppler = 1, 1000.0Hz
    ts = identical_l1ca_track_state(2, prn, carrier_doppler)
    set_preferred_num_code_blocks_to_integrate!(ts, 1, prn, 2, 2)
    @test get_preferred_num_code_blocks_to_integrate(get_sat_state(ts, prn), 2) == 2
    set_preferred_num_code_blocks_to_integrate!(ts, 1, prn, 1, 2)
    @test get_preferred_num_code_blocks_to_integrate(get_sat_state(ts, prn), 1) == 2
end

@testset "a supplied group delay moves the tracked code phase, by the right amount" begin
    gpsl1 = GPSL1CA()
    driver_alone = one_update(DISAGREEING_DRIVER_RECORD)
    passenger_alone = one_update(DISAGREEING_PASSENGER_RECORD)
    combined = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )

    # The code Doppler every record in these runs was evaluated at (the
    # satellite's, carrier-derived — `dll_disc` is fed the chunk-fixed value).
    code_doppler = 1000.0Hz * get_code_center_frequency_ratio(gpsl1)
    disc_driver = Tracking.dll_disc(gpsl1, DISAGREEING_DRIVER_RECORD, code_doppler, FS)
    disc_passenger =
        Tracking.dll_disc(gpsl1, DISAGREEING_PASSENGER_RECORD, code_doppler, FS)
    # Volts per chip of the code loop, measured rather than derived: the two
    # single-record runs above give two points on an affine map.
    per_chip =
        (_code_only(passenger_alone) - _code_only(driver_alone)) /
        (disc_passenger - disc_driver)

    δ = 4.0e-9s
    δchips = uconvert(NoUnits, δ * (get_code_frequency(gpsl1) + code_doppler))
    delayed = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_delay = δ,
    )

    # A larger group delay than the driver's says the passenger arrives later and
    # so sits at the *smaller* code phase, so `δchips` is added back to *its*
    # discriminator only — and the mean of two equally weighted records therefore
    # moves by half of it.
    expected =
        _code_only(driver_alone) +
        per_chip * ((disc_driver + disc_passenger + δchips) / 2 - disc_driver)
    @test _code_only(delayed) ≈ expected rtol = 1e-12
    @test _code_only(delayed) > _code_only(combined)
    # The carrier loops never see the delay.
    @test delayed[1] === combined[1]

    # Linear in the delay.
    twice = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_delay = 2δ,
    )
    @test _code_only(twice) - _code_only(combined) ≈
          2 * (_code_only(delayed) - _code_only(combined)) rtol = 1e-12

    # An unknown delay holds the passenger out of the code loop: the carrier
    # update is the combined one, the code-only part is the driver's.
    unknown = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_delay = nothing,
    )
    @test unknown[1] === combined[1]
    @test _code_only(unknown) ≈ _code_only(driver_alone) rtol = 1e-12
end

@testset "combining is a group property, and a mis-ordered group is accepted" begin
    # Combining is declared on the `SignalGroup`, next to its band and antenna
    # count, because its assumption — `signals[1]` integrates at least as long as
    # every passenger — is a property of the group's signal tuple, which every
    # satellite of the group shares. Nothing per-satellite carries it.
    combining_state = externally_fed_state(; combining = true)
    driver_only_state = externally_fed_state(; combining = false)
    @test only(Tuple(combining_state.groups)).discriminator_combining
    @test !only(Tuple(driver_only_state.groups)).discriminator_combining
    @test !hasfield(Tracking.SatConventionalPLLAndDLL, :discriminator_combining)
    @test !hasfield(Tracking.ConventionalPLLAndDLL, :discriminator_combining)

    # The default is off, everywhere.
    @test !only(Tuple(TrackState(; signal = GPSL1CA()).groups)).discriminator_combining
    @test !Tracking.SignalGroup((GPSL1C_P(), GPSL1C_D())).discriminator_combining
    # A pre-built group keeps its own flag; the `TrackState` kwarg is only the
    # default for bare signal tuples.
    prebuilt = TrackState(;
        signals = (
            a = Tracking.SignalGroup(
                (GPSL1C_P(), GPSL1C_D());
                discriminator_combining = true,
            ),
            b = (GPSL1CA(),),
        ),
    )
    @test prebuilt.groups.a.discriminator_combining
    @test !prebuilt.groups.b.discriminator_combining
    # The kwarg-update constructor carries it too.
    @test !Tracking.SignalGroup(prebuilt.groups.a; discriminator_combining = false).discriminator_combining
    @test Tracking.SignalGroup(prebuilt.groups.a).discriminator_combining

    # The ordering is a documented assumption and not an enforced one: a 1 ms
    # C/A driver with a 10 ms L1C-P passenger costs gain, never correctness, so
    # it is built rather than refused.
    @test TrackState(;
        signals = (gps_l1 = (GPSL1CA(), GPSL1C_P()),),
        discriminator_combining = true,
    ) isa TrackState
    @test Tracking.SignalGroup((GPSL1CA(), GPSL1C_P()); discriminator_combining = true) isa
          Tracking.SignalGroup
end

# --------------------------------------------------------------------------
# Vector tracking: the combination reaches exactly the loops still closed here
# --------------------------------------------------------------------------
#
# `VectorPLLAndDLL` combines every loop it closes itself and no others, so the
# reach follows `vt_on`:
#
#   * `vt_on = false` — scalar fallback, all three discriminators combine, and
#     the satellite must behave like the conventional estimator does.
#   * `vt_on = true` — the carrier phase loop alone. The code and carrier
#     frequency loops are the navigation filter's, and every signal hands it its
#     own raw measurement.

# The vector-tracking twin of `externally_fed_state`: the same two-signal L1 C/A
# satellite fed by an external producer, under `VectorPLLAndDLL`.
function vt_externally_fed_state(; combining = true, vt_on = false, passenger_delay = 0.0s)
    estimator = VectorPLLAndDLL()
    gpsl1 = GPSL1CA()
    signals = (
        TrackedSignal(
            gpsl1;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            cn0_estimator = NoCN0Estimator(),
            group_delay = 0.0s,
        ),
        TrackedSignal(
            gpsl1;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            cn0_estimator = NoCN0Estimator(),
            group_delay = passenger_delay,
        ),
    )
    sat = TrackedSat(signals, 1, 0.0, 1000.0Hz; doppler_estimator = estimator)
    ts = TrackState(
        gpsl1,
        sat;
        doppler_estimator = estimator,
        discriminator_combining = combining,
    )
    vt_on && enable_vt!(ts, [1])
    ts
end

get_vt_state(ts) = Tracking.get_doppler_estimator_state(get_sat_state(ts, 1))

# `one_update`'s twin: one estimate call over one driver record and, optionally,
# one passenger record ending at `passenger_sample_index`. Returns the Dopplers
# and the state the navigation filter would read.
function vt_one_update(
    driver_record;
    passenger_record = nothing,
    passenger_sample_index = RECORD_SAMPLES,
    kwargs...,
)
    ts = vt_externally_fed_state(; kwargs...)
    n = RECORD_SAMPLES
    isnothing(passenger_record) || append_correlator_output!(
        ts,
        CorrelatorOutput(passenger_record, n, passenger_sample_index),
        1,
        1,
        2,
    )
    append_correlator_output!(ts, CorrelatorOutput(driver_record, n, n), 1, 1, 1)
    estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
    (
        carrier_doppler = get_carrier_doppler(ts, 1),
        code_doppler = get_code_doppler(ts, 1),
        state = get_vt_state(ts),
    )
end

# Feed one record to each signal of the vector fixture and run one estimate
# call — `vt_one_update`'s repeatable form, for the properties that only show up
# across calls (an FLL measurement needs a previous prompt, so it needs two).
function vt_feed!(ts, driver_record, passenger_record)
    n = RECORD_SAMPLES
    isnothing(passenger_record) ||
        append_correlator_output!(ts, CorrelatorOutput(passenger_record, n, n), 1, 1, 2)
    isnothing(driver_record) ||
        append_correlator_output!(ts, CorrelatorOutput(driver_record, n, n), 1, 1, 1)
    estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
    get_vt_state(ts)
end

@testset "the scalar fallback combines every loop, as the conventional estimator does" begin
    driver_alone = vt_one_update(DISAGREEING_DRIVER_RECORD)
    passenger_alone = vt_one_update(DISAGREEING_PASSENGER_RECORD)
    combined = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )

    @test abs(driver_alone.carrier_doppler - passenger_alone.carrier_doppler) > 1.0Hz
    @test abs(driver_alone.code_doppler - passenger_alone.code_doppler) > 0.1Hz
    @test combined.carrier_doppler ≈
          (driver_alone.carrier_doppler + passenger_alone.carrier_doppler) / 2 rtol = 1e-12
    @test combined.code_doppler ≈
          (driver_alone.code_doppler + passenger_alone.code_doppler) / 2 rtol = 1e-12

    # The strongest statement available: the same satellite under
    # `ConventionalAssistedPLLAndDLL`, fed the same two records, lands on the same
    # Dopplers. The fallback is not "like" scalar tracking, it is scalar tracking
    # — same fold, same weighting, same filters.
    conventional = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )
    @test combined.carrier_doppler === conventional[1]
    @test combined.code_doppler === conventional[2]

    # With combining off, or with a record that does not coincide, the passenger
    # does not reach the loop at all.
    for kwargs in ((; combining = false), (; passenger_sample_index = RECORD_SAMPLES - 1))
        driver_only_run = vt_one_update(
            DISAGREEING_DRIVER_RECORD;
            passenger_record = DISAGREEING_PASSENGER_RECORD,
            kwargs...,
        )
        @test driver_only_run.carrier_doppler === driver_alone.carrier_doppler
        @test driver_only_run.code_doppler === driver_alone.code_doppler
    end
end

@testset "under vector closure only the carrier phase loop combines" begin
    driver_alone = vt_one_update(DISAGREEING_DRIVER_RECORD; vt_on = true)
    combined = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
    )

    @test combined.state.vt_on
    # Slot 1 is the driver's, and it reads exactly what it reads when it is the
    # only signal that reported.
    @test mean_code_discr(combined.state, 1) === mean_code_discr(driver_alone.state, 1)
    @test mean_carrier_discr(combined.state, 1) ===
          mean_carrier_discr(driver_alone.state, 1)
    @test combined.state.code_discr_acc[1] == (1, mean_code_discr(driver_alone.state, 1))
    # The code NCO follows the navigation filter's correction, which is zero here,
    # whatever the passenger measured: only the carrier aiding moves it.
    @test combined.code_doppler -
          combined.carrier_doppler * get_code_center_frequency_ratio(GPSL1CA()) ===
          driver_alone.code_doppler -
          driver_alone.carrier_doppler * get_code_center_frequency_ratio(GPSL1CA())

    # …while the carrier phase loop, which is still the satellite's own under
    # vector closure, did see the passenger.
    passenger_alone = vt_one_update(DISAGREEING_PASSENGER_RECORD; vt_on = true)
    @test combined.carrier_doppler ≈
          (driver_alone.carrier_doppler + passenger_alone.carrier_doppler) / 2 rtol = 1e-12
    # And not with combining off.
    separate = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
        combining = false,
    )
    @test separate.carrier_doppler === driver_alone.carrier_doppler
end

@testset "every signal dumps its own discriminators for the navigation filter" begin
    driver_alone = vt_one_update(DISAGREEING_DRIVER_RECORD; vt_on = true)
    passenger_alone = vt_one_update(DISAGREEING_PASSENGER_RECORD; vt_on = true)
    combined = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
    )

    # `passenger_alone` ran the passenger's record through the *driver's* slot,
    # so its slot-1 dump is what that record measures — the value the passenger
    # must report from slot 2.
    @test mean_code_discr(combined.state, 1) === mean_code_discr(driver_alone.state, 1)
    @test mean_code_discr(combined.state, 2) === mean_code_discr(passenger_alone.state, 1)
    @test mean_carrier_discr(combined.state, 2) ===
          mean_carrier_discr(passenger_alone.state, 1)
    # A passenger record that coincides with no driver record is still measured:
    # collection does not depend on combining.
    shifted = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_sample_index = RECORD_SAMPLES - 1,
        vt_on = true,
    )
    @test mean_code_discr(shifted.state, 2) === mean_code_discr(combined.state, 2)

    # …and the same values through the addressing ladder every other per-signal
    # accessor uses, including from the `TrackState`.
    ts = vt_externally_fed_state(; vt_on = true)
    vt_feed!(ts, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)
    @test mean_code_discr(ts, 1, 1, 2) === mean_code_discr(combined.state, 2)
    @test mean_code_discr(get_sat_state(ts, 1), 2) === mean_code_discr(combined.state, 2)
    @test mean_carrier_discr(ts, 1, 1, 2) === mean_carrier_discr(combined.state, 2)
    # Both signals are GPS L1 C/A here, so a type selector is ambiguous and must
    # say so rather than pick one — as it does for every other accessor.
    @test_throws ArgumentError mean_code_discr(get_sat_state(ts, 1), GPSL1CA)
    # An unqualified read of a multi-signal satellite is refused for the same
    # reason: answering with the driver's is the mistake per-signal accumulation
    # exists to remove.
    @test_throws ArgumentError mean_code_discr(get_vt_state(ts))
    @test_throws ArgumentError mean_carrier_discr(ts, 1)
    @test mean_code_discr(combined.state, 1) != mean_code_discr(combined.state, 2)

    # One record each, counted separately: the filter divides by these.
    @test first.(combined.state.code_discr_acc) == (1, 1)
    @test last.(combined.state.code_discr_acc) ==
          (mean_code_discr(combined.state, 1), mean_code_discr(combined.state, 2))
    @test first.(combined.state.carrier_discr_acc) == (1, 1)

    # A satellite still in the scalar fallback dumps nothing at all, on any slot.
    fallback = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )
    @test isnothing(mean_code_discr(fallback.state, 1))
    @test isnothing(mean_code_discr(fallback.state, 2))

    # And every slot resets together.
    @test !isnothing(mean_code_discr(get_vt_state(ts), 2))
    reset_code_discr_acc!(ts)
    reset_carrier_discr_acc!(ts)
    @test all(
        isnothing,
        (
            mean_code_discr(get_vt_state(ts), 1),
            mean_code_discr(get_vt_state(ts), 2),
            mean_carrier_discr(get_vt_state(ts), 1),
            mean_carrier_discr(get_vt_state(ts), 2),
        ),
    )
end

@testset "a dumped code measurement is raw, whatever group delay is supplied" begin
    # The group delay belongs to the *combined* code loop, and under `vt_on` that
    # loop is the navigation filter's — so a dump is what the signal measured,
    # and setting or clearing a delay must not move it.
    δt = 4.0e-9s
    with_delay = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
        passenger_delay = δt,
    )
    aligned = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
    )
    @test mean_code_discr(with_delay.state, 2) === mean_code_discr(aligned.state, 2)
    @test mean_code_discr(with_delay.state, 1) === mean_code_discr(aligned.state, 1)
    @test mean_carrier_discr(with_delay.state, 2) === mean_carrier_discr(aligned.state, 2)

    # And the value is the passenger's own discriminator, not the driver's.
    @test mean_code_discr(aligned.state, 2) ≈ Tracking.dll_disc(
        GPSL1CA(),
        DISAGREEING_PASSENGER_RECORD,
        get_code_doppler(get_sat_state(vt_externally_fed_state(), 1)),
        FS,
    )

    # A passenger with no delay supplied is folded like any other — nothing is
    # withheld, so `count` is the record count.
    unknown = vt_one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        vt_on = true,
        passenger_delay = nothing,
    )
    @test mean_code_discr(unknown.state, 2) === mean_code_discr(aligned.state, 2)
    @test unknown.state.code_discr_acc[2][1] == 1
    @test mean_carrier_discr(unknown.state, 2) === mean_carrier_discr(aligned.state, 2)
end

@testset "every record reaches both accumulators under vector tracking" begin
    # `fll_disc` answers 0 Hz when there is no previous prompt to difference
    # against, and that placeholder is counted like any other measurement — what
    # keeps it out of a navigation filter's mean is that `vt_on` is set long
    # after a signal's first record. Here the flag is set by hand on the first
    # integration, which is why the placeholder shows.
    ts = vt_externally_fed_state(; vt_on = true)
    after_first = vt_feed!(ts, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)
    @test first.(after_first.code_discr_acc) == (1, 1)
    @test first.(after_first.carrier_discr_acc) == (1, 1)
    @test mean_carrier_discr(after_first, 1) == 0.0Hz
    @test mean_carrier_discr(after_first, 2) == 0.0Hz

    # From the second record on there is a prompt pair on both signals, and the
    # two loops stay in step: one record in, one record counted, on every signal
    # and both loops. The records are swapped between the slots so each signal's
    # prompt really rotates.
    after_second = vt_feed!(ts, DISAGREEING_PASSENGER_RECORD, DISAGREEING_DRIVER_RECORD)
    @test first.(after_second.code_discr_acc) == (2, 2)
    @test first.(after_second.carrier_discr_acc) == (2, 2)
    @test all(!iszero, last.(after_second.carrier_discr_acc))
    for signal_index = 1:2
        count, discr_sum = after_second.carrier_discr_acc[signal_index]
        @test mean_carrier_discr(after_second, signal_index) == discr_sum / count
    end
end

@testset "measurements are collected under vt_on whether or not signals combine" begin
    separate = vt_externally_fed_state(; combining = false, vt_on = true)
    vt_feed!(separate, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)
    state = vt_feed!(separate, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)

    # Every slot filled, each with its own signal's measurement.
    @test first.(state.code_discr_acc) == (2, 2)
    @test first.(state.carrier_discr_acc) == (2, 2)
    @test mean_code_discr(state, 1) != mean_code_discr(state, 2)

    # The passenger's slot holds the passenger's own value — the same value that
    # record produces when it is the only signal reporting. Approximately: in the
    # reference run that record is the driver, so its own loops moved the code
    # Doppler every `dll_disc` is evaluated at by a few parts in 1e9.
    passenger_alone = vt_externally_fed_state(; combining = false, vt_on = true)
    vt_feed!(passenger_alone, DISAGREEING_PASSENGER_RECORD, nothing)
    passenger_state = vt_feed!(passenger_alone, DISAGREEING_PASSENGER_RECORD, nothing)
    @test mean_code_discr(state, 2) ≈ mean_code_discr(passenger_state, 1) rtol = 1e-6
    @test mean_carrier_discr(state, 2) ≈ mean_carrier_discr(passenger_state, 1) rtol = 1e-6
    @test !isapprox(mean_code_discr(state, 1), mean_code_discr(state, 2), rtol = 1e-6)

    # …and nothing the passengers measured reached the loops, which are still the
    # driver's alone.
    driver_alone = vt_externally_fed_state(; combining = false, vt_on = true)
    vt_feed!(driver_alone, DISAGREEING_DRIVER_RECORD, nothing)
    vt_feed!(driver_alone, DISAGREEING_DRIVER_RECORD, nothing)
    @test get_carrier_doppler(separate, 1) === get_carrier_doppler(driver_alone, 1)
    @test get_code_doppler(separate, 1) === get_code_doppler(driver_alone, 1)

    # A satellite outside the vector loop collects nothing, combining or not.
    fallback = vt_externally_fed_state(; combining = false, vt_on = false)
    vt_feed!(fallback, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)
    fallback_state =
        vt_feed!(fallback, DISAGREEING_DRIVER_RECORD, DISAGREEING_PASSENGER_RECORD)
    @test first.(fallback_state.code_discr_acc) == (0, 0)
    @test first.(fallback_state.carrier_discr_acc) == (0, 0)
end

@testset "a pre-sync-invalidated passenger record skips both slots and the combination" begin
    # A passenger record correlated with a pre-sync replica on a secondary-coded
    # signal carries no usable phase or delay, whether it is about to be averaged
    # or shipped: it is excluded from the combination and from both of its
    # accumulator slots together, so the two counts never diverge. A signal
    # without a secondary code is unaffected.
    n = RECORD_SAMPLES
    synced(b::Tracking.BitBuffer{B}) where {B} = Tracking.BitBuffer{B}(
        zero(B),
        0,
        true,
        0,
        Int8(+1),
        complex(0.0, 0.0),
        0,
        Float32[],
        Tracking.PhaseAccumulators(),
    )
    zero_measurements =
        Tracking.PerSignalMeasurements(((0, 0.0), (0, 0.0)), ((0, 0.0Hz), (0, 0.0Hz)))
    for (signal, excluded) in ((GalileoE1C(), true), (GPSL1CA(), false))
        base = TrackedSignal(
            signal;
            correlator = DISAGREEING_PASSENGER_RECORD,
            cn0_estimator = NoCN0Estimator(),
            group_delay = 0.0s,
        )
        # Synced *now*, not at the chunk boundary: the record follows a sync
        # detected earlier in this same fold.
        ts = TrackedSignal(
            base;
            bit_buffer = synced(Tracking.get_bit_buffer(base)),
            correlator_outputs = [CorrelatorOutput(DISAGREEING_PASSENGER_RECORD, n, n)],
        )
        ctx = RecordFoldContext(1, 2, FS, 0.0Hz, one(ComplexF64), 0.0, false, nothing, true)
        _, _, _, was_excluded, measurements = Tracking._apply_passenger_record(
            ts,
            ts.correlator_outputs[1],
            ctx,
            zero_measurements,
        )
        @test was_excluded == excluded
        @test first.(measurements.code_discr_acc) == (0, excluded ? 0 : 1)
        @test first.(measurements.carrier_discr_acc) == (0, excluded ? 0 : 1)
        _, _, contribution, _ = Tracking._advance_one_passenger_to(ts, 1, ctx, n, nothing)
        @test (contribution == DiscriminatorAccumulator()) == excluded
    end
end

@testset "vector combining reads the group's flag, like the conventional estimator" begin
    combining_state = vt_externally_fed_state(; combining = true, vt_on = true)
    driver_only_state = vt_externally_fed_state(; combining = false, vt_on = true)
    @test only(Tuple(combining_state.groups)).discriminator_combining
    @test !only(Tuple(driver_only_state.groups)).discriminator_combining
    @test !hasfield(Tracking.SatVectorPLLAndDLL, :discriminator_combining)
    @test !hasfield(Tracking.VectorPLLAndDLL, :discriminator_combining)
    # Nothing is parked between calls under vector tracking either.
    @test !hasfield(Tracking.SatVectorPLLAndDLL, :pending_combining_sums)
end

@testset "a single-signal vector-tracked sat is bit-identical with combining on" begin
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 20_000)
    for vt_on in (false, true)
        dopplers = map((false, true)) do combining
            ts = TrackState(;
                signal = GPSL1CA(),
                doppler_estimator = VectorPLLAndDLL(),
                discriminator_combining = combining,
            )
            ts = add_satellite!(ts; prn, code_phase = 0.0, carrier_doppler)
            vt_on && enable_vt!(ts, [prn])
            sat = get_sat_state(track(copy(buf), ts, FS), prn)
            (
                get_carrier_doppler(sat),
                get_code_doppler(sat),
                Tracking.get_code_phase(sat),
                Tracking.get_doppler_estimator_state(sat).code_discr_acc,
            )
        end
        @test dopplers[1] === dopplers[2]
    end
end

@testset "the vector combining fold is inferable and adds no allocations" begin
    prn, carrier_doppler = 1, 1000.0Hz
    het_fs = 15e6Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 45_000, het_fs)
    measurements = (L1 = Tracking.BandMeasurement(buf, het_fs, 0.0Hz),)

    for vt_on in (false, true)
        allocations = map((false, true)) do combining
            estimator = VectorPLLAndDLL()
            sigs = (
                TrackedSignal(
                    GPSL1CA();
                    correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                    group_delay = 0.0s,
                ),
                TrackedSignal(
                    GPSL1C_D();
                    correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                    group_delay = 1.0e-9s,
                ),
                TrackedSignal(
                    GPSL1C_P();
                    correlator = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
                ),
            )
            sat = TrackedSat(sigs, prn, 0.0, carrier_doppler; doppler_estimator = estimator)
            ts = TrackState(
                GPSL1CA(),
                sat;
                doppler_estimator = estimator,
                discriminator_combining = combining,
            )
            vt_on && enable_vt!(ts, [prn])
            dc = Tracking.CPUDownconvertAndCorrelator()
            Tracking.track!(measurements, ts; downconvert_and_correlator = dc)  # warmup
            Tracking.track!(measurements, ts; downconvert_and_correlator = dc)
            allocated = @allocated Tracking.track!(
                measurements,
                ts;
                downconvert_and_correlator = dc,
            )
            Tracking.downconvert_and_correlate!(
                dc,
                measurements,
                ts;
                chunk_index = 0,
                chunk_duration = 1e-3s,
                stop_before_partial = true,
            )
            @inferred Tracking.estimate_dopplers_and_filter_prompt!(ts, measurements)
            allocated
        end
        @test allocations[1] == allocations[2]
    end
end

@testset "the loop bandwidths are settable per satellite" begin
    # Seeded per satellite and preserved across a reset — and set through one
    # keyword-based updater rather than a function per field, so an omitted
    # keyword leaves that bandwidth as it was.
    for state_of in (externally_fed_state, vt_externally_fed_state)
        ts = state_of()
        set_loop_filter_bandwidths!(ts, 1, 1; carrier = 5.0Hz, code = 0.5Hz)
        state = Tracking.get_doppler_estimator_state(get_sat_state(ts, 1))
        @test state.carrier_loop_filter_bandwidth == 5.0Hz
        @test state.code_loop_filter_bandwidth == 0.5Hz
        # Integer hertz is accepted and floated, like every other dimensioned
        # setter in this package; and the omitted `code` is left alone.
        set_loop_filter_bandwidths!(ts, 1; carrier = 12Hz)
        state = Tracking.get_doppler_estimator_state(get_sat_state(ts, 1))
        @test state.carrier_loop_filter_bandwidth === 12.0Hz
        @test state.code_loop_filter_bandwidth == 0.5Hz
        # The override is what a reset preserves — the reason it lives on the
        # per-satellite state at all.
        Tracking.reset_loop_filters!(ts)
        @test Tracking.get_doppler_estimator_state(get_sat_state(ts, 1)).code_loop_filter_bandwidth ==
              0.5Hz
    end

    # A bandwidth is a frequency and carries its unit; both wrong shapes get a
    # sentence rather than a MethodError.
    ts = externally_fed_state()
    @test_throws ArgumentError set_loop_filter_bandwidths!(ts, 1; carrier = 5.0)
    @test_throws ArgumentError set_loop_filter_bandwidths!(ts, 1; code = 5.0s)
end

end
