module MultiSignalCombiningTest

# Tests for multi-signal discriminator combining: the passenger signals'
# (`signals[2:end]`) discriminator outputs are folded into the estimator-driver
# signal's before the loop filters see them.
#
# The four properties worth pinning, in order of how badly a regression would
# hurt:
#
#  1. **Bit-identity where nothing is combined.** A single-signal satellite, and
#     a multi-signal one with `signal_combining = false`, must produce exactly
#     the Dopplers they produced before combining existed. The whole feature is
#     safe-by-default only because of this.
#  2. **The weighted mean is a mean.** Two signals carrying the same information
#     must not double the loop gain — a plain sum would.
#  3. **The inter-signal correction is subtracted with the right sign**, and is
#     converted to chips at the code frequency in effect.
#  4. **An unset code bias gates the code loop only** — the carrier loops combine
#     from the first integration, with nothing supplied by the caller.

using Test: @allocated, @inferred, @test, @testset, @test_throws
using Unitful: Hz, NoUnits, s, uconvert
using GNSSSignals:
    get_relative_power,
    GPSL1CA,
    GPSL1C_D,
    GPSL1C_P,
    GalileoE1B,
    GalileoE1C,
    GPSL5I,
    GPSL5Q,
    GalileoE5aI,
    GalileoE5aQ,
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
    NumAnts,
    TrackState,
    TrackedSat,
    TrackedSignal,
    VeryEarlyPromptLateCorrelator,
    add_satellite!,
    dll_disc_noise_gain,
    get_carrier_doppler,
    get_code_doppler,
    get_code_bias,
    get_sat_state,
    set_code_bias!,
    track
using Tracking:
    PassengerFoldContext,
    _accumulate_passenger_discriminators,
    _combine_with_own,
    _nominal_record_snr

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
function identical_l1ca_track_state(n, prn, carrier_doppler; kwargs...)
    gpsl1 = GPSL1CA()
    # Combining is opt-in, so every fixture that means to exercise it says so.
    estimator = ConventionalAssistedPLLAndDLL(; signal_combining = true, kwargs...)
    tracked_signals = ntuple(
        _ -> TrackedSignal(
            gpsl1;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            # Two copies of one signal have no differential group delay between
            # them by construction, so pin 0.0 and the code loop combines. The
            # `nothing` default — and what it gates — is exercised separately
            # below.
            code_bias = 0.0,
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
    TrackState(gpsl1, sat; doppler_estimator = estimator)
end

@testset "single-signal sats are bit-identical with combining on" begin
    # A satellite with no passengers has nothing to combine, so with combining
    # switched on it must track *exactly* as with it off — not merely to within a
    # tolerance. `_process_signals` dispatches on the empty passenger tuple to
    # guarantee this at the type level.
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 20_000)

    ts_on = TrackState(;
        signal = GPSL1CA(),
        doppler_estimator = ConventionalAssistedPLLAndDLL(; signal_combining = true),
    )
    ts_on = add_satellite!(ts_on; prn, code_phase = 0.0, carrier_doppler)
    ts_off = TrackState(;
        signal = GPSL1CA(),
        doppler_estimator = ConventionalAssistedPLLAndDLL(; signal_combining = false),
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
    # discriminators equal the driver's. A weight-normalized mean of equal values
    # is that value, so the Dopplers must match the driver-only run closely; a
    # summing combiner would double the loop gain and visibly diverge.
    prn, carrier_doppler = 1, 1000.0Hz
    buf = l1ca_signal(prn, carrier_doppler, 0.0, 40_000)

    combined = track(copy(buf), identical_l1ca_track_state(2, prn, carrier_doppler), FS)
    driver_only = track(
        copy(buf),
        identical_l1ca_track_state(2, prn, carrier_doppler; signal_combining = false),
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

@testset "passenger records are consumed even with no driver record in the chunk" begin
    # The case that makes the cross-chunk accumulator necessary: the driver has
    # the *longest* primary code period, so most chunks (one shortest-code-period
    # each) hold passenger records and no driver record. Nothing may be left
    # behind in a passenger's `correlator_outputs` — the estimator owns clearing
    # them, and a stale record would be folded twice.
    gpsl1, l1cp = GPSL1CA(), GPSL1C_P()
    prn, carrier_doppler = 1, 1000.0Hz
    # Driver = L1C-P (10 ms primary period), passenger = L1 C/A (1 ms). The
    # default chunk is 1 ms, so nine chunks in ten see a C/A record only.
    ts = TrackState(;
        signals = (gps_l1 = (l1cp, gpsl1),),
        doppler_estimator = ConventionalAssistedPLLAndDLL(; signal_combining = true),
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

    # Those five records' discriminators are parked for the next chunk's driver
    # update rather than discarded.
    pending = Tracking.get_doppler_estimator_state(sat).pending_discriminators
    @test pending.pll_weight > 0
    # …but only the carrier ones: no code bias has been supplied for the C/A
    # passenger, so the code loop abstains. The carrier loop needs nothing from
    # the caller and contributes from the first integration — the whole point of
    # gating the two separately.
    @test pending.dll_weight == 0.0

    # Supply the passenger's bias and the code contribution appears, with no
    # other change. Only 3 more ms, so the L1C-P driver still has not completed a
    # record and the accumulator is observable — a driver record would consume and
    # zero it, which is the whole point of it but would make this assertion
    # vacuous. Only the passenger needs a value: the driver's own bias is 0.0 by
    # definition.
    set_code_bias!(ts, :gps_l1, prn, GPSL1CA, 2.0e-9)
    ts = track(l1ca_signal(prn, carrier_doppler, 0.0, 45_000, fs), ts, fs)
    sat = get_sat_state(ts, prn)
    @test isempty(Tracking.get_filtered_prompts(sat.signals[1]))
    pending = Tracking.get_doppler_estimator_state(sat).pending_discriminators
    @test pending.dll_weight > 0
end

@testset "code bias: default, plumbing, and addressing" begin
    # Every signal starts unknown, whatever its constellation: this package holds
    # no table of per-signal values, because deciding one needs constellation
    # knowledge and usually a decoded navigation message. A signal defaulting to
    # `0.0` here would be this package guessing on the caller's behalf, in the
    # one direction that is not safe — see `set_code_bias!`.
    for signal in (
        GPSL1CA(),
        GPSL1C_D(),
        GPSL1C_P(),
        GPSL5I(),
        GPSL5Q(),
        GalileoE1B(),
        GalileoE1C(),
        GalileoE5aI(),
        GalileoE5aQ(),
    )
        @test get_code_bias(TrackedSignal(signal)) === nothing
    end

    l1cp, l1cd, l1ca = GPSL1C_P(), GPSL1C_D(), GPSL1CA()
    estimator = ConventionalAssistedPLLAndDLL()
    ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd, l1ca),))
    ts = add_satellite!(
        ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )

    # Addressed by signal type … (L1C-D is `signals[2]`; the driver is excluded
    # from carrying a bias at all — see below.)
    set_code_bias!(ts, :gps_l1, 7, GPSL1C_D, -1.5e-9)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === -1.5e-9
    # … and by index, and the other signals are untouched.
    set_code_bias!(ts, :gps_l1, 7, 3, 3.0e-9)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), 3) === 3.0e-9
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), GPSL1C_P) === nothing

    # `nothing` marks it unknown again — the one case a plain `Maybe` kwarg
    # could not express, hence the `Some` wrapper inside.
    set_code_bias!(ts, :gps_l1, 7, GPSL1C_D, nothing)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === nothing

    # An integer is accepted for a value in seconds and lands as a Float64.
    set_code_bias!(ts, :gps_l1, 7, GPSL1CA, 0)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), GPSL1CA) === 0.0
end

@testset "a non-zero code bias on the driver is refused, not ignored" begin
    # A bias is a differential *against the driver*, so the driver's own is 0 by
    # definition and the combining fold never reads it. Accepting a non-zero one
    # silently is the dangerous direction: it is what a caller does when it
    # derives every signal's group delay against some external datum and sets
    # them all, and then every passenger's bias is short by the driver's — a
    # constant, invisible offset in the shared code phase, which is what the bias
    # exists to prevent in the first place.
    l1cp, l1cd, l1ca = GPSL1C_P(), GPSL1C_D(), GPSL1CA()
    ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd, l1ca),))
    ts = add_satellite!(
        ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )

    @test_throws ArgumentError set_code_bias!(ts, :gps_l1, 7, GPSL1C_P, 1.2e-9)
    @test_throws ArgumentError set_code_bias!(ts, :gps_l1, 7, 1, -1.2e-9)

    # `0.0` and `nothing` restate the definition rather than contradicting it, so
    # a caller that writes a bias for every signal of the group still works as
    # long as the driver's is zero.
    set_code_bias!(ts, :gps_l1, 7, GPSL1C_P, 0.0)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), GPSL1C_P) === 0.0
    set_code_bias!(ts, :gps_l1, 7, 1, nothing)
    @test get_code_bias(get_sat_state(ts, :gps_l1, 7), 1) === nothing

    # A single-signal satellite is its own driver, so the same rule applies to
    # the selector-free form.
    single = add_satellite!(
        TrackState(; signal = GPSL1CA());
        prn = 3,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )
    @test_throws ArgumentError set_code_bias!(single, 3, 1.2e-9)
end

# Build one completed correlator record plus the pieces
# `_accumulate_passenger_discriminators` needs, for a passenger whose correlator
# sits at a known code error. Using the accumulator directly is what lets the
# sign of the group-delay correction be pinned exactly, rather than inferred
# from a converged loop.
function passenger_contribution(
    tracked_signal,
    correlator,
    integrated_samples,
    passenger_bias;
    code_doppler = 0.0Hz,
    # Default: driver and passenger co-phased, so the PLL de-rotation is the
    # identity and the accumulated phase equals the bare discriminator. Not
    # `0.0` — GPS L1 C/A's own carrier phase offset is −π/2, so a literal zero
    # here would inject a 90° rotation and mean "the driver is an L1C component".
    driver_carrier_phase_offset = get_carrier_phase_offset(
        Tracking.get_signal(tracked_signal),
    ),
    previous_prompt = complex(0.0, 0.0),
    correlated_pre_sync = false,
)
    signal = Tracking.get_signal(tracked_signal)
    ts = TrackedSignal(signal; correlator, code_bias = passenger_bias)
    output = CorrelatorOutput(correlator, integrated_samples, integrated_samples)
    ctx = PassengerFoldContext(1, FS, code_doppler, Float64(driver_carrier_phase_offset))
    _accumulate_passenger_discriminators(
        zero(DiscriminatorAccumulator),
        ts,
        output,
        correlator,
        previous_prompt,
        ctx,
        correlated_pre_sync,
    )
end

@testset "the code bias is subtracted, in chips, with the right sign" begin
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
    snr = _nominal_record_snr(gpsl1, n)

    # A zero bias asserts "this signal's code phase *is* the driver's", so the
    # discriminator enters as measured, weighted by SNR over the noise gain.
    same = passenger_contribution(base, correlator, n, 0.0)
    @test same.dll_weight ≈ snr / gain
    @test same.dll_sum / same.dll_weight ≈ raw

    # A positive bias means the passenger sits at a *larger* code phase than the
    # driver, so its raw discriminator over-reads the shared error by
    # `bias · f_code` chips — which must be subtracted.
    δt = 4.0e-9                      # 4 ns of extra passenger group delay
    # Computed independently of the implementation, from seconds × chips/second.
    δchips = uconvert(NoUnits, δt * 1.0s * get_code_frequency(gpsl1))
    later = passenger_contribution(base, correlator, n, δt)
    @test later.dll_sum / later.dll_weight ≈ raw - δchips
    # …and the opposite sign for a passenger that leads the driver.
    earlier = passenger_contribution(base, correlator, n, -δt)
    @test earlier.dll_sum / earlier.dll_weight ≈ raw + δchips

    # Converted at the code frequency actually in effect, code Doppler included —
    # not at the nominal chip rate.
    code_doppler = 1.0Hz * 4000
    doppler_shifted = passenger_contribution(base, correlator, n, δt; code_doppler)
    raw_shifted = Tracking.dll_disc(gpsl1, correlator, code_doppler, FS)
    @test doppler_shifted.dll_sum / doppler_shifted.dll_weight ≈
          raw_shifted -
          uconvert(NoUnits, δt * 1.0s * (get_code_frequency(gpsl1) + code_doppler))
end

@testset "an unset code bias gates the code loop only" begin
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    base = TrackedSignal(gpsl1; correlator)
    raw = Tracking.dll_disc(gpsl1, correlator, 0.0Hz, FS)

    # No bias supplied: the passenger stays out of the code loop …
    unset = passenger_contribution(base, correlator, n, nothing)
    @test unset.dll_weight == 0.0
    # … but its carrier phase contribution is unaffected: carrier phase carries
    # no inter-signal code bias, so it combines from the first integration with
    # nothing supplied at all.
    @test unset.pll_weight > 0
    @test unset.pll_sum / unset.pll_weight ≈ Tracking.pll_disc(gpsl1, correlator)

    # `0.0` is a different statement from `nothing`, and the only one of the two
    # that lets the code loop combine: it asserts the passenger's code phase is
    # the driver's, rather than leaving it unknown.
    zeroed = passenger_contribution(base, correlator, n, 0.0)
    @test zeroed.dll_weight > 0
    @test zeroed.dll_sum / zeroed.dll_weight ≈ raw
end

@testset "weights: SNR, discriminator noise gain, and the FLL's T²" begin
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    base = TrackedSignal(gpsl1; correlator)

    # SNR is `|P|² · N`: the sample-normalized prompt's magnitude does not grow
    # with the record length, but its noise variance falls as 1/N.
    # Nominal SNR is `power share × N` — the prompt does not enter it at all.
    @test _nominal_record_snr(gpsl1, 100) === get_relative_power(gpsl1) * 100
    @test _nominal_record_snr(gpsl1, 200) === get_relative_power(gpsl1) * 200
    # GPS L1C's 75/25 pilot/data split is the one asymmetric case, and it is
    # exactly the ratio a joint L1C loop needs.
    @test get_relative_power(GPSL1C_P()) / get_relative_power(GPSL1C_D()) === 3.0
    # Galileo's intra-band pairs split evenly, so combining E1B+E1C (or E5aI+Q)
    # at equal integration length is a plain average of two independent
    # discriminators — a √2 noise reduction.
    @test get_relative_power(GalileoE1B()) === get_relative_power(GalileoE1C()) === 0.5
    @test get_relative_power(GalileoE5aI()) === get_relative_power(GalileoE5aQ()) === 0.5

    # Code-loop weight scales with the record length …
    short = passenger_contribution(base, correlator, 5000, 0.0)
    long = passenger_contribution(base, correlator, 50_000, 0.0)
    @test long.dll_weight ≈ 10 * short.dll_weight

    # A zeroed previous prompt means "no frequency measurement yet", not "zero
    # frequency error", so it must not pull the combined FLL towards 0 Hz.
    @test short.fll_weight == 0.0
    # With a previous prompt the FLL weight carries the extra `N²`: dividing the
    # inter-prompt phase difference by 2π·T scales its noise by 1/T.
    with_prev = passenger_contribution(
        base,
        correlator,
        5000,
        0.0;
        previous_prompt = complex(0.9, 0.1),
    )
    with_prev_long = passenger_contribution(
        base,
        correlator,
        50_000,
        0.0;
        previous_prompt = complex(0.9, 0.1),
    )
    @test with_prev.fll_weight ≈ _nominal_record_snr(gpsl1, 5000) * 5000^2
    @test with_prev_long.fll_weight ≈ 1000 * with_prev.fll_weight

    # A BOC VEML discriminator is several times more accurate than a 1-chip BPSK
    # early-late one at equal SNR, so its noise gain must be correspondingly
    # smaller — that is what lets a pilot outvote a legacy signal.
    veml = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    epl = EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    g_veml = dll_disc_noise_gain(GPSL1C_P(), veml, 0.0Hz, FS)
    g_epl = dll_disc_noise_gain(GPSL1CA(), epl, 0.0Hz, FS)
    @test 0 < g_veml < g_epl
    @test 4 < g_epl / g_veml < 16

    # An unmodelled correlator type falls back to the 1-chip early-late gain
    # rather than erroring — mis-weighting costs efficiency, never bias.
    @test dll_disc_noise_gain(GPSL1CA(), veml, 0.0Hz, FS) !==
          dll_disc_noise_gain(GPSL1CA(), epl, 0.0Hz, FS)
end

@testset "_combine_with_own is a normalized mean with an exact driver-only path" begin
    # The bit-identity guarantee of the whole feature reduces to this branch.
    @test _combine_with_own(0.0, 0.0, 0.37, 12.5) === 0.37
    # Passenger-only (driver record carries no usable measurement of its own).
    @test _combine_with_own(6.0, 3.0, 99.0, 0.0) === 2.0
    # Equal values combine to that value regardless of the weights — the
    # property that keeps the loop gain fixed as contributors come and go.
    @test _combine_with_own(5.0 * 2.0, 5.0, 2.0, 3.0) ≈ 2.0
    # And a genuine weighted mean otherwise.
    @test _combine_with_own(4.0 * 1.0, 4.0, 3.0, 4.0) ≈ 2.0
end

@testset "quadrature passengers are de-rotated before the PLL" begin
    # GPS L5 is a quadrature pair: with the pilot driving, the data component
    # sits on the imaginary axis. Without de-rotation `atan(imag/real)` reads
    # ±π/2 — a fabricated phase error that would drag the combined PLL off lock.
    l5i, l5q = GPSL5I(), GPSL5Q()
    driver_phase = get_carrier_phase_offset(l5q)
    # A data prompt as the driver's real-axis lock leaves it: rotated by the
    # component's own carrier phase offset relative to the driver.
    on_axis = complex(1.0, 0.0)
    quadrature = on_axis * cis(get_carrier_phase_offset(l5i) - driver_phase)
    correlator =
        EarlyPromptLateCorrelator([complex(0.7, 0.0), quadrature, complex(0.7, 0.0)], 0.5)

    rotated = Tracking.pll_disc(
        l5i,
        correlator,
        Tracking._carrier_phase_derotation(driver_phase, l5i),
    )
    @test abs(rotated) < 1e-12
    @test abs(Tracking.pll_disc(l5i, correlator)) > 1.0     # the un-rotated trap

    # The three-argument method is bit-identical to the two-argument one for a
    # co-phased component, so the driver itself is untouched.
    @test Tracking.pll_disc(l5q, correlator, complex(1.0, 0.0)) ===
          Tracking.pll_disc(l5q, correlator)
end

@testset "only a pre-sync-invalidated record is excluded" begin
    # Nominal weights say what a component's power share should be; they cannot
    # tell whether it is being received, and there is deliberately no check that
    # they agree with the measured prompt — a signal is in a satellite's tuple
    # only because the caller declared it. So a record contributes on its nominal
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

    # A prompt 40 dB down still contributes, and with exactly the same weight as
    # a strong one: the weight is nominal, so the prompt magnitude never enters it.
    a = passenger_contribution(base, strong, n, 0.0)
    b = passenger_contribution(base, weak, n, 0.0)
    @test a.pll_weight === b.pll_weight === _nominal_record_snr(gpsl1, n)
    @test a.dll_weight ≈ b.dll_weight

    # The one exclusion: a record correlated with a pre-sync replica on a
    # secondary-coded signal, whose coherent sum partially cancels without the
    # secondary-code wipe-off, so its discriminators measure nothing.
    e1c = TrackedSignal(GalileoE1C(); correlator = strong)
    @test passenger_contribution(e1c, strong, n, 0.0; correlated_pre_sync = true) ==
          zero(DiscriminatorAccumulator)

    # A signal with no secondary code correlates identically either side of the
    # sync instant, so it is unaffected.
    @test passenger_contribution(base, strong, n, 0.0; correlated_pre_sync = true).pll_weight >
          0
end

@testset "the combining fold is inferable and adds no allocations" begin
    # The interleaved fold threads a tuple of signals, a tuple of cursors and a
    # tuple of chunk-boundary sync flags through a recursive walk. That has to
    # fold at inference time like the rest of the hot path, and must not put the
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
                signal_combining = combining,
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
        # Combining is free: identical allocation to the driver-only path, so the
        # residue is `track!`'s pre-existing per-call overhead and nothing the
        # fold added.
        @test allocations[1] == allocations[2]
    end
end

@testset "reset_loop_filters! drops pending passenger discriminators" begin
    # Pending measurements belong to the pre-reset cadence and Doppler; carrying
    # them past a deliberate loop reset is exactly the history the reset exists
    # to discard.
    prn, carrier_doppler = 1, 1000.0Hz
    ts = identical_l1ca_track_state(2, prn, carrier_doppler)
    ts = track(l1ca_signal(prn, carrier_doppler, 0.0, 20_000), ts, FS)
    Tracking.reset_loop_filters!(ts)
    pending =
        Tracking.get_doppler_estimator_state(get_sat_state(ts, prn)).pending_discriminators
    @test pending == zero(DiscriminatorAccumulator)
end

end
