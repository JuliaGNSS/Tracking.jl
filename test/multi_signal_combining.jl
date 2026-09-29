module MultiSignalCombiningTest

# Tests for multi-signal discriminator combining: the passenger signals'
# (`signals[2:end]`) discriminator outputs are folded into the estimator-driver
# signal's before the loop filters see them.
#
# The four properties worth pinning, in order of how badly a regression would
# hurt:
#
#  1. **Bit-identity where nothing is combined.** A single-signal satellite, and
#     a multi-signal one with `discriminator_combining = false`, must produce exactly
#     the Dopplers they produced before combining existed. The whole feature is
#     safe-by-default only because of this.
#  2. **The weighted mean is a mean.** Two signals carrying the same information
#     must not double the loop gain — a plain sum would.
#  3. **The inter-signal correction is subtracted with the right sign**, and is
#     converted to chips at the code frequency in effect.
#  4. **An unset differential group delay gates the code loop only** — the carrier loops combine
#     from the first integration, with nothing supplied by the caller.

using Test: @allocated, @inferred, @test, @test_logs, @testset, @test_throws
using Unitful: @u_str, DimensionError, Hz, NoUnits, s, uconvert
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
    NoCN0Estimator,
    NumAnts,
    TrackState,
    TrackedSat,
    TrackedSignal,
    VectorPLLAndDLL,
    VeryEarlyPromptLateCorrelator,
    add_satellite!,
    append_correlator_output!,
    dll_disc_noise_gain,
    enable_vt!,
    estimate_dopplers_and_filter_prompt!,
    get_carrier_doppler,
    get_code_doppler,
    get_group_delay,
    get_preferred_num_code_blocks_to_integrate,
    get_sat_state,
    set_group_delay!,
    set_preferred_num_code_blocks_to_integrate!,
    track
using Tracking:
    RecordFoldContext, _accumulate_record, _normalized_mean, _relative_record_snr

const FS = 5e6Hz

# A correlator type with no `dll_disc_noise_gain` method: combining refuses it rather
# than guessing a gain.
struct UnmodelledCorrelator <: Tracking.AbstractCorrelator{1} end

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
    # Combining is opt-in, so every fixture that means to exercise it says so.
    estimator = ConventionalAssistedPLLAndDLL(; kwargs...)
    tracked_signals = ntuple(
        _ -> TrackedSignal(
            gpsl1;
            correlator = EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),
            # Two copies of one signal have no differential group delay between
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
    # tolerance. `_seed_and_fold` branches on the empty passenger tuple's type to
    # guarantee this at compile time.
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
    # discriminators equal the driver's. A weight-normalized mean of equal values
    # is that value, so the Dopplers must match the driver-only run closely; a
    # summing combiner would double the loop gain and visibly diverge.
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

    # Those five records' discriminators are parked for the next chunk's driver
    # update rather than discarded.
    pending = Tracking.get_doppler_estimator_state(sat).pending_combining_sums
    @test pending.pll_weight > 0
    # …but only the carrier ones: no group delay has been supplied for the C/A
    # passenger, so the code loop abstains. The carrier loop needs nothing from
    # the caller and contributes from the first integration — the whole point of
    # gating the two separately.
    @test pending.dll_weight == 0.0

    # Supply the passenger's delay and the code contribution appears, with no
    # other change. Only 3 more ms, so the L1C-P driver still has not completed a
    # record and the accumulator is observable — a driver record would consume and
    # zero it, which is the whole point of it but would make this assertion
    # vacuous. Both ends are needed — the difference is what the loop applies, so
    # the datum has to be stated as much as the passenger's own value.
    set_group_delay!(ts, :gps_l1, prn, 1, 0.0s)
    set_group_delay!(ts, :gps_l1, prn, 1, 0.0s)
    set_group_delay!(ts, :gps_l1, prn, GPSL1CA, 2.0e-9s)
    ts = track(l1ca_signal(prn, carrier_doppler, 0.0, 45_000, fs), ts, fs)
    sat = get_sat_state(ts, prn)
    @test isempty(Tracking.get_filtered_prompts(sat.signals[1]))
    pending = Tracking.get_doppler_estimator_state(sat).pending_combining_sums
    @test pending.dll_weight > 0
end

@testset "differential group delay: default, plumbing, and addressing" begin
    # Every signal starts unknown, whatever its constellation: this package holds
    # no table of per-signal values, because deciding one needs constellation
    # knowledge and usually a decoded navigation message. A signal defaulting to
    # `0.0` here would be this package guessing on the caller's behalf, in the
    # one direction that is not safe — see `set_group_delay!`.
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
        @test get_group_delay(TrackedSignal(signal)) === nothing
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

    # Addressed by signal type … (L1C-D is `signals[2]`; the driver is not
    # settable at all — see below.)
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, -1.5e-9s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === -1.5e-9s
    # … and readable straight off the `TrackState`, like every other per-signal
    # setting: a value that can be set through one addressing form and only read
    # back through another is the kind of asymmetry nothing but a test notices.
    @test get_group_delay(ts, :gps_l1, 7, GPSL1C_D) === -1.5e-9s
    # … and by index, and the other signals are untouched.
    set_group_delay!(ts, :gps_l1, 7, 3, 3.0e-9s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), 3) === 3.0e-9s
    # Slot 1 is an ordinary slot: untouched by the writes above, it is still the
    # `nothing` every slot starts at, and it is settable like any other.
    @test isnothing(get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_P))
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_P, 0.0s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_P) === 0.0s

    # `nothing` marks it unknown again — the one case a plain `Maybe` kwarg
    # could not express, hence the `Some` wrapper inside.
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, nothing)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1C_D) === nothing

    # Any time unit is accepted and normalized to seconds as a Float64 quantity,
    # so the field type stays concrete however the caller spells the value.
    set_group_delay!(ts, :gps_l1, 7, GPSL1CA, -1200u"ps")
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1CA) === -1.2e-9s
    set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 0s)
    @test get_group_delay(get_sat_state(ts, :gps_l1, 7), GPSL1CA) === 0.0s

    # A bare number, and a quantity that is not a time, are both refused: the
    # field carries its unit like every dimensioned quantity here, and at
    # sub-nanosecond scale an assumed unit is a metre-scale mistake. `convert`
    # does the refusing — one conversion point, no bespoke messages.
    @test_throws DimensionError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 1.2e-9)
    @test_throws DimensionError TrackedSignal(GPSL1CA(); group_delay = 1.2e-9)
    @test_throws DimensionError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 1.2Hz)
    # A length is the *likely* wrong input rather than an exotic one: an SSR/PPP
    # "code bias" is a pseudorange correction in metres.
    # …except a length, which gets a sentence naming the conversion and the datum
    # trap that comes with it: an SSR/PPP "code bias" in metres is the likely
    # wrong input here, not an exotic one.
    @test_throws ArgumentError set_group_delay!(ts, :gps_l1, 7, GPSL1CA, 0.45u"m")
    @test_throws ArgumentError TrackedSignal(GPSL1CA(); group_delay = 0.45u"m")
    @test occursin("Unitful.c0", sprint(showerror, try
        Tracking._as_group_delay(0.45u"m")
    catch e
        e
    end))
end

@testset "every slot starts unknown, and the datum is whichever one you state" begin
    # Only *differences* between a satellite's signals are used, so the datum is
    # the caller's to choose and `signals[1]` is an ordinary slot: it starts at
    # `nothing` like every other, and it is written like every other.
    l1cp, l1cd, l1ca = GPSL1C_P(), GPSL1C_D(), GPSL1CA()
    ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd, l1ca),))
    ts = add_satellite!(
        ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )
    @test isnothing(get_group_delay(ts, :gps_l1, 7, 1))

    # Every slot is settable, the driver's included, and reads back what it was
    # given — the stored value is the caller's datum, not a referred one.
    set_group_delay!(ts, :gps_l1, 7, GPSL1C_D, 1.2e-9s)
    set_group_delay!(ts, :gps_l1, 7, 1, 0.0s)
    @test get_group_delay(ts, :gps_l1, 7, 1) === 0.0s
    @test get_group_delay(ts, :gps_l1, 7, GPSL1C_D) == 1.2e-9s

    # Shifting every value by one constant changes no difference, which is the
    # whole property — and what the fold applies is the difference.
    offset = 5.0e-9s
    shifted = TrackState(; signals = (gps_l1 = (l1cp, l1cd, l1ca),))
    shifted = add_satellite!(
        shifted;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )
    set_group_delay!(shifted, :gps_l1, 7, 1, offset)
    set_group_delay!(shifted, :gps_l1, 7, GPSL1C_D, 1.2e-9s + offset)
    @test get_group_delay(shifted, :gps_l1, 7, 1) == offset

    # The two satellites' passengers carry the same *difference*, so the chips
    # each context hands the fold are equal — the datum has cancelled.
    chips(t) = Tracking.RecordFoldContext(
        get_sat_state(t, :gps_l1, 7).signals[2],
        7,
        FS,
        0.0Hz,
        0.0,
        get_group_delay(t, :gps_l1, 7, 1),
        Tracking.NoiseFloor(nothing, true),
    ).group_delay_chips
    @test chips(ts) ≈ chips(shifted)

    # Clearing either end withholds that passenger's code contribution.
    set_group_delay!(shifted, :gps_l1, 7, 1, nothing)
    @test isnothing(chips(shifted))
    @test isnothing(
        Tracking.RecordFoldContext(
            get_sat_state(ts, :gps_l1, 7).signals[3],
            7,
            FS,
            0.0Hz,
            0.0,
            get_group_delay(ts, :gps_l1, 7, 1),
            Tracking.NoiseFloor(nothing, true),
        ).group_delay_chips,
    )
end

@testset "a pre-built signals tuple states its own datum" begin
    # The `TrackedSignal`-tuple constructor is the full-control form and seeds
    # nothing, so a caller using it owns the datum as it owns the correlators.
    l1cp, l1ca = GPSL1C_P(), GPSL1CA()
    sat = TrackedSat((TrackedSignal(l1cp), TrackedSignal(l1ca)), 1, 0.0, 1000.0Hz)
    @test isnothing(get_group_delay(sat, 1))
    @test isnothing(get_group_delay(sat, 2))

    # Slot 1 is an ordinary settable slot, and a rebuild may put a datum value on
    # it — under the old driver-must-be-zero rule this was an error to catch.
    sat = Tracking._set_sat_group_delay(sat, 2.0e-9s, 2)
    @test get_group_delay(sat, 2) == 2.0e-9s
    rebuilt = TrackedSat(
        sat;
        # `Some` because `nothing` is a legal value of this field, so the
        # copy-update constructor cannot read a bare `nothing` as "clear it".
        signals = (
            TrackedSignal(
                sat.signals[1];
                group_delay = Some{Union{Nothing,typeof(1.0s)}}(1.0e-9s),
            ),
            sat.signals[2],
        ),
    )
    @test get_group_delay(rebuilt, 1) == 1.0e-9s
end

@testset "a value on slot 1 is the datum, at construction as anywhere else" begin
    # No normalisation and nothing refused: `signals[1]` carries a datum value like
    # any other slot, and the constructors seed neither.
    l1cp, l1cd = GPSL1C_P(), GPSL1C_D()
    sat = TrackedSat(
        (
            TrackedSignal(l1cp; group_delay = 1.0e-9s),
            TrackedSignal(l1cd; group_delay = -1.0e-9s),
        ),
        7,
        0.0,
        1000.0Hz,
    )
    @test get_group_delay(sat, GPSL1C_P) == 1.0e-9s
    @test get_group_delay(sat, GPSL1C_D) == -1.0e-9s

    fresh = TrackedSat((l1cp, l1cd), 7, 0.0, 1000.0Hz)
    @test isnothing(get_group_delay(fresh, GPSL1C_P))
    @test isnothing(get_group_delay(fresh, GPSL1C_D))
end

@testset "setting a delay under vector tracking is stored and read back" begin
    # No warning here, by design: the setter is about the satellite's signals,
    # not about which estimator will read the value, and `VectorPLLAndDLL`
    # leaving it unread is a documentation matter — see `set_group_delay!`.
    l1cp, l1cd = GPSL1C_P(), GPSL1C_D()
    vector_ts = TrackState(;
        signals = (gps_l1 = (l1cp, l1cd),),
        doppler_estimator = Tracking.VectorPLLAndDLL(),
    )
    vector_ts = add_satellite!(
        vector_ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )
    @test_logs set_group_delay!(vector_ts, :gps_l1, 7, GPSL1C_D, 1.0e-9s)
    @test get_group_delay(vector_ts, :gps_l1, 7, GPSL1C_D) == 1.0e-9s
    # The conventional estimator behaves identically.
    scalar_ts = TrackState(; signals = (gps_l1 = (l1cp, l1cd),))
    scalar_ts = add_satellite!(
        scalar_ts;
        group = :gps_l1,
        prn = 7,
        code_phase = 0.0,
        carrier_doppler = 1000.0Hz,
    )
    @test_logs set_group_delay!(scalar_ts, :gps_l1, 7, GPSL1C_D, 1.0e-9s)
end

# Build one completed correlator record plus the pieces
# `_accumulate_record` needs, for a passenger whose correlator
# sits at a known code error. Using the accumulator directly is what lets the
# sign of the group-delay correction be pinned exactly, rather than inferred
# from a converged loop.
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
    correlated_pre_sync = false,
    # The driver is the datum unless a test says otherwise, so the passenger's
    # stored value is directly its referred bias.
    driver_delay = 0.0s,
)
    signal = Tracking.get_signal(tracked_signal)
    ts = TrackedSignal(signal; correlator, group_delay = passenger_delay)
    output = CorrelatorOutput(correlator, integrated_samples, integrated_samples)
    # The per-passenger context derives its de-rotation and chunk-boundary sync
    # state from the signal itself, so build it through that constructor rather
    # than by listing fields — the point of the test is the accumulate step, not
    # the context's field order. No noise estimator here, hence a `nothing` floor.
    ctx = RecordFoldContext(
        ts,
        1,
        FS,
        code_doppler,
        Float64(driver_carrier_phase_offset),
        driver_delay,
        Tracking.NoiseFloor(nothing, true),
    )
    _accumulate_record(
        DiscriminatorAccumulator(),
        ts,
        output,
        correlator,
        previous_prompt,
        ctx,
        correlated_pre_sync,
    )
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
    # the shared error by `δ · f_code` chips, which the referral adds back. The
    # driver is the datum here (`driver_delay = 0.0s`), so the passenger's stored
    # value is the difference.
    δt = 4.0e-9s                     # 4 ns more group delay than the driver
    # Computed independently of the implementation, from seconds × chips/second.
    δchips = uconvert(NoUnits, δt * get_code_frequency(gpsl1))
    later = passenger_contribution(base, correlator, n, δt)
    @test later.dll_sum / later.dll_weight ≈ raw + δchips
    # …and the opposite sign for a passenger that leads the driver.
    earlier = passenger_contribution(base, correlator, n, -δt)
    @test earlier.dll_sum / earlier.dll_weight ≈ raw - δchips

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

    # `0.0` is a different statement from `nothing`, and the only one of the two
    # that lets the code loop combine: it asserts the passenger's code phase is
    # the driver's, rather than leaving it unknown.
    zeroed = passenger_contribution(base, correlator, n, 0.0s)
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
    @test _relative_record_snr(gpsl1, 100) === get_relative_power(gpsl1) * 100
    @test _relative_record_snr(gpsl1, 200) === get_relative_power(gpsl1) * 200
    # GPS L1C's 75/25 pilot/data split is the one asymmetric case, and it is
    # exactly the ratio a joint L1C loop needs.
    @test get_relative_power(GPSL1C_P()) / get_relative_power(GPSL1C_D()) === 3.0
    # Galileo's intra-band pairs split evenly, so combining E1B+E1C (or E5aI+Q)
    # at equal integration length is a plain average of two independent
    # discriminators — a √2 noise reduction.
    @test get_relative_power(GalileoE1B()) === get_relative_power(GalileoE1C()) === 0.5
    @test get_relative_power(GalileoE5aI()) === get_relative_power(GalileoE5aQ()) === 0.5

    # Code-loop weight scales with the record length …
    short = passenger_contribution(base, correlator, 5000, 0.0s)
    long = passenger_contribution(base, correlator, 50_000, 0.0s)
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

    # A BOC VEML discriminator is several times more accurate than a 1-chip BPSK
    # early-late one at equal SNR, so its noise gain must be correspondingly
    # smaller — that is what lets a pilot outvote a legacy signal.
    veml = VeryEarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    epl = EarlyPromptLateCorrelator(; num_ants = NumAnts(1))
    g_veml = dll_disc_noise_gain(GPSL1C_P(), veml, 0.0Hz, FS)
    g_epl = dll_disc_noise_gain(GPSL1CA(), epl, 0.0Hz, FS)
    @test 0 < g_veml < g_epl

    # The documented values, not just their ordering. `d / 4` is Kaplan & Hegarty's
    # tracking-jitter constant for the noncoherent early-minus-late envelope
    # discriminator with early-late spacing `d` chips, so it is pinned at the two
    # sampling frequencies where the preferred chip shift lands on a whole number
    # of samples and `d` is therefore exactly what it was asked for — which also
    # pins the linear `d` dependence the docstring claims.
    @test dll_disc_noise_gain(
        GPSL1CA(),
        EarlyPromptLateCorrelator(; num_ants = NumAnts(1)),   # ±0.5 chips
        0.0Hz,
        2.046e6Hz,                                            # 0.5 chips = 1 sample
    ) === 0.25
    @test dll_disc_noise_gain(
        GPSL1CA(),
        EarlyPromptLateCorrelator([complex(0.0), complex(0.0), complex(0.0)], 0.25),
        0.0Hz,
        4.092e6Hz,                                            # 0.25 chips = 1 sample
    ) === 0.125
    # …and the VEML gain for the default ±0.15/±0.6 taps, the `≈ 0.03` the
    # docstring quotes — about 8× the precision of the 1-chip early-late layout.
    @test g_veml ≈ 0.03125 rtol = 1e-12
    @test dll_disc_noise_gain(GPSL1C_P(), veml, 0.0Hz, 20e6Hz) ≈ 0.03125 rtol = 1e-9

    # A degenerate tap layout — both pairs past the 1-chip correlation support —
    # carries no delay information, and must come out as `Inf` so the record gets
    # weight *zero*. Not `NaN`: the envelope sum is what the S-curve slope divides
    # by, so `-0.0 / 0.0` is one line away, and a `NaN` weight propagates out of
    # the accumulator into the *driver's* code loop filter and wrecks a satellite
    # that was tracking fine (`iszero` compares equal to nothing).
    # Real accumulator values, because an all-zero correlator is a separate `0/0`
    # (nothing was correlated at all) and would confuse the two cases.
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
    @test dll_disc_noise_gain(GalileoE1C(), degenerate, 0.0Hz, 20e6Hz) === Inf
    # The paired guard in the discriminator itself: an uncalibrated raw value
    # rather than a division by a zero (or NaN) slope.
    @test isfinite(Tracking.dll_disc(GalileoE1C(), degenerate, 0.0Hz, 20e6Hz))
    # And such a record is simply excluded from the code loop rather than
    # poisoning it — `Inf` weight is what makes that happen.
    excluded = passenger_contribution(
        TrackedSignal(GalileoE1C(); correlator = degenerate),
        degenerate,
        5000,
        0.0s,
    )
    @test excluded.dll_weight === 0.0
    # Zero weight is zero contribution, and reads as a zero rather than `0/0`.
    @test _normalized_mean(excluded.dll_sum, excluded.dll_weight) === 0.0
    @test excluded.pll_weight > 0   # the carrier loops are unaffected

    # No generic gain: a correlator type without a `dll_disc_noise_gain` method is a
    # `MethodError`, like one without `dll_disc`.
    @test_throws MethodError dll_disc_noise_gain(
        GPSL1CA(),
        UnmodelledCorrelator(),
        0.0Hz,
        FS,
    )
end

@testset "_normalized_mean divides, and reads zero weight as zero" begin
    # A plain weighted mean.
    @test _normalized_mean(4.0 * 2.0, 4.0) === 2.0
    # Equal values give that value regardless of the weights — the property that
    # keeps the loop gain fixed as contributors come and go.
    @test _normalized_mean(5.0 * 2.0 + 3.0 * 2.0, 5.0 + 3.0) ≈ 2.0
    # Zero weight answers a zero of the sum's own type, never `0/0`. Not a
    # degenerate case: it is the FLL on a signal's first record, and the code
    # loop on a record whose group delay could not be referred.
    @test _normalized_mean(0.0, 0.0) === 0.0
    @test _normalized_mean(0.0Hz, 0.0) === 0.0Hz
    @test !isnan(_normalized_mean(0.0, 0.0))
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

    # The identity rotation leaves the driver's own discriminator as it was.
    @test Tracking.pll_disc(
        l5q,
        Tracking._derotate_correlator(correlator, complex(1.0, 0.0)),
    ) === Tracking.pll_disc(l5q, correlator)

    # …and the passenger fold really passes that rotation in. Asserting only on
    # the helper above would leave the *wiring* untested: swap
    # `_carrier_phase_derotation(…)` for `one(ComplexF64)` inside
    # `_accumulate_record` and every assertion so far still
    # passes, while the combined PLL reads a fabricated ±π/2. The delay is left
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
    a = passenger_contribution(base, strong, n, 0.0s)
    b = passenger_contribution(base, weak, n, 0.0s)
    @test a.pll_weight === b.pll_weight === _relative_record_snr(gpsl1, n)
    @test a.dll_weight ≈ b.dll_weight

    # The one exclusion: a record correlated with a pre-sync replica on a
    # secondary-coded signal, whose coherent sum partially cancels without the
    # secondary-code wipe-off, so its discriminators measure nothing.
    e1c = TrackedSignal(GalileoE1C(); correlator = strong)
    @test passenger_contribution(e1c, strong, n, 0.0s; correlated_pre_sync = true) ==
          DiscriminatorAccumulator()

    # A signal with no secondary code correlates identically either side of the
    # sync instant, so it is unaffected.
    @test passenger_contribution(base, strong, n, 0.0s; correlated_pre_sync = true).pll_weight >
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
        # Combining is free: identical allocation to the driver-only path, so the
        # residue is `track!`'s pre-existing per-call overhead and nothing the
        # fold added.
        @test allocations[1] == allocations[2]
    end

    # The same pin on a **heterogeneous** group, which is the realistic one and
    # the only one that can catch a regression in the tuple walk: three distinct
    # signal types, two distinct correlator types, and per-passenger contexts
    # whose noise-density type parameters differ. `n` copies of one signal make
    # every element of `signals`, `cursors` and `contexts` identically typed, so
    # a walk that accidentally widened to `Vector` or boxed a context would still
    # look free there.
    het_fs = 15e6Hz
    het_buf = l1ca_signal(prn, carrier_doppler, 0.0, 45_000, het_fs)
    het_measurements = (L1 = Tracking.BandMeasurement(het_buf, het_fs, 0.0Hz),)
    het_allocations = map((false, true)) do combining
        estimator = ConventionalAssistedPLLAndDLL()
        # An L1 C/A driver with L1C passengers is exactly the ordering
        # `ConventionalPLLAndDLL` warns against for a *converged* loop, but the
        # allocation contract is about the fold's types, not its numerics, and
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
                # Left unknown on purpose: the carrier-only passenger is the
                # branch whose weight types differ from its neighbours'.
                group_delay = nothing,
            ),
        )
        sat = TrackedSat(sigs, prn, 0.0, carrier_doppler; doppler_estimator = estimator)
        ts = TrackState(GPSL1CA(), sat; doppler_estimator = estimator)
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
    # record it folds. Nothing observable in a tracking run distinguishes the two
    # (the floors of co-tracked signals are similar, and no C/N₀ assertion is
    # made on them), so this checks the routing where it is decided.
    gpsl1 = GPSL1CA()
    correlator = EarlyPromptLateCorrelator(
        [complex(0.9, 0.0), complex(1.0, 0.0), complex(0.6, 0.0)],
        0.5,
    )
    n = 5000
    # Two steps: `correlator_outputs` is a kwarg of the copy-update constructor
    # only, so the record is planted on an already-built signal.
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
    noise = (
        Tracking.NoiseFloor(driver_density, true),
        Tracking.NoiseFloor(passenger_density, true),
    )
    new_driver, new_passengers, _, _, _ = Tracking._fold_satellite_chunk(
        sat.signals[1],
        Base.tail(sat.signals),
        sat,
        state,
        FS,
        noise,
        get_carrier_phase_offset(gpsl1),
        state.pending_combining_sums,
    )
    @test Tracking.get_cn0_estimator(new_driver).density ≈ driver_density
    @test Tracking.get_cn0_estimator(new_passengers[1]).density ≈ passenger_density
end

@testset "reset_loop_filters! drops pending passenger discriminators" begin
    # Pending measurements belong to the pre-reset cadence and Doppler; carrying
    # them past a deliberate loop reset is exactly the history the reset exists
    # to discard.
    #
    # This needs the long-driver fixture rather than two identical signals: with
    # equal code periods every passenger record is consumed — and the accumulator
    # zeroed — by the driver record that closes its window, so the accumulator is
    # *already* empty at the end of the chunk and asserting it is empty after the
    # reset would pass with the drop removed. A 10 ms L1C-P driver with a 1 ms
    # C/A passenger leaves records genuinely pending, which is what makes the
    # assertion mean something.
    gpsl1, l1cp = GPSL1CA(), GPSL1C_P()
    prn, carrier_doppler = 1, 1000.0Hz
    fs = 15e6Hz                                  # L1C-P TMBOC needs ≥ 12.276 MHz
    ts = TrackState(;
        signals = (gps_l1 = (l1cp, gpsl1),),
        doppler_estimator = ConventionalAssistedPLLAndDLL(),
        discriminator_combining = true,
    )
    ts = add_satellite!(ts; group = :gps_l1, prn, code_phase = 0.0, carrier_doppler)
    set_group_delay!(ts, :gps_l1, prn, 1, 0.0s)
    set_group_delay!(ts, :gps_l1, prn, GPSL1CA, 2.0e-9s)
    ts = track(l1ca_signal(prn, carrier_doppler, 0.0, 75_000, fs), ts, fs)   # 5 ms

    # Five C/A records are parked, in all three loops, with no driver record yet.
    pending_before =
        Tracking.get_doppler_estimator_state(get_sat_state(ts, prn)).pending_combining_sums
    @test pending_before.pll_weight > 0
    @test pending_before.dll_weight > 0
    @test pending_before.fll_weight > 0

    Tracking.reset_loop_filters!(ts)
    pending =
        Tracking.get_doppler_estimator_state(get_sat_state(ts, prn)).pending_combining_sums
    @test pending == DiscriminatorAccumulator()
end

# --------------------------------------------------------------------------
# The combined value actually reaches the loop filters
# --------------------------------------------------------------------------
#
# Every other testset in this file pins the arithmetic at the accumulator level,
# or compares a combined run against a driver-only one and expects them to
# *agree* — which two copies of one signal make true whether the combiner runs or
# not. So none of them would notice `_fold_satellite_chunk` pushing the
# driver's own discriminator into the filters instead of the combined one.
#
# These do, by feeding the two signals records that disagree and exploiting the
# fact that the loop filters and `aid_dopplers` are affine in the discriminator
# on a fresh filter state: if `d` is what the driver's record alone produces and
# `p` is what the passenger's record alone produces (installed on the driver, so
# the same code path evaluates it), then a weight-normalized mean of two
# equally-weighted records must land exactly halfway between them. No loop
# internals are reproduced here, which is what keeps the assertion honest.
#
# Records are appended through `append_correlator_output!` rather than correlated
# from samples, so each run is one record per signal at a known sample index and
# the arithmetic is exact rather than convergence-limited.

# Two records of one satellite that disagree about both the code error (E/L
# asymmetry, opposite ways round) and the carrier phase error (prompt argument,
# opposite signs). Same correlator layout and same sample count, so the two carry
# equal weight in all three loops and the mean is the midpoint.
const DISAGREEING_DRIVER_RECORD = EarlyPromptLateCorrelator(
    [complex(0.90, 0.00), complex(1.00, 0.10), complex(0.60, 0.00)],
    0.5,
)
const DISAGREEING_PASSENGER_RECORD = EarlyPromptLateCorrelator(
    [complex(0.55, 0.00), complex(1.00, -0.25), complex(0.95, 0.00)],
    0.5,
)
const RECORD_SAMPLES = 5000

# A two-signal L1 C/A satellite whose records arrive from an external producer.
# `NoCN0Estimator` keeps the noise-density plumbing out of it — this is a test
# about the loop update, and the C/N₀ path has its own.
function externally_fed_state(; combining = true, passenger_delay = 0.0s)
    estimator = ConventionalAssistedPLLAndDLL()
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
    set_group_delay!(ts, 1, 1, 2, passenger_delay)
    ts
end

# One estimate call over one driver record and, optionally, one passenger record
# that completed at the same sample. Returns the satellite's Dopplers after it.
function one_update(driver_record; passenger_record = nothing, kwargs...)
    ts = externally_fed_state(; kwargs...)
    n = RECORD_SAMPLES
    isnothing(passenger_record) ||
        append_correlator_output!(ts, CorrelatorOutput(passenger_record, n, n), 1, 1, 2)
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

    # The two records really do disagree — by ~18 Hz of carrier Doppler and
    # ~0.7 Hz of code Doppler — so the midpoint below is nowhere near either end
    # and no tolerance can hide a combiner that dropped one of them.
    @test abs(driver_alone[1] - passenger_alone[1]) > 1.0Hz
    @test abs(driver_alone[2] - passenger_alone[2]) > 0.1Hz

    # Equal weights, so the combined update is exactly the midpoint. This is the
    # assertion that fails if the fold pushes `own_pll` / `own_dll` into the
    # filters instead of `_normalized_mean`'s output, or if the driver's own
    # contribution is dropped from the mean.
    @test combined[1] ≈ (driver_alone[1] + passenger_alone[1]) / 2 rtol = 1e-12
    @test combined[2] ≈ (driver_alone[2] + passenger_alone[2]) / 2 rtol = 1e-12

    # …and with combining off the passenger's record must not reach the loop at
    # all, bit-identically to the run where it was never appended. It is still
    # consumed as a record — the estimator owns clearing those buffers.
    driver_only_run = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        combining = false,
    )
    @test driver_only_run === driver_alone
end

@testset "a supplied group delay moves the tracked code phase, by the right amount" begin
    gpsl1 = GPSL1CA()
    driver_alone = one_update(DISAGREEING_DRIVER_RECORD)
    passenger_alone = one_update(DISAGREEING_PASSENGER_RECORD)
    combined = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )

    # The initial code Doppler every record in these runs was evaluated at (the
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
    # The direction, stated without reference to the arithmetic above.
    @test _code_only(delayed) > _code_only(combined)
    # The carrier loops never see the delay.
    @test delayed[1] === combined[1]

    # Linear in the delay, i.e. the conversion really is `δ · f_code` chips and
    # not a saturating or once-only correction.
    twice = one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
        passenger_delay = 2δ,
    )
    @test _code_only(twice) - _code_only(combined) ≈
          2 * (_code_only(delayed) - _code_only(combined)) rtol = 1e-12
end

@testset "a pending accumulator is consumed by the next call's driver record" begin
    # The half of the cross-chunk accumulator that the `track`-based testset
    # above cannot reach: it stops before the long-period driver ever completes a
    # record, so nothing there shows the parked sums *arriving* at a loop update.
    # Here the two records are split across two estimate calls, and the result
    # must equal the single call that saw both — which is exactly the claim that
    # `pending_combining_sums` is seeded into the fold rather than started at
    # zero, and that a driver record consumes it.
    n = RECORD_SAMPLES
    ts = externally_fed_state()
    append_correlator_output!(
        ts,
        CorrelatorOutput(DISAGREEING_PASSENGER_RECORD, n, n),
        1,
        1,
        2,
    )
    estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))

    # No driver record, so the Dopplers hold and the passenger's discriminators
    # are parked rather than discarded.
    @test get_carrier_doppler(ts, 1) === 1000.0Hz
    parked =
        Tracking.get_doppler_estimator_state(get_sat_state(ts, 1)).pending_combining_sums
    @test parked.pll_weight > 0
    @test parked.dll_weight > 0

    append_correlator_output!(
        ts,
        CorrelatorOutput(DISAGREEING_DRIVER_RECORD, n, n),
        1,
        1,
        1,
    )
    estimate_dopplers_and_filter_prompt!(ts, (L1 = FS,))
    split_over_two_calls = (get_carrier_doppler(ts, 1), get_code_doppler(ts, 1))

    @test split_over_two_calls === one_update(
        DISAGREEING_DRIVER_RECORD;
        passenger_record = DISAGREEING_PASSENGER_RECORD,
    )
    # …and the driver record consumed the accumulator on its way through.
    @test Tracking.get_doppler_estimator_state(get_sat_state(ts, 1)).pending_combining_sums ==
          DiscriminatorAccumulator()
end

@testset "combining is a group property, and a mis-ordered group is accepted" begin
    # Combining is declared on the `SignalGroup`, next to its band and antenna
    # count, because its one precondition — `signals[1]` integrates longest — is a
    # property of the group's signal tuple, which every satellite of the group
    # shares. Nothing per-satellite carries it.
    combining_state = externally_fed_state(; combining = true)
    driver_only_state = externally_fed_state(; combining = false)
    @test only(Tuple(combining_state.groups)).discriminator_combining
    @test !only(Tuple(driver_only_state.groups)).discriminator_combining
    @test !hasfield(Tracking.SatConventionalPLLAndDLL, :discriminator_combining)
    @test !hasfield(Tracking.ConventionalPLLAndDLL, :discriminator_combining)

    # The ordering is a documented precondition and not an enforced one. A 1 ms
    # C/A driver with a 10 ms L1C-P passenger is the ordering that throws away
    # most of the gain — but it costs gain, never correctness, so it is built
    # rather than refused.
    @test TrackState(;
        signals = (gps_l1 = (GPSL1CA(), GPSL1C_P()),),
        discriminator_combining = true,
    ) isa TrackState
    @test Tracking.SignalGroup((GPSL1CA(), GPSL1C_P()); discriminator_combining = true) isa
          Tracking.SignalGroup
    # As are the right way round, either ordering with combining off, and equal
    # integration lengths in any order.
    @test TrackState(;
        signals = (gps_l1 = (GPSL1C_P(), GPSL1CA()),),
        discriminator_combining = true,
    ) isa TrackState
    @test TrackState(; signals = (gps_l1 = (GPSL1CA(), GPSL1C_P()),)) isa TrackState
    @test TrackState(;
        signals = (gps_l1 = (GPSL1C_D(), GPSL1C_P()),),
        discriminator_combining = true,
    ) isa TrackState
end

@testset "lengthening a passenger past the driver is accepted too" begin
    # The other way the ordering can break is a later
    # `set_preferred_num_code_blocks_to_integrate!`, and it is not refused either:
    # the setter is about this signal's integration length and says nothing about
    # what combining will make of it.
    prn, carrier_doppler = 1, 1000.0Hz
    ts = identical_l1ca_track_state(2, prn, carrier_doppler)
    set_preferred_num_code_blocks_to_integrate!(ts, 1, prn, 2, 2)
    @test get_preferred_num_code_blocks_to_integrate(get_sat_state(ts, prn), 2) == 2
    # Raising the driver is the direction that keeps the gain, and is no more or
    # less accepted.
    set_preferred_num_code_blocks_to_integrate!(ts, 1, prn, 1, 2)
    @test get_preferred_num_code_blocks_to_integrate(get_sat_state(ts, prn), 1) == 2

    # And with combining off the ordering is nobody's business at all.
    off =
        identical_l1ca_track_state(2, prn, carrier_doppler; discriminator_combining = false)
    set_preferred_num_code_blocks_to_integrate!(off, 1, prn, 2, 2)
    @test get_preferred_num_code_blocks_to_integrate(get_sat_state(off, prn), 2) == 2
end

end
