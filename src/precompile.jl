# Precompile workload (PrecompileTools). The first `track!` of a session costs
# 2.6 s of compilation for the float backend and another 1.7 s for the Int16
# backend on a workstation, three to four times that on an embedded ARM host,
# and a live receiver pays it on its first tracked satellite, i.e. after the
# stream has started (GNSSReceiver.jl#107). The tracking kernels are
# specialised on the signal type, so every signal this package can track is
# tracked here — one satellite through a few code periods, with both sample
# element types and the threaded backends a receiver uses — so the loop
# filters, prompt filter, C/N₀ estimator and bit buffer of each compile at
# precompile time instead of on the first satellite of a run.
using PrecompileTools: @setup_workload, @compile_workload

# The signals whose kernels are compiled ahead of time. Not every signal
# Tracking can track: Galileo E5a/E5b/E6 and BeiDou B1C/B2a/B2b have default
# correlators too, and each one added costs precompile time and image size. They
# are left out until a receiver reports paying for them at run time; B2b has a
# second reason, below.
const _PRECOMPILE_SIGNALS = (
    GPSL1CA(),
    GPSL1C_D(),
    GPSL1C_P(),
    GPSL2CM(),
    GPSL2CL(),
    GPSL5I(),
    GPSL5Q(),
    GalileoE1B(),
    GalileoE1B_BOC11(),
    GalileoE1C(),
    GalileoE1C_BOC11(),
    BeiDouB1I(),
    # BeiDou B2b is left out for now: tracking a clean B2b replica yields an
    # all-zero prompt after the first code period, the discriminators turn that
    # into a NaN Doppler and the next `track!` throws `InexactError` — a
    # robustness bug in its own right, reported separately.
    BeiDouB3I(),
)

# Sampling frequency for one signal's workload. The BOC-modulated signals (GPS
# L1C, Galileo E1, BeiDou B1C) need at least their sub-chip rate — twelve times
# the chip rate for TMBOC/CBOC/QMBOC — so every signal at the 1.023 MHz chip
# rate is sampled at 24 chips per sample period; the 5.115 and 10.23 MHz BPSK
# codes get four samples per chip.
# Float64-valued, like the `5e6Hz` a receiver passes: the kernels specialise on
# the frequency's element type, and an `Int`-valued rate would compile a
# specialisation nobody calls.
function _precompile_sampling_frequency(system)
    code_frequency = get_code_frequency(system)
    code_frequency <= 1.1e6Hz ? 24.0 * code_frequency : 4.0 * code_frequency
end

# One code period of a clean replica at four samples per chip, with a small
# carrier Doppler, so the loops have something to pull on.
function _precompile_signal(system)
    sampling_frequency = _precompile_sampling_frequency(system)
    num_samples = round(
        Int,
        get_code_length(system) * sampling_frequency / get_code_frequency(system),
    )
    carrier_doppler = 200.0Hz
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(system) +
        get_code_frequency(system)
    samples = 0:(num_samples-1)
    signal_f32 = ComplexF32.(
        cis.(2π .* 200.0 .* samples ./ ustrip(Hz, sampling_frequency)) .*
        gen_code(num_samples, system, 1, sampling_frequency, code_frequency, 0.0),
    )
    signal_f32, Complex{Int16}.(round.(signal_f32 .* 512)), sampling_frequency
end

# One multi-signal satellite's workload under one estimator. Nothing in the
# single-signal workload below reaches the shared multi-signal fold, any of the
# accumulator types, or `VectorPLLAndDLL` at all — so without this a receiver
# compiles all of that on its first pilot/data satellite, which is exactly the
# "after the stream has started" cost the rest of this file exists to remove.
#
# Galileo E1B+E1C because it is the shortest path through the machinery that has
# both a data and a pilot component: one band, one sampling frequency, a real
# secondary code on the pilot, and a passenger whose records a driver record has
# to window. `vt_on` is left unset — a meaningful `vt_on` pass needs a navigation
# filter's corrections, and the `VectorPLLAndDLL` call below already compiles the
# state type and the traversal its accumulators seed.
#
# Called once per estimator rather than looped over them, for the same reason
# `_precompile_track` is: a loop over a heterogeneous tuple dispatches
# dynamically and the methods compiled that way are not all cached.
function _precompile_track_multi_signal(estimator, signal_f32, sampling_frequency)
    signals = (GalileoE1B(), GalileoE1C())
    sat = TrackedSat(signals, 1, 0.0, 180.0Hz; doppler_estimator = estimator)
    state = TrackState(first(signals), sat; doppler_estimator = estimator)
    for _ = 1:3
        track!(signal_f32, state, sampling_frequency)
    end
    # Via the satellite rather than the `TrackState`, so the signal selector is
    # unambiguously a selector and not a group index.
    sat = get_sat_state(state, 1)
    estimate_cn0(sat, GalileoE1C)
    get_soft_bits(sat, GalileoE1B)
    get_code_phase(sat)
    nothing
end

# A clean E1B+E1C composite: both components carry real power, so neither
# component's prompt collapses to the all-zero correlation that turns a
# discriminator into a NaN (see the BeiDou B2b exclusion above).
function _precompile_multi_signal(sampling_frequency)
    num_samples = round(
        Int,
        get_code_length(GalileoE1B()) * sampling_frequency /
        get_code_frequency(GalileoE1B()),
    )
    samples = 0:(num_samples-1)
    carrier = cis.(2π .* 200.0 .* samples ./ ustrip(Hz, sampling_frequency))
    component(system) = gen_code(
        num_samples,
        system,
        1,
        sampling_frequency,
        200.0Hz * get_code_center_frequency_ratio(system) + get_code_frequency(system),
        0.0,
    )
    ComplexF32.(carrier .* (component(GalileoE1B()) .+ component(GalileoE1C())))
end

# One signal's workload, as a function so every call inside it is statically
# dispatched on the concrete signal type — the specialisations then land in
# the package image. Iterating the signal tuple in a loop instead dispatches
# dynamically, and the methods compiled that way were not all cached (the
# first `track!` of a session still cost ~1 s).
function _precompile_track(system, signal_f32, signal_i16, sampling_frequency, backend16)
    state = TrackState(system, [TrackedSat(system, 1, 0.0, 180.0Hz)])
    for _ = 1:3
        track!(signal_f32, state, sampling_frequency)
    end
    state16 = TrackState(system, [TrackedSat(system, 1, 0.0, 180.0Hz)])
    for _ = 1:3
        track!(
            signal_i16,
            state16,
            sampling_frequency;
            downconvert_and_correlator = backend16,
        )
    end
    estimate_cn0(state, 1)
    get_soft_bits(state, 1)
    get_carrier_doppler(state, 1)
    get_code_phase(state, 1)
    nothing
end

@setup_workload begin
    signals = map(_precompile_signal, _PRECOMPILE_SIGNALS)
    backend16 = Int16ThreadedDownconvertAndCorrelator(2^12)
    e1_sampling_frequency = _precompile_sampling_frequency(GalileoE1B())
    e1_signal_f32 = _precompile_multi_signal(e1_sampling_frequency)
    @compile_workload begin
        # `map` over tuples unrolls, so each call below is a static dispatch.
        map(
            _PRECOMPILE_SIGNALS,
            signals,
        ) do system, (signal_f32, signal_i16, sampling_frequency)
            _precompile_track(system, signal_f32, signal_i16, sampling_frequency, backend16)
        end
        _precompile_track_multi_signal(
            ConventionalAssistedPLLAndDLL(; discriminator_combining = true),
            e1_signal_f32,
            e1_sampling_frequency,
        )
        _precompile_track_multi_signal(
            VectorPLLAndDLL(; discriminator_combining = true),
            e1_signal_f32,
            e1_sampling_frequency,
        )
    end
end
