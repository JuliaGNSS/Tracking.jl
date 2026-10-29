# Precompile workload. The first `track!` of a session costs 2.6 s of compilation
# for the float backend and 1.7 s more for Int16 on a workstation, three to four
# times that on an embedded ARM host, and a live receiver pays it after its
# stream has started (GNSSReceiver.jl#107). The kernels specialise on the signal
# type, so each signal below tracks one satellite for a few code periods, with
# both sample types and the threaded backends a receiver uses.
using PrecompileTools: @setup_workload, @compile_workload

# The precompiled signals. Galileo E5/E6 and BeiDou B1C/B2a are tracked too
# (docs/src/signals.md) but not precompiled.
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
    # BeiDou B2b is left out: a clean replica yields an all-zero prompt, a NaN
    # Doppler and then an `InexactError` (#234).
    BeiDouB3I(),
)

# The BOC signals need at least their sub-chip rate (12× the chip rate for
# TMBOC/CBOC/QMBOC), so 1.023 MHz codes get 24 samples per chip and the faster
# BPSK codes four. Float64, like the `5e6Hz` a receiver passes: the kernels
# specialise on the frequency's element type.
function _precompile_sampling_frequency(signal)
    code_frequency = get_chip_rate(signal)
    code_frequency <= 1.1e6Hz ? 24.0 * code_frequency : 4.0 * code_frequency
end

# One code period of a clean replica with a small carrier Doppler for the loops
# to pull on.
function _precompile_samples(signal)
    sampling_frequency = _precompile_sampling_frequency(signal)
    num_samples =
        round(Int, get_code_length(signal) * sampling_frequency / get_chip_rate(signal))
    carrier_doppler = 200.0Hz
    code_frequency =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(signal) +
        get_chip_rate(signal)
    indices = 0:(num_samples-1)
    samples_f32 = ComplexF32.(
        cis.(2π .* 200.0 .* indices ./ ustrip(Hz, sampling_frequency)) .*
        gen_code(num_samples, signal, 1, sampling_frequency, code_frequency, 0.0),
    )
    samples_f32, Complex{Int16}.(round.(samples_f32 .* 512)), sampling_frequency
end

# A function, so every call is statically dispatched on the signal type; methods
# compiled through dynamic dispatch were not all cached.
function _precompile_track(signal, samples_f32, samples_i16, sampling_frequency, backend16)
    state = TrackState(signal, [TrackedSat(signal, 1, 0.0, 180.0Hz)])
    for _ = 1:3
        track!(samples_f32, state, sampling_frequency)
    end
    state16 = TrackState(signal, [TrackedSat(signal, 1, 0.0, 180.0Hz)])
    for _ = 1:3
        track!(
            samples_i16,
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
    samples = map(_precompile_samples, _PRECOMPILE_SIGNALS)
    backend16 = Int16ThreadedDownconvertAndCorrelator(2^12)
    @compile_workload begin
        # `map` over tuples unrolls, so each call below is a static dispatch.
        map(
            _PRECOMPILE_SIGNALS,
            samples,
        ) do signal, (samples_f32, samples_i16, sampling_frequency)
            _precompile_track(
                signal,
                samples_f32,
                samples_i16,
                sampling_frequency,
                backend16,
            )
        end
    end
end
