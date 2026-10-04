# Precompile workload (PrecompileTools). The first `track!` of a session costs
# seconds of compilation per backend (several times that on embedded ARM), and a
# live receiver pays it after the stream has started (GNSSReceiver.jl#107). The
# kernels are specialised on the signal type, so every trackable signal is
# tracked here for a few code periods with both sample element types.
using PrecompileTools: @setup_workload, @compile_workload

# Every signal Tracking has a default correlator for. (The Galileo E5/E6 and
# BeiDou B1C/B2a signals are defined by GNSSSignals but not tracked yet.)
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

# Sampling frequency for one signal's workload: 24 samples per chip at the
# 1.023 MHz chip rate (BOC signals need their sub-chip rate, 12× for
# TMBOC/CBOC/QMBOC), four for the 5.115 and 10.23 MHz BPSK codes. Float64-valued
# like a receiver's `5e6Hz`, since the kernels specialise on the element type.
function _precompile_sampling_frequency(system)
    code_frequency = get_code_frequency(system)
    code_frequency <= 1.1e6Hz ? 24.0 * code_frequency : 4.0 * code_frequency
end

# One code period of a clean replica with a small carrier Doppler, so the loops
# have something to pull on.
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

# One signal's workload, as a function so every call inside it is statically
# dispatched and lands in the package image; a dynamically dispatched loop over
# the signals left methods uncached.
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
    @compile_workload begin
        # `map` over tuples unrolls, so each call below is a static dispatch.
        map(
            _PRECOMPILE_SIGNALS,
            signals,
        ) do system, (signal_f32, signal_i16, sampling_frequency)
            _precompile_track(system, signal_f32, signal_i16, sampling_frequency, backend16)
        end
    end
end
