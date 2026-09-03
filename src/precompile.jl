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

# Sampling frequency for one signal's workload. The BOC-modulated signals (GPS
# L1C, Galileo E1, BeiDou B1C) need at least their sub-chip rate — twelve times
# the chip rate for TMBOC/CBOC/QMBOC — so every signal at the 1.023 MHz chip
# rate is sampled at 24 chips per sample period; the 5.115 and 10.23 MHz BPSK
# codes get four samples per chip.
function _precompile_sampling_frequency(system)
    code_frequency = get_code_frequency(system)
    code_frequency <= 1.1e6Hz ? 24code_frequency : 4code_frequency
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

@setup_workload begin
    signals = map(_precompile_signal, _PRECOMPILE_SIGNALS)
    backend16 = Int16ThreadedDownconvertAndCorrelator(2^12)
    @compile_workload begin
        for (system, (signal_f32, signal_i16, sampling_frequency)) in
            zip(_PRECOMPILE_SIGNALS, signals)
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
        end
    end
end
