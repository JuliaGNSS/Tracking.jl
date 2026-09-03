# Precompile workload (PrecompileTools). The first `track!` of a session costs
# 2.6 s of compilation for the float backend and another 1.7 s for the Int16
# backend on a workstation, three to four times that on an embedded ARM host,
# and a live receiver pays it on its first tracked satellite, i.e. after the
# stream has started (GNSSReceiver.jl#107). Run the shapes a receiver uses —
# both sample element types, the threaded backends, one satellite tracked
# through a few code periods so the loop filters, prompt filter, C/N₀ estimator
# and bit buffer all execute — so `Pkg.precompile` pays for them instead.
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    system = GPSL1CA()
    sampling_frequency = 4e6Hz
    carrier_doppler = 200.0Hz
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(system) + get_code_frequency(system)
    samples = 0:3999
    signal_f32 = ComplexF32.(
        cis.(2π .* 200.0 .* samples ./ 4e6) .*
        gen_code(4000, system, 1, sampling_frequency, code_frequency, 100.0),
    )
    signal_i16 = Complex{Int16}.(round.(signal_f32 .* 512))
    @compile_workload begin
        state = TrackState(system, [TrackedSat(system, 1, 100.0, 180.0Hz)])
        for _ = 1:3
            track!(signal_f32, state, sampling_frequency)
        end
        state16 = TrackState(system, [TrackedSat(system, 1, 100.0, 180.0Hz)])
        backend = Int16ThreadedDownconvertAndCorrelator(2^12)
        for _ = 1:3
            track!(signal_i16, state16, sampling_frequency; downconvert_and_correlator = backend)
        end
        estimate_cn0(state, 1)
        get_soft_bits(state, 1)
        get_carrier_doppler(state, 1)
        get_code_phase(state, 1)
    end
end
