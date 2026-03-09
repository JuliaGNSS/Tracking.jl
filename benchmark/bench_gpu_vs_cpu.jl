# GPU (CUDA, AMDGPU) vs CPU downconvert-and-correlate benchmarks
#
# Returns a BenchmarkGroup with entries for each available backend.

using BenchmarkTools
using GNSSSignals: GPSL1, GalileoE1B
using Unitful: Hz
using Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate

if !@isdefined(HAS_CUDA)
    const HAS_CUDA = try
        @eval using CUDA: CUDA, cu, CuArray
        CUDA.functional()
    catch
        false
    end

    const HAS_AMDGPU = try
        @eval using AMDGPU: AMDGPU, ROCArray
        AMDGPU.functional()
    catch
        false
    end

    const AMD_GPU_CUS = HAS_AMDGPU ? try
        AMDGPU.HIP.properties(AMDGPU.device()).multiProcessorCount
    catch
        0
    end : 0
end

function gpu_suite()
    suite = BenchmarkGroup()
    gpsl1 = GPSL1()
    gal = GalileoE1B()

    configs = [
        ("L1 4sat/5K", (gpsl1,), [4], 5e6Hz, 5000, 32, 1000.0),
        ("L1 8sat/5K", (gpsl1,), [8], 5e6Hz, 5000, 32, 1000.0),
        ("E1B 4sat/25K", (gal,), [4], 25e6Hz, 25000, 50, 100.0),
        ("E1B 8sat/50K", (gal,), [8], 25e6Hz, 50000, 50, 100.0),
        ("8L1+8E1B/25K", (gpsl1, gal), [8, 8], 25e6Hz, 25000, 50, 100.0),
        ("8L1+8E1B/50K", (gpsl1, gal), [8, 8], 25e6Hz, 50000, 50, 100.0),
    ]

    # Limit GPU configs on low-CU GPUs to avoid freezing the display
    max_gpu_configs = AMD_GPU_CUS <= 8 ? 2 : length(configs)
    if AMD_GPU_CUS > 0 && max_gpu_configs < length(configs)
        @info "AMD GPU has only $AMD_GPU_CUS CUs — limiting GPU benchmarks to first $max_gpu_configs configs"
    end

    for (cfg_idx, (label, systems, nsats_list, sfreq, nsamp, prn_max, code_dop)) in enumerate(configs)
        all_sss = []
        total_sats = 0
        for (si, sys) in enumerate(systems)
            ns = nsats_list[min(si, length(nsats_list))]
            pm = sys isa GPSL1 ? 32 : prn_max
            cd = sys isa GPSL1 ? 1000.0 : code_dop
            sats = [
                SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i = 1:ns
            ]
            push!(all_sss, SystemSatsState(sys, sats))
            total_sats += ns
        end
        ts = TrackState(Tuple(all_sss))
        signal_cpu = rand(ComplexF32, nsamp)

        # CPU baseline
        cpu_dc = CPUDownconvertAndCorrelator(Val(sfreq))
        suite[label]["CPU"] = @benchmarkable downconvert_and_correlate(
            $cpu_dc,
            $signal_cpu,
            $ts,
            1,
            $sfreq,
            $(0.0Hz),
        )

        # CPU threaded (only available on branches with the threaded correlator)
        if isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
            cpu_threaded_dc = Tracking.CPUThreadedDownconvertAndCorrelator(
                systems,
                Val(sfreq);
                max_sats = max(total_sats, 4),
                max_num_samples = nsamp,
            )
            suite[label]["CPU-Threaded"] = @benchmarkable downconvert_and_correlate(
                $cpu_threaded_dc,
                $signal_cpu,
                $ts,
                1,
                $sfreq,
                $(0.0Hz),
            )
        end

        # KA + CUDA
        if HAS_CUDA && isdefined(Tracking, :KADownconvertAndCorrelator)
            signal_cu = cu(signal_cpu)
            ka = Tracking.KADownconvertAndCorrelator(systems, CuArray; max_sats = max(total_sats, 4))
            suite[label]["KA-CUDA"] = @benchmarkable downconvert_and_correlate(
                $ka,
                $signal_cu,
                $ts,
                1,
                $sfreq,
                $(0.0Hz),
            )
        end

        # KA + AMDGPU (skip heavy configs on low-CU GPUs)
        if HAS_AMDGPU && isdefined(Tracking, :KADownconvertAndCorrelator) && cfg_idx <= max_gpu_configs
            signal_roc = ROCArray(signal_cpu)
            ka =
                Tracking.KADownconvertAndCorrelator(systems, ROCArray; max_sats = max(total_sats, 4))
            suite[label]["KA-AMD"] = @benchmarkable downconvert_and_correlate(
                $ka,
                $signal_roc,
                $ts,
                1,
                $sfreq,
                $(0.0Hz),
            )
        end
    end

    suite
end
