# GPU (CUDA, AMDGPU) vs CPU downconvert-and-correlate benchmarks
#
# Returns a BenchmarkGroup with entries for each available backend.

using BenchmarkTools
using GNSSSignals: GPSL1, GalileoE1B
using Unitful: Hz
using Tracking
using Tracking:
    CPUDownconvertAndCorrelator,
    KADownconvertAndCorrelator,
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
end

function gpu_suite()
    suite = BenchmarkGroup()
    gpsl1 = GPSL1()
    gal = GalileoE1B()

    configs = [
        ("L1 1sat/5K",    (gpsl1,),       [1],    5e6Hz,  5000,  32, 1000.0),
        ("L1 8sat/5K",    (gpsl1,),       [8],    5e6Hz,  5000,  32, 1000.0),
        ("E1B 8sat/25K",  (gal,),         [8],    25e6Hz, 25000, 50, 100.0),
        ("E1B 8sat/100K", (gal,),         [8],    25e6Hz, 100000,50, 100.0),
        ("8L1+8E1B/25K",  (gpsl1, gal),   [8, 8], 25e6Hz, 25000, 50, 100.0),
    ]

    for (label, systems, nsats_list, sfreq, nsamp, prn_max, code_dop) in configs
        all_sss = []
        total_sats = 0
        for (si, sys) in enumerate(systems)
            ns = nsats_list[min(si, length(nsats_list))]
            pm = sys isa GPSL1 ? 32 : prn_max
            cd = sys isa GPSL1 ? 1000.0 : code_dop
            sats = [SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i in 1:ns]
            push!(all_sss, SystemSatsState(sys, sats))
            total_sats += ns
        end
        ts = TrackState(Tuple(all_sss))
        signal_cpu = rand(ComplexF32, nsamp)

        # CPU baseline
        cpu_dc = CPUDownconvertAndCorrelator(Val(sfreq))
        suite[label]["CPU"] = @benchmarkable downconvert_and_correlate(
            $cpu_dc, $signal_cpu, $ts, 1, $sfreq, $(0.0Hz),
        )

        # KA + CUDA
        if HAS_CUDA
            signal_cu = cu(signal_cpu)
            ka = KADownconvertAndCorrelator(systems, CuArray; max_sats=max(total_sats, 4))
            suite[label]["KA-CUDA"] = @benchmarkable downconvert_and_correlate(
                $ka, $signal_cu, $ts, 1, $sfreq, $(0.0Hz),
            )
        end

        # KA + AMDGPU
        if HAS_AMDGPU
            signal_roc = ROCArray(signal_cpu)
            ka = KADownconvertAndCorrelator(systems, ROCArray; max_sats=max(total_sats, 4))
            suite[label]["KA-AMD"] = @benchmarkable downconvert_and_correlate(
                $ka, $signal_roc, $ts, 1, $sfreq, $(0.0Hz),
            )
        end
    end

    suite
end
