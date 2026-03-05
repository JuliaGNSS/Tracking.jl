# Benchmark: GPU backends (CUDA, AMDGPU) vs CPU
#
# Runs all available backends and prints a comparison table.
# Usage: julia --project benchmark/bench_gpu_vs_cpu.jl

using Tracking
using Tracking:
    KADownconvertAndCorrelator,
    CPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    downconvert_and_correlate
using GNSSSignals: GPSL1, GalileoE1B
using Unitful: Hz
using BenchmarkTools

# Detect available GPU backends
const HAS_CUDA = try
    @eval using CUDA: CUDA, cu, CuArray
    CUDA.functional()
catch
    false
end

const HAS_AMDGPU = try
    @eval using AMDGPU: AMDGPU, ROCArray
    true
catch
    false
end

# Access old CUDA extension if available
const HAS_CUDA_EXT = HAS_CUDA && try
    ext = Base.get_extension(Tracking, :TrackingCUDAExt)
    @eval const GPUDownconvertAndCorrelator = $ext.GPUDownconvertAndCorrelator
    true
catch
    false
end

function run_benchmarks()
    gpsl1 = GPSL1()
    gal = GalileoE1B()
    intermediate_frequency = 0.0Hz

    println("=" ^ 90)
    println("Tracking.jl — GPU vs CPU benchmark")
    println()
    println("Backends:")
    println("  CPU:    LoopVectorization")
    HAS_CUDA     && println("  CUDA:   ", CUDA.name(CUDA.device()))
    HAS_CUDA_EXT && println("  CUDA (old extension): texture memory")
    HAS_AMDGPU   && println("  AMDGPU: ", AMDGPU.device())
    println("=" ^ 90)
    println()

    # Build column headers dynamically
    cols = ["CPU"]
    HAS_CUDA_EXT && push!(cols, "CUDA-ext")
    HAS_CUDA     && push!(cols, "KA-CUDA-I64", "KA-CUDA-I32")
    HAS_AMDGPU   && push!(cols, "KA-AMD-I64", "KA-AMD-I32")

    header = rpad("Config", 24) * join(rpad.(cols .* " (μs)", 16))
    println(header)
    println("-" ^ length(header))

    configs = [
        ("L1 4sat/5K",     (gpsl1,), [4],  5e6Hz,   5000,  32, 1000.0),
        ("L1 16sat/5K",    (gpsl1,), [16], 5e6Hz,   5000,  32, 1000.0),
        ("L1 16sat/25K",   (gpsl1,), [16], 25e6Hz, 25000,  32, 1000.0),
        ("E1B 4sat/25K",   (gal,),   [4],  25e6Hz, 25000,  50, 100.0),
        ("E1B 16sat/25K",  (gal,),   [16], 25e6Hz, 25000,  50, 100.0),
        ("E1B 4sat/100K",  (gal,),   [4],  25e6Hz, 100000, 50, 100.0),
        ("E1B 16sat/100K", (gal,),   [16], 25e6Hz, 100000, 50, 100.0),
    ]

    multi_configs = [
        ("4L1+4E1B/25K",   (gpsl1, gal), [4,  4],  25e6Hz, 25000),
        ("8L1+8E1B/25K",   (gpsl1, gal), [8,  8],  25e6Hz, 25000),
        ("8L1+8E1B/100K",  (gpsl1, gal), [8,  8],  25e6Hz, 100000),
        ("16L1+16E1B/25K", (gpsl1, gal), [16, 16], 25e6Hz, 25000),
    ]

    for (label, systems, nsats_list, sfreq, nsamp, prn_max, code_dop) in configs
        results = bench_config(systems, nsats_list, sfreq, nsamp, prn_max, code_dop)
        print_row(label, results, cols)
    end

    println()
    println("--- Multi-system (GPSL1 + GalileoE1B) ---")
    println(header)
    println("-" ^ length(header))

    for (label, systems, nsats_list, sfreq, nsamp) in multi_configs
        results = bench_config(systems, nsats_list, sfreq, nsamp, 50, 100.0)
        print_row(label, results, cols)
    end
end

function bench_config(systems, nsats_list, sfreq, nsamp, prn_max, code_dop)
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

    results = Dict{String,Float64}()

    # CPU
    cpu_dc = CPUDownconvertAndCorrelator(Val(sfreq))
    downconvert_and_correlate(cpu_dc, signal_cpu, ts, 1, sfreq, 0.0Hz)
    b = @benchmark downconvert_and_correlate($cpu_dc, $signal_cpu, $ts, 1, $sfreq, $(0.0Hz)) samples=100
    results["CPU"] = round(median(b).time / 1000, digits=1)

    # Old CUDA extension
    if HAS_CUDA_EXT
        signal_cu = cu(signal_cpu)
        cuda_dc = GPUDownconvertAndCorrelator(Tuple(all_sss), nsamp)
        downconvert_and_correlate(cuda_dc, signal_cu, ts, 1, sfreq, 0.0Hz)
        b = @benchmark downconvert_and_correlate($cuda_dc, $signal_cu, $ts, 1, $sfreq, $(0.0Hz)) samples=100
        results["CUDA-ext"] = round(median(b).time / 1000, digits=1)
    end

    # KA + CUDA
    if HAS_CUDA
        signal_cu = cu(signal_cpu)
        for (plabel, ptype) in [("KA-CUDA-I64", Int64), ("KA-CUDA-I32", Int32)]
            ka = KADownconvertAndCorrelator(systems, CuArray; phase_type=ptype, max_sats=max(total_sats, 4))
            downconvert_and_correlate(ka, signal_cu, ts, 1, sfreq, 0.0Hz)
            b = @benchmark downconvert_and_correlate($ka, $signal_cu, $ts, 1, $sfreq, $(0.0Hz)) samples=100
            results[plabel] = round(median(b).time / 1000, digits=1)
        end
    end

    # KA + AMDGPU
    if HAS_AMDGPU
        signal_roc = ROCArray(signal_cpu)
        for (plabel, ptype) in [("KA-AMD-I64", Int64), ("KA-AMD-I32", Int32)]
            ka = KADownconvertAndCorrelator(systems, ROCArray; phase_type=ptype, max_sats=max(total_sats, 4))
            downconvert_and_correlate(ka, signal_roc, ts, 1, sfreq, 0.0Hz)
            b = @benchmark downconvert_and_correlate($ka, $signal_roc, $ts, 1, $sfreq, $(0.0Hz)) samples=100
            results[plabel] = round(median(b).time / 1000, digits=1)
        end
    end

    results
end

function print_row(label, results, cols)
    row = rpad(label, 24)
    for col in cols
        val = get(results, col, nothing)
        row *= rpad(val === nothing ? "-" : "$val", 16)
    end
    println(row)
end

run_benchmarks()
