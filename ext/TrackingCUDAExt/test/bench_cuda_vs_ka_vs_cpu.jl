# Benchmark: old CUDA extension vs new KernelAbstractions.jl (CUDA backend) vs CPU
#
# Runs on Buildkite CI with CUDA GPUs to compare all three implementations.

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
using CUDA: CUDA, cu, CuArray, functional
using KernelAbstractions
using BenchmarkTools

# Access extension types
const TrackingCUDAExt = Base.get_extension(Tracking, :TrackingCUDAExt)
const GPUDownconvertAndCorrelator = TrackingCUDAExt.GPUDownconvertAndCorrelator

function run_benchmarks()
    if !functional()
        println("CUDA not functional, skipping benchmarks")
        return
    end

    gpsl1 = GPSL1()
    gal = GalileoE1B()
    intermediate_frequency = 0.0Hz

    println("=" ^ 80)
    println("CUDA vs KernelAbstractions.jl (CUDA backend) vs CPU")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("=" ^ 80)
    println()

    header = rpad("Config", 24) *
        rpad("CUDA (μs)", 14) *
        rpad("KA (μs)", 14) *
        rpad("CPU (μs)", 14) *
        "KA/CUDA"
    println(header)
    println("-" ^ length(header))

    for (label, system, nsats, sfreq, nsamp, prn_max, code_dop) in [
        ("L1 4sat/5K",    gpsl1, 4,  5e6Hz,   5000,  32, 1000.0),
        ("L1 16sat/5K",   gpsl1, 16, 5e6Hz,   5000,  32, 1000.0),
        ("L1 16sat/25K",  gpsl1, 16, 25e6Hz, 25000,  32, 1000.0),
        ("E1B 4sat/25K",  gal,   4,  25e6Hz, 25000,  50, 100.0),
        ("E1B 16sat/25K", gal,   16, 25e6Hz, 25000,  50, 100.0),
        ("E1B 4sat/100K", gal,   4,  25e6Hz, 100000, 50, 100.0),
        ("E1B 16sat/100K",gal,   16, 25e6Hz, 100000, 50, 100.0),
    ]
        sat_states = [SatState(system, mod1(i, prn_max), 10.5 + i * 0.1, (code_dop + i * 10) * Hz) for i in 1:nsats]
        sss = SystemSatsState(system, sat_states)
        ts = TrackState((sss,))

        signal_cpu = rand(ComplexF32, nsamp)
        signal_cu = cu(signal_cpu)

        # --- Old CUDA extension ---
        cuda_dc = GPUDownconvertAndCorrelator((sss,), nsamp)
        downconvert_and_correlate(cuda_dc, signal_cu, ts, 1, sfreq, intermediate_frequency)
        b_cuda = @benchmark downconvert_and_correlate($cuda_dc, $signal_cu, $ts, 1, $sfreq, $intermediate_frequency) samples=100
        t_cuda = round(median(b_cuda).time / 1000, digits=1)

        # --- KA with CUDA backend ---
        ka = KADownconvertAndCorrelator((system,), CuArray; max_sats=max(nsats, 4))
        downconvert_and_correlate(ka, signal_cu, ts, 1, sfreq, intermediate_frequency)
        b_ka = @benchmark downconvert_and_correlate($ka, $signal_cu, $ts, 1, $sfreq, $intermediate_frequency) samples=100
        t_ka = round(median(b_ka).time / 1000, digits=1)

        # --- CPU ---
        cpu_dc = CPUDownconvertAndCorrelator(Val(sfreq))
        downconvert_and_correlate(cpu_dc, signal_cpu, ts, 1, sfreq, intermediate_frequency)
        b_cpu = @benchmark downconvert_and_correlate($cpu_dc, $signal_cpu, $ts, 1, $sfreq, $intermediate_frequency) samples=100
        t_cpu = round(median(b_cpu).time / 1000, digits=1)

        ratio = round(t_ka / t_cuda, digits=2)
        println(
            rpad(label, 24),
            rpad("$t_cuda", 14),
            rpad("$t_ka", 14),
            rpad("$t_cpu", 14),
            "$(ratio)x",
        )
    end

    # Multi-system benchmarks
    println()
    println("--- Multi-system (GPSL1 + GalileoE1B) ---")
    println(header)
    println("-" ^ length(header))

    for (nsats_l1, nsats_gal, sfreq, nsamp) in [
        (4,  4,  25e6Hz, 25000),
        (8,  8,  25e6Hz, 25000),
        (8,  8,  25e6Hz, 100000),
        (16, 16, 25e6Hz, 25000),
    ]
        label = "$(nsats_l1)L1+$(nsats_gal)E1B/$(Int(nsamp÷1000))K"
        sat_l1 = [SatState(gpsl1, mod1(i, 32), 10.5 + i * 0.1, (1000.0 + i * 10) * Hz) for i in 1:nsats_l1]
        sat_gal = [SatState(gal, mod1(i, 50), 11.0 + i * 0.1, (100.0 + i * 10) * Hz) for i in 1:nsats_gal]
        sss_l1 = SystemSatsState(gpsl1, sat_l1)
        sss_gal = SystemSatsState(gal, sat_gal)
        ts = TrackState((sss_l1, sss_gal))
        total_sats = nsats_l1 + nsats_gal

        signal_cpu = rand(ComplexF32, nsamp)
        signal_cu = cu(signal_cpu)

        # --- Old CUDA extension ---
        # Multi-system with different GNSS types is unsupported by the old extension
        # (GPUDownconvertAndCorrelator requires homogeneous NTuple type)
        t_cuda = "N/A"

        # --- KA with CUDA backend ---
        ka = KADownconvertAndCorrelator((gpsl1, gal), CuArray; max_sats=total_sats)
        downconvert_and_correlate(ka, signal_cu, ts, 1, sfreq, intermediate_frequency)
        b_ka = @benchmark downconvert_and_correlate($ka, $signal_cu, $ts, 1, $sfreq, $intermediate_frequency) samples=100
        t_ka = round(median(b_ka).time / 1000, digits=1)

        # --- CPU ---
        cpu_dc = CPUDownconvertAndCorrelator(Val(sfreq))
        downconvert_and_correlate(cpu_dc, signal_cpu, ts, 1, sfreq, intermediate_frequency)
        b_cpu = @benchmark downconvert_and_correlate($cpu_dc, $signal_cpu, $ts, 1, $sfreq, $intermediate_frequency) samples=100
        t_cpu = round(median(b_cpu).time / 1000, digits=1)

        println(
            rpad(label, 24),
            rpad("$t_cuda", 14),
            rpad("$t_ka", 14),
            rpad("$t_cpu", 14),
            "N/A",
        )
    end
end

run_benchmarks()
