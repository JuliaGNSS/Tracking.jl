using BenchmarkTools

include("bench_cpu.jl")
const SUITE = BenchmarkGroup()

for (k, v) in cpu_suite()
    SUITE[k] = v
end

# GPU vs CPU benchmarks are run separately via Buildkite (benchmark/run_gpu_benchmarks.jl),
# not through AirspeedVelocity, since they require a GPU runner.

# ── Multi-satellite benchmarks (threaded if available, CPU fallback) ──────

function _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max=32, code_dop=1000.0)
    all_sss = []
    total_sats = 0
    for (si, sys) in enumerate(systems)
        ns = nsats_list[min(si, length(nsats_list))]
        pm = sys isa GPSL1 ? 32 : prn_max
        cd = sys isa GPSL1 ? 1000.0 : code_dop
        sats = [SatState(sys, mod1(i, pm), 10.5 + i * 0.1, (cd + i * 10) * Hz) for i = 1:ns]
        push!(all_sss, SystemSatsState(sys, sats))
        total_sats += ns
    end
    TrackState(Tuple(all_sss)), rand(ComplexF32, nsamp), total_sats
end

function bench_multi_sat(; systems, nsats_list, sfreq, nsamp, prn_max=32, code_dop=1000.0)
    ts, signal, total_sats = _make_multi_sat_state(; systems, nsats_list, nsamp, prn_max, code_dop)
    if isdefined(Tracking, :CPUThreadedDownconvertAndCorrelator)
        dc = Tracking.CPUThreadedDownconvertAndCorrelator(
            systems, Val(sfreq); max_sats=max(total_sats, 4), max_num_samples=nsamp,
        )
    else
        dc = CPUDownconvertAndCorrelator(Val(sfreq))
    end
    @benchmarkable downconvert_and_correlate($dc, $signal, $ts, 1, $sfreq, $(0.0Hz))
end

let gpsl1 = GPSL1(), gal = GalileoE1B()
    SUITE["multi-sat"]["L1 8sat/5K"] =
        bench_multi_sat(; systems=(gpsl1,), nsats_list=[8], sfreq=5e6Hz, nsamp=5000)
    SUITE["multi-sat"]["E1B 4sat/25K"] =
        bench_multi_sat(; systems=(gal,), nsats_list=[4], sfreq=25e6Hz, nsamp=25000, prn_max=50, code_dop=100.0)
    SUITE["multi-sat"]["8L1+8E1B/25K"] =
        bench_multi_sat(; systems=(gpsl1, gal), nsats_list=[8, 8], sfreq=25e6Hz, nsamp=25000, prn_max=50, code_dop=100.0)
end
