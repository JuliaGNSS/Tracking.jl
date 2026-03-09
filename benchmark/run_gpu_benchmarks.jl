# Run GPU vs CPU benchmarks and print results as a markdown table.
# Usage: julia --project benchmark/run_gpu_benchmarks.jl

using Pkg
Pkg.instantiate()

include("bench_gpu_vs_cpu.jl")

suite = gpu_suite()
results = run(suite; verbose=true, seconds=5)

configs = sort(collect(results), by=first)

# Find all GPU backends present (exclude CPU and CPU-Threaded)
gpu_backends = String[]
has_threaded = false
for (_, backends) in configs
    for (name, _) in backends
        if name == "CPU-Threaded"
            has_threaded = true
        elseif name != "CPU" && name ∉ gpu_backends
            push!(gpu_backends, name)
        end
    end
end
sort!(gpu_backends)

function to_us(t_ns)
    round(t_ns / 1e3, digits=1)
end

println("\n## Benchmark Results (GPU vs CPU)\n")
println("Julia threads: $(Threads.nthreads())\n")

for gpu in gpu_backends
    if has_threaded
        println("| Config | $gpu (μs) | CPU-Threaded (μs) | CPU (μs) | vs CPU | vs Threaded |")
        println("|--------|--------:|--------:|--------:|--------:|--------:|")
    else
        println("| Config | $gpu (μs) | CPU (μs) | vs CPU |")
        println("|--------|--------:|--------:|--------:|")
    end
    for (config, backends) in configs
        haskey(backends, "CPU") || continue
        haskey(backends, gpu) || continue
        cpu_us = to_us(median(backends["CPU"]).time)
        gpu_us = to_us(median(backends[gpu]).time)
        vs_cpu = round(cpu_us / gpu_us, digits=2)
        if has_threaded
            if haskey(backends, "CPU-Threaded")
                thr_us = to_us(median(backends["CPU-Threaded"]).time)
                vs_thr = round(thr_us / gpu_us, digits=2)
                println("| $config | $gpu_us | $thr_us | $cpu_us | $(vs_cpu)x | $(vs_thr)x |")
            else
                println("| $config | $gpu_us | — | $cpu_us | $(vs_cpu)x | — |")
            end
        else
            println("| $config | $gpu_us | $cpu_us | $(vs_cpu)x |")
        end
    end
    println()
end
