# Run GPU vs CPU benchmarks and print results as a markdown table.
# Usage: julia --project benchmark/run_gpu_benchmarks.jl

using Pkg
Pkg.instantiate()

include("bench_gpu_vs_cpu.jl")

suite = gpu_suite()
results = run(suite; verbose=true, seconds=5)

configs = sort(collect(results), by=first)

# Find all GPU backends present
gpu_backends = String[]
for (_, backends) in configs
    for (name, _) in backends
        if name != "CPU" && name ∉ gpu_backends
            push!(gpu_backends, name)
        end
    end
end
sort!(gpu_backends)

function to_us(t_ns)
    round(t_ns / 1e3, digits=1)
end

println("\n## Benchmark Results (GPU vs CPU)\n")

for gpu in gpu_backends
    println("| Config | GPU (μs) | CPU (μs) | GPU/CPU |")
    println("|--------|--------:|--------:|--------:|")
    for (config, backends) in configs
        haskey(backends, "CPU") || continue
        haskey(backends, gpu) || continue
        cpu_us = to_us(median(backends["CPU"]).time)
        gpu_us = to_us(median(backends[gpu]).time)
        ratio = round(gpu_us / cpu_us, digits=2)
        println("| $config | $gpu_us | $cpu_us | $(ratio)x |")
    end
    println()
end
