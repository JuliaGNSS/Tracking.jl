# Run GPU vs CPU benchmarks and print results as a markdown table.
# Usage: julia --project benchmark/run_gpu_benchmarks.jl

using Pkg
Pkg.instantiate()

include("bench_gpu_vs_cpu.jl")

suite = gpu_suite()
results = run(suite; verbose=true, seconds=5)

# Collect config order from the suite (sorted alphabetically)
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

println("\n## Benchmark Results (GPU vs CPU)\n")

# Build header
let header = "| Config | CPU (median)", separator = "|--------|-------------"
    for gpu in gpu_backends
        header *= " | $gpu (median) | Speedup"
        separator *= "|--------------:|-------:"
    end
    println(header, " |")
    println(separator, "|")
end

for (config, backends) in configs
    haskey(backends, "CPU") || continue
    cpu_med = median(backends["CPU"]).time
    row = "| $config | $(BenchmarkTools.prettytime(cpu_med))"
    for gpu in gpu_backends
        if haskey(backends, gpu)
            gpu_med = median(backends[gpu]).time
            speedup = cpu_med / gpu_med
            row *= " | $(BenchmarkTools.prettytime(gpu_med)) | $(round(speedup, digits=1))x"
        else
            row *= " | N/A | N/A"
        end
    end
    println(row, " |")
end
