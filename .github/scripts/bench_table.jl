# Build a PR-comment markdown table from benchpkg result JSONs, reporting the
# *minimum* time per benchmark instead of AirspeedVelocity's median. The minimum
# is robust to noisy-neighbour contention on shared CI runners (it captures the
# least-disturbed sample), so it doesn't manufacture phantom regressions when the
# head build happens to land on a busier runner window.
#
# The comment has two parts:
#
#   1. **Backends vs Float32** — a head-to-head of the alternative `track!`
#      backends (integer `Int16`, and the bit-wise `OneBit` / `TwoBit` /
#      `TwoBitC2`) against the Float32 default. These backends are NEW on the PR
#      (no base counterpart), so a plain base-vs-head diff can't show their
#      speedup. Instead we pair each scenario's `Float32`, `Int16`, `OneBit`,
#      `TwoBit` and `TwoBitC2` leaves — registered side by side under
#      `track! Int16 vs Float32/<scenario>/{Float32,Int16,OneBit,TwoBit,TwoBitC2}`
#      — into one row and report Float32 / backend (>1 ⇒ that backend is faster).
#      The bit-wise backends are BPSK-only, so they have no cell for CBOC
#      (Galileo E1B) scenarios. Headers are acronyms (see the in-comment legend).
#
#   2. **Regression table** — every other benchmark, base rev vs PR head, the
#      usual "did this PR get slower?" view.
#
# Usage:  julia bench_table.jl <input_dir> <pkg> <base_rev> <head_rev> [label]
# The optional <label> (e.g. the runner OS) is appended to the headings so a
# multi-platform matrix can post one distinct comment per platform.
using JSON3

input_dir, pkg, base_rev, head_rev = ARGS[1], ARGS[2], ARGS[3], ARGS[4]
label = length(ARGS) >= 5 ? ARGS[5] : ""

readrev(rev) = open(joinpath(input_dir, "results_$(pkg)@$(rev).json"), "r") do io
    JSON3.read(io, Dict{String,Any})
end

# Recursively collect leaf trials (nodes carrying a "times" array) keyed by "a/b/c".
function leaves!(acc, node, prefix = "")
    if haskey(node, "times")
        acc[prefix] = node
    elseif haskey(node, "data")
        for (k, v) in node["data"]
            name = isempty(prefix) ? String(k) : prefix * "/" * String(k)
            leaves!(acc, v, name)
        end
    end
    acc
end

mintime(node) = minimum(Float64.(node["times"]))   # ns

function fmt_time(ns)
    unit, div = ns < 1e3 ? ("ns", 1.0) :
                ns < 1e6 ? ("μs", 1e3) :
                ns < 1e9 ? ("ms", 1e6) : ("s", 1e9)
    string(round(ns / div; sigdigits = 3), " ", unit)
end

fmt_mem(node) = string(round(Int, node["allocs"]), " allocs: ", round(Int, node["memory"]), " B")

# Ratio cell with a fast/slow marker. `faster_when` picks which direction is good:
# for regressions base/head (>1 ⇒ PR faster); for backends-vs-Float32 Float32/backend
# (>1 ⇒ that backend faster). ✅ ≥ 5 % faster, ⚠️ ≥ 5 % slower.
function fmt_ratio(r)
    s = string(round(r; sigdigits = 3))
    r >= 1.05 ? "$s ✅" : r <= 0.95 ? "$s ⚠️" : s
end

base = leaves!(Dict{String,Any}(), readrev(base_rev))
head = leaves!(Dict{String,Any}(), readrev(head_rev))

# The Int16-vs-Float32 backend head-to-head group (see benchmark/benchmarks.jl).
const INT16_GROUP = "track! Int16 vs Float32"

shortrev(rev) = length(rev) >= 16 && !occursin('/', rev) ? rev[1:8] * "…" : rev
headlbl = shortrev(head_rev)
baselbl = shortrev(base_rev)

io = IOBuffer()
suffix = isempty(label) ? "" : " — $label"
println(io, "## Benchmark Results (minimum time)$suffix")
println(io)
println(io, "Reporting the **minimum** over all samples (robust to shared-runner ",
            "contention), not the median.")
println(io)

# ── Section 1: Int16 vs Float32 backend head-to-head (PR head only) ──────────
# Pair "<INT16_GROUP>/<scenario>/Float32" with ".../Int16" from the head rev.
pairs = Dict{String,Dict{String,Any}}()   # scenario => Dict("Float32"=>node, "Int16"=>node)
order = String[]                           # first-seen scenario order
for k in sort(collect(keys(head)))
    startswith(k, INT16_GROUP * "/") || continue
    rest = k[lastindex(INT16_GROUP)+2:end]           # "<scenario>/<backend>"
    slash = findlast('/', rest)
    slash === nothing && continue
    scenario = rest[1:prevind(rest, slash)]
    backend = rest[nextind(rest, slash):end]
    d = get!(pairs, scenario) do
        push!(order, scenario)
        Dict{String,Any}()
    end
    d[backend] = head[k]
end

if !isempty(order)
    println(io, "### Alternative backends vs Float32 (`track!`, PR head)")
    println(io)
    # Compact acronym headers (a wide 5-backend table otherwise overflows the PR
    # comment); this legend is the key.
    println(io, "**Legend** — backends: `F32` Float32 (default) · `I16` Int16 · ",
                "`1b` OneBit · `2b` TwoBit (2-bit meas, 1-bit carrier) · ",
                "`2c` TwoBit (2-bit meas, 2-bit carrier). Time columns are the minimum ",
                "`track!` time; **`×B` = F32 / B** (so **>1 ⇒ backend B is faster** than ",
                "Float32), ✅ ≥ 5 % faster, ⚠️ ≥ 5 % slower. `1b`/`2b`/`2c` are BPSK-only, ",
                "so their cells are blank for CBOC (Galileo E1B) scenarios.")
    println(io)
    println(io, "| Scenario | F32 | I16 | 1b | 2b | 2c | ×I16 | ×1b | ×2b | ×2c |")
    println(io, "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for scenario in order
        d = pairs[scenario]
        f = get(d, "Float32", nothing)
        i = get(d, "Int16", nothing)
        b = get(d, "OneBit", nothing)
        t = get(d, "TwoBit", nothing)
        c = get(d, "TwoBitC2", nothing)
        cell(n) = n === nothing ? "" : fmt_time(mintime(n))
        ratio(n) = (f !== nothing && n !== nothing) ? fmt_ratio(mintime(f) / mintime(n)) : ""
        println(io, "| $scenario | $(cell(f)) | $(cell(i)) | $(cell(b)) | $(cell(t)) | ",
                    "$(cell(c)) | $(ratio(i)) | $(ratio(b)) | $(ratio(t)) | $(ratio(c)) |")
    end
    println(io)
end

# ── Section 2: base-vs-head regression table for every other benchmark ───────
# Stable ordering over the UNION of both revs' benchmarks (so rows new on the PR
# still appear), excluding the Int16-vs-Float32 group (covered above). Sort
# alphabetically, then push time_to_load last.
isint16(n) = startswith(n, INT16_GROUP * "/")
names = sort(collect(union(keys(base), keys(head))))
filter!(n -> !isint16(n) && n != "time_to_load", names)
(haskey(base, "time_to_load") || haskey(head, "time_to_load")) && push!(names, "time_to_load")

println(io, "<details open><summary>Time benchmarks (base vs PR head)</summary>")
println(io)
println(io, "Ratio = $baselbl / $headlbl: **>1 means the PR is faster**. ✅ ≥ 5 % faster, ",
            "⚠️ ≥ 5 % slower. A blank cell means the benchmark exists on only one revision ",
            "(🆕 = new on the PR, 🗑 = removed).")
println(io)
println(io, "|  | $baselbl | $headlbl | $baselbl / $headlbl |")
println(io, "|:--|--:|--:|--:|")
for n in names
    b = get(base, n, nothing); h = get(head, n, nothing)
    bcell = b === nothing ? "" : fmt_time(mintime(b))
    hcell = h === nothing ? "" : fmt_time(mintime(h))
    if b !== nothing && h !== nothing
        rcell = fmt_ratio(mintime(b) / mintime(h))
    else
        rcell = h === nothing ? "🗑" : "🆕"
    end
    println(io, "| $n | $bcell | $hcell | $rcell |")
end
println(io)
println(io, "</details>")
println(io)

# --- memory table ---
println(io, "<details><summary>Memory benchmarks (base vs PR head)</summary>")
println(io)
println(io, "|  | $baselbl | $headlbl |")
println(io, "|:--|--:|--:|")
for n in names
    b = get(base, n, nothing); h = get(head, n, nothing)
    bcell = b === nothing ? "" : fmt_mem(b)
    hcell = h === nothing ? "" : fmt_mem(h)
    println(io, "| $n | $bcell | $hcell |")
end
println(io)
println(io, "</details>")

print(String(take!(io)))
