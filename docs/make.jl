using Documenter, Tracking, GNSSSignals, TrackingLoopFilters, TrackingLoops
using Documenter: Remotes

DocMeta.setdocmeta!(Tracking, :DocTestSetup, :(using Tracking, GNSSSignals, TrackingLoopFilters; using Tracking: Hz, NumAnts; using Unitful: u_str); recursive=true)

# TrackingLoops is resolved from its repository rather than developed, so its
# checkout carries no `.git` and Documenter cannot work out where to point the
# "source" link of a docstring that lives there. Name the remote explicitly.
const TRACKINGLOOPS_ROOT = dirname(dirname(pathof(TrackingLoops)))

makedocs(
    sitename="Tracking.jl",
    remotes = Dict(
        TRACKINGLOOPS_ROOT =>
            (Remotes.GitHub("JuliaGNSS", "TrackingLoops.jl"), "main"),
    ),
    format = Documenter.HTML(prettyurls = false),
    # The per-record loop arithmetic — correlators, discriminators, the bit
    # buffer, the C/N0 estimators, the Doppler estimators — lives in
    # TrackingLoops and is re-exported here, so the docstrings this manual is
    # built from are its. A page naming `EarlyPromptLateCorrelator` means the
    # same binding either way.
    modules = [Tracking, TrackingLoops],
    doctest = true,
    checkdocs = :exports,  # only complain about undocumented *exported* symbols
    # TrackingLoops has to be in `modules` for its docstrings to be spliced into
    # the pages below at all — Documenter resolves a `@docs` entry only against
    # the listed modules. That also puts its exports under `checkdocs`, and this
    # manual is not the place that documents them: they are a package of their
    # own, re-exported here for compatibility. `checkdocs_ignored_modules` does
    # not help, because it only stops recursion into *sub*modules of what is
    # listed. So the missing-docs check reports and does not fail.
    warnonly = [:missing_docs],
    pages = [
        "index.md",
        "signals.md",
        "track.md",
        "tracking_state.md",
        "bit_sync.md",
        "loop_filter.md",
        "custom_doppler_estimator.md",
        "vector_tracking.md",
        "correlator.md",
        "cn0_estimator.md",
        "noise_estimator.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaGNSS/Tracking.jl.git",
    push_preview = true,
)
