using Documenter, Tracking, GNSSSignals, TrackingLoopFilters

DocMeta.setdocmeta!(Tracking, :DocTestSetup, :(using Tracking, GNSSSignals, TrackingLoopFilters; using Tracking: Hz, NumAnts; using Unitful: u_str); recursive=true)

makedocs(
    sitename="Tracking.jl",
    format = Documenter.HTML(prettyurls = false),
    modules = [Tracking],
    doctest = true,
    warnonly = [:missing_docs],
    pages = [
        "index.md",
        "track.md",
        "tracking_state.md",
        "loop_filter.md",
        "correlator.md",
        "cn0_estimator.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaGNSS/Tracking.jl.git",
)
