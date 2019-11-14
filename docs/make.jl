using Documenter, Tracking, GNSSSignals

makedocs(
    sitename="Tracking.jl",
    format = Documenter.HTML(prettyurls = false),
    modules = [Tracking, GNSSSignals],
    pages = [
        "index.md",
        "track.md",
        "tracking_state.md",
        "tracking_results.md",
        "loop_filter.md",
        "correlator.md",
        "cn0_estimator.md"
        ]
)

deploydocs(
    repo = "github.com/JuliaGNSS/Tracking.jl.git",
)
