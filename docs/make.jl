using Documenter, Tracking

makedocs(
    sitename="Tracking.jl",
    format = Documenter.HTML(prettyurls = false),
    modules = [Tracking, GNSSSignals],
    pages = [
        "index.md",
        "track.md",
        "tracking_state.md",
        "tracking_results.md",
        "correlator.md",
        "cn0_estimator.md"
        ]
)
