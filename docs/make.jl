using Documenter, DocumenterInterLinks, Tracking, GNSSSignals, TrackingLoopFilters, TrackingLoops

const DOCTEST_SETUP = :(using Tracking, TrackingLoops, GNSSSignals, TrackingLoopFilters; using Tracking: Hz; using Unitful: u_str)
DocMeta.setdocmeta!(Tracking, :DocTestSetup, DOCTEST_SETUP; recursive=true)

# The per-record loop arithmetic — correlators, discriminators, the bit buffer,
# the C/N0 estimators, the Doppler estimators — is TrackingLoops' API, which
# users load next to Tracking and which Tracking does not re-export. Its own
# manual documents it; these pages link there with `@extref` rather than
# repeating its docstrings.
const LINKS = InterLinks(
    "TrackingLoops" => "https://juliagnss.github.io/TrackingLoops.jl/stable/",
)

makedocs(
    sitename="Tracking.jl",
    format = Documenter.HTML(prettyurls = false),
    modules = [Tracking],
    plugins = [LINKS],
    doctest = true,
    checkdocs = :exports,  # only complain about undocumented *exported* symbols
    pages = [
        "index.md",
        "signals.md",
        "track.md",
        "tracking_state.md",
        "bit_sync.md",
        "loop_filter.md",
        "custom_doppler_estimator.md",
        "correlator.md",
        "cn0_estimator.md",
        "noise_estimator.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaGNSS/Tracking.jl.git",
    push_preview = true,
)
