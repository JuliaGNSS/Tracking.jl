using Documenter, Tracking, GNSSSignals, TrackingLoopFilters, TrackingLoops
using Documenter: Remotes

const DOCTEST_SETUP = :(using Tracking, TrackingLoops, GNSSSignals, TrackingLoopFilters; using Tracking: Hz; using Unitful: u_str)
DocMeta.setdocmeta!(Tracking, :DocTestSetup, DOCTEST_SETUP; recursive=true)
DocMeta.setdocmeta!(TrackingLoops, :DocTestSetup, DOCTEST_SETUP; recursive=true)

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
    # buffer, the C/N0 estimators, the Doppler estimators — is TrackingLoops'
    # API, which users load next to Tracking and which Tracking does not
    # re-export. Its docstrings are still what these pages explain, and
    # Documenter only splices in docstrings from the listed modules, so both
    # are listed. That puts TrackingLoops' exports under `checkdocs` too; the
    # ones no topic page covers are collected in `trackingloops.md`.
    modules = [Tracking, TrackingLoops],
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
        "vector_tracking.md",
        "correlator.md",
        "cn0_estimator.md",
        "noise_estimator.md",
        "trackingloops.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaGNSS/Tracking.jl.git",
    push_preview = true,
)
