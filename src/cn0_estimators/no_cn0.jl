"""
$(SIGNATURES)

A CN0 estimator that estimates nothing: it keeps no state, does no work per
record, and reports `-Inf dB-Hz`. Use it to say **"do not measure this signal's
C/N₀"** — either to skip the per-record work on a signal whose C/N₀ nobody reads,
or, more importantly, to avoid *publishing a number that cannot be trusted*.

Two places where that matters:

  - the non-driver signals of a multi-signal [`TrackedSat`](@ref) whose C/N₀ is
    never read: `cn0_estimator = (NWPRCN0Estimator(), NoCN0Estimator())` (the
    saving is small — the per-record update is well under a per cent of a
    correlation — so reach for this for the reason below, not for speed);
  - as [`NWPRCN0Estimator`](@ref)'s `fallback`, where it replaces the
    [`MomentsCN0Estimator`](@ref)'s ~27.6 dB-Hz noise floor with an honest "no
    estimate" for the signals and phases that admit no coherent window (a
    secondary-coded signal before sync, or GPS L1C-D / Galileo E1B at any time):
    `NWPRCN0Estimator(; fallback = NoCN0Estimator())`.

`-Inf dB-Hz` rather than `NaN dB-Hz`, even though "not measured" is what is meant:
`NaN dB-Hz >= threshold` is `true` for every threshold with Unitful's `Level`
comparison, so a `NaN` would clear every lock detector it met. `-Inf` compares
`false` against any finite threshold, which is the safe answer to "is this signal
locked?".
"""
struct NoCN0Estimator <: AbstractCN0Estimator end

"""
$(SIGNATURES)

Ignore the prompt and return the estimator unchanged — [`NoCN0Estimator`](@ref)
holds no state.
"""
update(estimator::NoCN0Estimator, prompt) = estimator

"""
$(SIGNATURES)

Always `-Inf dB-Hz`; see [`NoCN0Estimator`](@ref) for why that value and not
`NaN`.
"""
estimate_cn0(::NoCN0Estimator, integration_time) = dBHz(0.0 / integration_time)
