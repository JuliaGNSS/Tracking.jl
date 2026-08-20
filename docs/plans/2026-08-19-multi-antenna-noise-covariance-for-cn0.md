# Multi-antenna noise covariance for per-satellite C/N₀

*2026-08-19. Follows `2026-08-06-per-band-noise-estimation-for-cn0.md`, which shipped the
per-signal noise reference in v7.0.0 and left the antenna-array case explicitly deferred.*

## The problem

The prompt a C/N₀ estimator divides is **post**-`AbstractPostCorrFilter`; the noise floor it
divides by was measured **pre**-filter, on a single antenna column. `update_noise!` collapsed
the sample matrix through `_noise_reference_samples` — `view(samples, :, size(samples, 2))` —
and the ratio was only valid because `DefaultPostCorrFilter` happened to select that same
last column and have unity gain.

Any real beamformer broke it silently. An averaging combiner over `M` uncorrelated antennas
suppresses the noise by `1/M` while the signal adds coherently, so its true C/N₀ is
`10log₁₀(M)` dB above the single-antenna one — but the denominator did not move at all, so
the reported number was simply wrong, with nothing in the output to say so. The v7.0.0 docs
recorded this as a known gap ("the post-correlation filter is assumed to have unity noise
gain… there is no `noise_gain` hook yet").

Falling back to a self-contained estimator instead was considered and rejected.
`MomentsCN0Estimator` carries a ~27.6 dB-Hz noise-floor bias that makes it unusable below
that C/N₀; `NWPRCN0Estimator` needs bit sync and does not work on GPS L1C-D, Galileo E1B or
presync secondary-coded signals. Steering beamforming users onto those would strip exactly
the receivers that beamform to fight jamming of the one estimator that is good in that
regime — which is why the library's default moved from NWPR to `NoiseRefCN0Estimator` in the
first place.

## The design

**Measure a covariance, reduce it per satellite.** The noise reference despreads *every*
antenna column and pools the taps into the array's `M×M` spatial covariance
`R̂ = Σ_k b_k·b_kᴴ` instead of the scalar `Σ_k |b_k|²`. Its diagonal is each antenna's own
`N₀`; its off-diagonals are the antennas' noise correlation.

This stays **one window per signal**, shared by every satellite tracking it — exactly as
satellite-agnostic as the scalar version, and costing one despread per signal per chunk
rather than one per satellite. The per-satellite step happens at read time:

```
N₀ = wᴴ R̂ w
```

with `w` that satellite's own current combining weights. This is exact, not approximate, for
any fixed `w`, because `E[|wᴴn|²] = wᴴRw` — so no per-satellite noise state is stored
anywhere. Only the final read is per-satellite.

**`AbstractPostCorrFilter` now declares weights instead of being a callable.** `wᴴR̂w`
requires knowing `w`, so the filter's contract changed from "any callable applied to a tap"
to a required `get_weights(filter, ::NumAnts{M})`. Tracking combines the taps itself. The
contract is now deliberately **linear in the antennas**; a non-linear combiner has no closed
form for its own noise floor, and would be back to a silently wrong ratio. Nothing in the
repo needed one — both existing implementations (`DefaultPostCorrFilter` and the docs'
`MyBeamformer`) were already linear.

`get_default_post_corr_filter` was removed at the same time. It had no caller outside its own
test, and its multi-antenna method (`AbstractCorrelator{T} where T<:SVector`) could never
dispatch, because that type parameter is the antenna *count*, an `Int`.

**Hardware sources report per-antenna too.** All three public builders gained the
per-antenna form of the input they already take — `noise_observation` an `SVector` of
accumulations, `noise_observation_from_correlator` and `noise_observation_from_samples` a
pre-summed `M×M` covariance. An `M`-element front end has `M` numbers; making it collapse
them would recreate the same bug on the hardware path, since a scalar cannot answer `wᴴR̂w`.

**Antenna counts are checked at construction.** `TrackState` auto-provisions each
`CorrelatorNoiseEstimator` at its signal group's own `num_ants` (the key-walk now carries
`(signal_id, num_ants)` pairs rather than bare symbols). An explicitly-passed estimator whose
count disagrees throws an `ArgumentError` naming the fix, rather than failing as a shape
error inside the per-chunk fold. The check reads the count off `noise_density_type` alone, so
it holds for any `AbstractNoiseEstimator`, a custom hardware one included.

## Choices worth recording

**A bare `SMatrix`, not `Hermitian`.** The window's totals need only `zero`, `+` and `/`,
which StaticArrays supplies, and the result stays isbits — which is what keeps
`NoiseWindowTotals` inline in its `Ref` and the whole measurement allocation-free. A
`Hermitian` wrapper or a `Matrix` would cost that for nothing; `R̂` is never inverted or
factorised.

**"Positive" means a positive diagonal, not positive definiteness.**
`_noise_density_and_ready` rejects an unusable floor, and for a covariance the usable test is
that each antenna measured some power. A near-singular `R̂` is legitimate rather than faulty —
strongly correlated antennas (in the limit, one signal fed to every element) give a rank-1
covariance that is still the right answer for every `wᴴR̂w` anyone will ask of it. Only a
genuinely dead element matters, and that shows on the diagonal.

**But rank *deficiency from too few looks* is a different thing, and is gated.** The
reference pools `num_taps` — three — rank-1 outer products per observation, so a window
holding one observation spans at most three of `M` dimensions. That is full rank at `M ≤ 3`
and deficient above it, and it is deficient *by construction* rather than because of what
the samples were. Measured on a one-observation window at `M = 4`, a beamformer whose
weights lie near the unmeasured subspace reads a floor down to 2 % of the true one — a C/N₀
about 16 dB optimistic — and the diagonal stays positive throughout, so the check above
cannot see it. So `_noise_density_and_ready` also requires the window to hold `M` looks
(`noise_window_looks`, `_sufficient_looks`), which withholds the first `⌈M/3⌉` observations.
`M ≤ 3` is never gated, and a chunk spanning several code periods clears it within the first
call; only a caller streaming one-code-period buffers at `M ≥ 4` ever sees it, and then for
one or two chunks.

While gated, the density is withheld exactly as an empty window's is: the fold skips the
C/N₀ update and `estimate_cn0` reports `-Inf dB-Hz`, which every downstream consumer already
handles (GNSSReceiver's `CodeLockDetector` linearises it to `0 Hz` and needs 80 ms of warm-up
plus 200 ms of dwell before it drops a satellite, so a one-chunk transient is invisible to
it). The one thing that is *not* reused is the warning: a filling window is a normal
transient, not a misconfiguration, so `_noise_window_filling` keeps it quiet. It is phrased
as "everything else about this density is fine and only the look count is short" rather than
as "the window is non-empty", so a measured floor of zero — a dead input — still reaches the
diagnostic it is meant to reach.

**`M = 1` is bit-identical to before, by construction.** `DefaultPostCorrFilter`'s scalar
weight is exactly `1.0 + 0.0im`, so `conj(w)·tap` and `abs2(w)·N₀` are IEEE-exact no-ops. The
regression lock in `test/track.jl` is stronger than a tolerance: a three-antenna run and a
single-antenna run of the same last column are asserted to report the *identical* C/N₀.

**The combine step is its own function, and that is load-bearing.**
`_apply_correlator_output` is large enough that a closure capturing the weights is built on
the heap there — 256 B per record, which the `track!` allocation guards catch immediately.
`_combine_correlator` keeps the capture inside a function small enough to inline. The same
trap is why `_num_ants_val` (dispatching on the correlator's type parameter) exists rather
than `NumAnts(get_num_ants(c))`, which routes the count through a runtime `Int` and loses
inference; that helper moved from `downconvert_and_correlate_int16.jl` to
`correlators/correlator.jl` so both callers share one definition.

## What broke, and for whom

- Custom `AbstractPostCorrFilter`s must define `get_weights` and stop defining a call
  operator. Both doctest beamformers were migrated; the mean-of-antennas example they show is
  now *correct* rather than caveated.
- `get_default_post_corr_filter` is gone.
- A multi-antenna signal group with an explicitly-passed single-antenna noise estimator now
  throws at `TrackState` construction. That combination previously "worked" by silently
  measuring one column.
