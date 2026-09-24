# Loop core reference (TrackingLoops.jl)

The per-record arithmetic of a tracking loop — correlators, discriminators, the
bit buffer and its sync detectors, the C/N₀ estimators, the noise window, the
Doppler estimators and the NCO timeline — lives in
[TrackingLoops.jl](https://github.com/JuliaGNSS/TrackingLoops.jl). Tracking.jl
builds on it and re-exports only the functions that read its results
(`estimate_cn0`, `get_prompt`, `get_early`, `get_late`, `get_accumulators`,
`get_soft_bits`, `has_bit_or_secondary_code_been_found`, `get_num_ants`,
`append_noise_observation!`). To configure a correlator or an estimator, load
TrackingLoops next to Tracking:

```julia
using Tracking, TrackingLoops, GNSSSignals
```

```@docs
TrackingLoops
```

The topic pages of this manual document the parts of TrackingLoops that
configure a `TrackState` — correlators, C/N₀ and noise estimators, Doppler
estimators. The rest of its API, which a custom estimator, an external
correlator producer or a hardware loop process reaches for, is listed here.

## Correlator taps

```@docs
TrackingLoops.get_very_early
TrackingLoops.get_very_late
TrackingLoops.get_prompt_index
TrackingLoops.normalize
```

## Discriminators and loop-filter steps

```@docs
TrackingLoops.pll_disc
TrackingLoops.fll_disc
TrackingLoops.dll_disc
TrackingLoops.aid_dopplers
TrackingLoops.calc_num_code_blocks_to_integrate
TrackingLoops.calc_num_code_blocks_for_bit_buffer
TrackingLoops.estimator_state_type
```

## The per-record fold

```@docs
TrackingLoops.update
TrackingLoops.SignalLoopState
TrackingLoops.apply_record
TrackingLoops.restart_bit_clock
TrackingLoops.reset_signal_state
```

## Bit and secondary-code synchronisation

```@docs
TrackingLoops.BitBuffer
TrackingLoops.PhaseAccumulators
TrackingLoops.SyncResult
TrackingLoops.detect_bit_or_secondary_code_sync
TrackingLoops.get_code_block_buffer_type
TrackingLoops.uses_soft_bit_edge_detection
TrackingLoops.uses_soft_secondary_code_detection
TrackingLoops.get_bit_edge_detection_confidence
TrackingLoops.get_bit_edge_or_secondary_code_tolerance
```

## Noise window

```@docs
TrackingLoops.NoiseEstimators
TrackingLoops.noise_window_looks
```

## NCO timeline

```@docs
TrackingLoops.NO_LANDING_SAMPLE
TrackingLoops.reset_timeline!
TrackingLoops.schedule_word!
TrackingLoops.promote_words!
TrackingLoops.reschedule_word!
TrackingLoops.word_changes_within
TrackingLoops.nco_word_at
```
