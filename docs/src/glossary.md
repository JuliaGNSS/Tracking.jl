# Glossary

Tracking.jl uses the terms of the
[TrackingLoops glossary](https://juliagnss.github.io/TrackingLoops.jl/stable/glossary/),
which ties them to the textbooks (Kaplan & Hegarty, *Understanding GPS: Principles and
Applications*, 2nd ed., 2006; Misra & Enge, *Global Positioning System*; Parkinson &
Spilker (eds.), *Global Positioning System: Theory and Applications*, 1996): an **integration** is a coherent integration (Kaplan and Hegarty's
predetection integration) over whole primary code periods, its length `T` the
*coherent integration time*; the **prompt** is the on-time correlator; the
**post-correlation filter** combines an antenna array's correlator outputs into the
*filtered* ones the discriminators read; the **driver** is the signal whose
integrations close a satellite's carrier and code loops, a **passenger** a signal
riding on its replica; `last_x` is `x` of the latest integration and `previous_x` that
of the one before. The terms below are this package's own.

- **Samples**: one buffer of complex baseband or IF samples, one column per antenna
  element, as one [`track!`](@ref) call takes it; a *chunk* in TrackingLoops'
  glossary. [`BandSamples`](@ref) wraps one band's samples with their sampling and
  intermediate frequencies; a multi-band call takes a
  [`SamplesByBand`](@ref Tracking.SamplesByBand), one per band id.
- **NCO update span**: the part of the samples between two NCO updates,
  `nco_update_interval` long (one primary code period by default): one correlate pass
  at a fixed NCO Doppler, then one loop update. `nco_update_index` counts the spans of
  a call (see [NCO updates](track.md#NCO-updates)).
- **Sample index**: a sample's position in the current samples, counted from 1;
  `sample_range` is a span of them. A satellite's
  [`get_next_sample_index`](@ref) is the first sample it has not yet correlated. The
  integrations' end sample indices run on across calls once the loop update adds the
  band's `sample_index_offsets` (see [`TrackState`](@ref)), the shared time grid the
  vector loop needs.
- **Correlator bank**: a signal's early, prompt and late (and more) correlators, the
  correlator outputs of one integration; a
  [`CorrelatorBankOutput`](@extref TrackingLoops.CorrelatorBankOutput) is one completed
  integration's, tagged with its sample count and end sample index.
- **Tracked satellite**: a [`TrackedSat`](@ref), the carrier and code NCO state a
  satellite's signals share, with one [`TrackedSignal`](@ref) per signal; its first
  signal is the driver.
- **Signal group**: a [`SignalGroup`](@ref), the satellites tracking one tuple of
  signals on one band and antenna array, addressed by its `group_key`.
- **Replica carrier phase**: [`get_carrier_phase`](@ref), the phase of the carrier
  replica in cycles, folded to `[-0.5, 0.5)`; not the accumulated carrier-phase
  measurement (RINEX `L`), for which a host adds
  [`calc_carrier_phase_correction`](@ref TrackingLoops.calc_carrier_phase_correction(::TrackState, ::Union{Integer, Symbol, Val}, ::Any))
  and counts the cycles itself.

## Sign conventions

Tracking reports its quantities with
[TrackingLoops' signs](https://juliagnss.github.io/TrackingLoops.jl/stable/glossary/#Sign-conventions):
[`get_carrier_doppler`](@ref) is positive for a closing satellite, the received carrier
above the intermediate frequency, and [`get_code_doppler`](@ref) has the same sign,
scaled by the chip rate's ratio to the carrier frequency. The replica carrier phase
([`get_carrier_phase`](@ref)) advances by the intermediate frequency plus the carrier
Doppler, so a closing satellite turns it forward, and the code phase ([`get_code_phase`](@ref)) by the chip
rate plus the code Doppler.
