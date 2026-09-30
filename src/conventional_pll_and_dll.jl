# Target `BL · Δt` for a loop update interval of `Δt` — ~10× margin from the
# `BL · Δt < 0.18` practical stability edge of the bilinear third-order carrier
# filter. Sizes the carrier default against the primary code period, and caps
# the code loop against each record's actual integration time (see
# `effective_code_loop_filter_bandwidth`). For the code loop the shared product
# is a conservative reuse: its filter is the *second*-order bilinear one, which
# destabilizes only around the classic `BL · Δt ≈ 0.4` of transform-designed
# digital loops (Stephens & Thomas 1995, "Controlled-Root Formulation for
# Digital Phase-Locked Loops", IEEE Trans. AES 31(1) — this implementation's
# linearized edge lands there numerically too), so the cap only adds margin.
const MAX_LOOP_BANDWIDTH_TIME_PRODUCT = 0.018

"""
$(SIGNATURES)

Recommended carrier-loop-filter bandwidth for `signal`'s primary integration
period. Sized so that the PLL time-bandwidth product `BL * T` lands at
about 0.018 (≈10× margin from the 0.18 stability edge of the bilinear
third-order filter). Used by [`TrackState(; signal=…)`](@ref) when the
user doesn't pass an explicit `doppler_estimator`.

Override by defining a method for your signal type, or by constructing
[`ConventionalAssistedPLLAndDLL`](@ref) yourself with explicit
`carrier_loop_filter_bandwidth =` / `code_loop_filter_bandwidth =` kwargs.

```julia
T = get_code_length(signal) / get_code_frequency(signal)   # primary period
BL = 0.018 / T                                              # this default
```

`T` here is the **primary**-code period, not the chosen coherent
integration length. For GPS L1 C/A (T = 1 ms) and GPS L5I (T = 1 ms, a
10230-chip code at 10.23 MHz) this returns 18 Hz — matching the historical
hand-picked default. For L1C-D / L1C-P (T = 10 ms) it returns 1.8 Hz, and
for Galileo E1B (T = 4 ms) 4.5 Hz — the well-inside-stability values the
multi-signal flagship use case needs.

This value is the **reference** bandwidth for a one-primary-code-period
integration; it is not the bandwidth that ends up in the loop when you
integrate longer. Coherently integrating `N` primary blocks grows the loop
update interval to `N·T`, which would push `BL·N·T` toward the ~0.18
stability edge of the bilinear filter. To avoid that, the conventional
estimator **automatically scales the effective loop bandwidth by
`1/N`** at filter time (see [`ConventionalPLLAndDLL`](@ref)), holding the
`BL·Δt` stability product fixed at its single-period value. So you set this
reference bandwidth once and the loop stays stable at any integration length
— no manual `1/N` adjustment is needed.
"""
function default_carrier_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    # T = the primary code period — one code block, not the chosen coherent
    # integration length. The estimator's bandwidth fields are typed
    # `typeof(1.0Hz)`, so explicitly land on Hz (otherwise `1/s` propagates and
    # trips the typed field assignment).
    primary_period = get_code_length(signal) / get_code_frequency(signal)
    uconvert(Hz, MAX_LOOP_BANDWIDTH_TIME_PRODUCT / primary_period)
end

"""
$(SIGNATURES)

Recommended code-loop-filter (DLL) bandwidth for `signal`: a flat 1 Hz for every
signal.

A *carrier-aided* DLL has almost no dynamic stress to track — the code Doppler
is handed to it by the PLL (see `aid_dopplers`) — so its bandwidth is a
thermal-noise-versus-pull-in trade that scales with neither the symbol rate
(the old 18:1 carrier:code ratio starved the long-primary signals' pull-in) nor
the coherent integration length. 1 Hz sits inside the 0.25–2 Hz the reference
software receivers (GNSS-SDR, SoftGNSS, PocketSDR) use across signals.

Unlike the carrier bandwidth this is an **absolute** value, not a
per-primary-code-period reference. Only the loop's own `BL · Δt` stability
product caps it, at filter time, against each record's actual integration time
— see [`effective_code_loop_filter_bandwidth`](@ref); the cap binds only past
18 ms (0.9 Hz for a 20 ms L2 CM integration, 0.012 Hz for a 1.5 s L2 CL one).

Override by defining a method for your signal type.
"""
function default_code_loop_filter_bandwidth(signal::AbstractGNSSSignal)
    1.0Hz
end

"""
$(SIGNATURES)

Effective code-loop bandwidth for a record that integrated for
`integration_time`: the configured bandwidth, capped so the code loop's
`BL · Δt` product stays inside `MAX_LOOP_BANDWIDTH_TIME_PRODUCT`.

The carrier loop takes a `1/N` scaling instead, because its configured
bandwidth is a per-primary-code-period *reference* — see
[`ConventionalPLLAndDLL`](@ref). The DLL's is an absolute value: carrier-aided,
it has no dynamic stress that grows with the integration length, and neither its
pull-in time nor its thermal-noise floor does either, so integrating longer must
not narrow it. Only stability may, and stability depends on the update interval
the record actually had — hence the cap against `integration_time` rather than a
scaling by the block count. A `1/N` here would take a 20 ms L1 C/A integration
down to 0.05 Hz where stability allows 0.9 Hz, re-introducing through the
integration length exactly the pull-in sag that sizing the DLL off the carrier
default used to cause by signal.

For a single-block integration of any signal at or below the 18 ms period where
the cap starts to bind, this returns the configured bandwidth unchanged.
"""
@inline function effective_code_loop_filter_bandwidth(bandwidth, integration_time)
    min(bandwidth, uconvert(Hz, MAX_LOOP_BANDWIDTH_TIME_PRODUCT / integration_time))
end

"""
Per-satellite state for the conventional PLL and DLL Doppler estimator.
Holds initial Doppler values and loop filter states.
"""
@kwdef struct SatConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
end

function SatConventionalPLLAndDLL(
    sat::TrackedSat,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

function SatConventionalPLLAndDLL(
    sat_conventional_pll_and_dll::SatConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatConventionalPLLAndDLL{CA,CO}(
        sat_conventional_pll_and_dll.init_carrier_doppler,
        sat_conventional_pll_and_dll.init_code_doppler,
        isnothing(carrier_loop_filter) ? sat_conventional_pll_and_dll.carrier_loop_filter :
        carrier_loop_filter,
        isnothing(code_loop_filter) ? sat_conventional_pll_and_dll.code_loop_filter :
        code_loop_filter,
        isnothing(carrier_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.carrier_loop_filter_bandwidth :
        carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        sat_conventional_pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Conventional Phase-Locked Loop (PLL) and Delay-Locked Loop (DLL) Doppler
estimator. Configuration-only — per-satellite state lives in each
[`TrackedSat`](@ref) wrapper, produced via [`init_estimator_state`](@ref).

Type parameters `CA` and `CO` select the carrier and code loop filter types;
the bandwidth fields configure the loop bandwidths used when seeding new
satellites. Each bandwidth field is `Maybe{typeof(1.0Hz)}`: a `nothing`
field (the default) means **auto** — [`init_estimator_state`](@ref) sizes the
bandwidth per satellite from that sat's estimator-driver signal (`signals[1]`)
via [`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref): the carrier loop is sized for the
signal's own integration period (18 Hz for GPS L1 C/A, 4.5 Hz for Galileo E1B,
1.8 Hz for L1C-D / L1C-P, …), the code loop takes a flat 1 Hz. Pass an explicit
bandwidth to override the auto-sizing for every satellite this estimator seeds.

The two bandwidths are referenced differently, because only the carrier loop's
tuning tracks the update rate. The **carrier** bandwidth is referenced to a
one-primary-code-period integration: when a signal coherently integrates `N`
primary blocks (its per-[`TrackedSignal`](@ref)
`preferred_num_code_blocks_to_integrate`, set via
[`set_preferred_num_code_blocks_to_integrate!`](@ref)), it is automatically
scaled to `BL/N` at filter time so the loop's `BL·Δt` stability product stays at
its single-period value. This keeps the loop stable across integration lengths
without the caller re-tuning the bandwidth — e.g. a 1 ms→10 ms switch needs no
bandwidth change. The **code** bandwidth is an absolute value that longer
integration does not narrow; it is only capped by the same stability product
against the record's actual integration time — see
[`effective_code_loop_filter_bandwidth`](@ref).

## Multi-signal discriminator combining

On a satellite tracking several signals, `signals[1]` is still the
estimator-driver signal — it sets the loop cadence, the bandwidths, the carrier
phase reference and the datum of the shared `code_phase` — but in a
[`SignalGroup`](@ref) built with `discriminator_combining = true` the
**passenger** signals `signals[2:end]` aid it: a passenger record that ended on
the same sample as a driver record (equal `sample_index`) is folded into that
driver record's loop update, as a weighted mean of the PLL, FLL and DLL
discriminators. Passenger records that coincide with no driver record are
applied to their own signal as always and reach no loop; nothing is carried
across driver records or chunks. A single-signal satellite, and any satellite of
a group with the flag off, is unaffected — bit-identically so.

Each record's weight is its nominal post-integration SNR,
`get_relative_power(signal) * integrated_samples`, over the loop's
discriminator noise gain: [`dll_disc_noise_gain`](@ref Tracking.dll_disc_noise_gain)
for the code loop, a signal-independent constant for the carrier phase loop,
and an extra `integrated_samples²` for the FLL. A passenger in carrier
quadrature with the driver is de-rotated onto its phase frame before the PLL
discriminator; a passenger's code discriminator is referred to the driver's code
phase with the two signals' [`get_group_delay`](@ref) values, and withheld from
the code loop while either is unknown.

Combining assumes the driver integrates at least as long as every passenger —
for full aiding, equally long, which every intra-band pilot/data pair has at the
default integration length. Nothing enforces it; see
[Multi-signal discriminator combining](@ref Multi-signal-discriminator-combining)
in the manual for what a mismatch costs, for why the weights are nominal rather
than measured, and for why the fusion happens at the discriminator level rather
than by summing prompts.
"""
struct ConventionalPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
end

function ConventionalPLLAndDLL(
    ::Type{CA} = ThirdOrderBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(carrier_loop_filter_bandwidth, code_loop_filter_bandwidth)
end

"""
$(SIGNATURES)

Create a ConventionalPLLAndDLL with FLL-assisted carrier tracking. This is the
default Doppler estimator used by TrackState. Uses a ThirdOrderAssistedBilinearLF
for the carrier loop filter which combines PLL and FLL discriminators for
improved tracking under high dynamics.

Bandwidths default to `nothing` (auto): each satellite is seeded with the
loop bandwidth recommended for its own estimator-driver signal — see
[`ConventionalPLLAndDLL`](@ref). Pass explicit bandwidths to override.
"""
function ConventionalAssistedPLLAndDLL(
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL(
        ThirdOrderAssistedBilinearLF,
        CO;
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Kwarg-update constructor for tweaking bandwidths in place.
function ConventionalPLLAndDLL(
    pll_and_dll::ConventionalPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    ConventionalPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Build the per-satellite estimator state stored in a [`TrackedSat`](@ref) for a
satellite tracked under [`ConventionalPLLAndDLL`](@ref).

Auto bandwidths (`nothing` on the estimator) are resolved here, per satellite,
from the sat's estimator-driver signal (`signals[1]`): each sat gets the loop
bandwidth recommended for the signal that actually drives its loop, so a
multi-group / multi-constellation [`TrackState`](@ref) ends up with the right
bandwidth per group even though it carries one shared estimator. An explicit
bandwidth on the estimator is used verbatim for every satellite.
"""
function init_estimator_state(
    estimator::ConventionalPLLAndDLL{CA,CO},
    sat::TrackedSat,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    carrier_loop_filter = constructorof(CA)()
    code_loop_filter = constructorof(CO)()
    driver_signal = first(sat.signals).signal
    carrier_loop_filter_bandwidth =
        isnothing(estimator.carrier_loop_filter_bandwidth) ?
        default_carrier_loop_filter_bandwidth(driver_signal) :
        estimator.carrier_loop_filter_bandwidth
    code_loop_filter_bandwidth =
        isnothing(estimator.code_loop_filter_bandwidth) ?
        default_code_loop_filter_bandwidth(driver_signal) :
        estimator.code_loop_filter_bandwidth
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Re-seed hook used by `reset_loop_filters!`. The generic fallback simply
# rebuilds the per-sat state from scratch via `init_estimator_state`; custom
# estimators may specialize to preserve per-sat configuration across the
# reset.
_reset_estimator_state(estimator::AbstractDopplerEstimator, sat::TrackedSat) =
    init_estimator_state(estimator, sat)

# Conventional PLL/DLL: zero the loop-filter integrators and re-seed the
# init Dopplers from the sat's current (converged) Dopplers, but keep the
# bandwidths from the EXISTING per-sat state — a per-sat
# `SatConventionalPLLAndDLL` bandwidth override must survive the reset
# (going through `init_estimator_state` would silently revert it to the
# estimator-level defaults).
function _reset_estimator_state(
    ::ConventionalPLLAndDLL,
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatConventionalPLLAndDLL},
)
    state = sat.doppler_estimator_state
    SatConventionalPLLAndDLL(
        sat.carrier_doppler,
        sat.code_doppler,
        constructorof(typeof(state.carrier_loop_filter))(),
        constructorof(typeof(state.code_loop_filter))(),
        state.carrier_loop_filter_bandwidth,
        state.code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(
    signal::AbstractGNSSSignal,
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
)
    carrier_doppler = carrier_freq_update
    code_doppler =
        code_freq_update + carrier_doppler * get_code_center_frequency_ratio(signal)
    init_carrier_doppler + carrier_doppler, init_code_doppler + code_doppler
end

# Per-sat update for the conventional PLL/DLL estimator. Pure: takes a
# TrackedSat and returns the updated TrackedSat. Shared by the immutable
# `estimate_dopplers_and_filter_prompt` and the in-place
# `estimate_dopplers_and_filter_prompt!` so the two cannot drift.
#
# `noise` is one `(density, ready)` pair per signal, in `sat.signals` order (see
# `_signal_noise_densities`) — so the driver takes `first(noise)` and the
# passengers `Base.tail(noise)`, and no signal is ever handed another's floor.
#
# `discriminator_combining` is the satellite's `SignalGroup` flag: whether the
# passengers' discriminators may be folded into the driver's loop update (see
# `_seed_and_fold`).
function _update_tracked_sat_doppler(
    sat::TrackedSat,
    sampling_frequency,
    noise::Tuple,
    discriminator_combining::Bool,
)
    # Walk all signals. For each one whose integration completed this
    # iteration, normalize/filter its prompt, advance CN0 and bit buffer,
    # and move its correlator to `last_fully_integrated_*`. `signals[1]` (the
    # estimator-driver signal) additionally runs the PLL/DLL and updates the
    # sat-shared carrier/code Doppler — fed, where the group combines, by the
    # passenger records that coincide with each driver record as well as by its
    # own (see `_fold_satellite_chunk`). Each signal's coherent-integration length
    # comes from its own `preferred_num_code_blocks_to_integrate`.
    state = sat.doppler_estimator_state
    driver = first(sat.signals)
    passengers = Base.tail(sat.signals)

    # The loops lock the driver (`signals[1]`) onto the real axis; every signal's
    # bit-buffer prompt is de-rotated by its carrier-phase offset from the driver
    # so a quadrature component (QPSK data/pilot, e.g. GPS L5 / Galileo E5a) does
    # not decode off the collapsed real part. The offset is the nominal per-signal
    # one from `get_carrier_phase_offset`, not the satellite's tracked
    # `carrier_phase`. A combined PLL de-rotates by the very same factor, for the
    # same reason.
    driver_carrier_phase_offset = get_carrier_phase_offset(driver.signal)

    new_driver,
    new_passengers,
    new_doppler_estimator_state,
    new_carrier_doppler,
    new_code_doppler = _seed_and_fold(
        driver,
        passengers,
        sat,
        state,
        sampling_frequency,
        noise,
        driver_carrier_phase_offset,
        discriminator_combining,
    )

    # Phase-snap fallback chain. Picks the synced signal with the
    # longest `(primary × secondary)` code length, and uses its
    # secondary-code phase to anchor `sat.code_phase` to the right
    # secondary-chip window.
    #
    # This is a *one-time* anchoring applied only on the iteration a
    # signal transitions `found == false → true`. It preserves the
    # within-primary-block phase (`mod(code_phase, primary)`) so the loop
    # keeps the current chunk-bounded position; re-running it on later
    # iterations would wedge the satellite (see issue #117). After sync,
    # `update`'s `mod(…, current_code_wrap)` maintains the alignment.
    #
    # Because sync is detected in this estimate pass — *after* the whole
    # chunk was correlated — any in-flight partial integration for this
    # chunk was accumulated at the pre-snap code phase (and, pre-sync, with
    # no secondary-code overlay). Once the snap jumps `code_phase` into the
    # secondary window that partial's data is phase-inconsistent with the
    # new alignment (for a flipped NH chip it can even cancel to zero and
    # feed a 0/0 into the discriminators). So on the snap we also reset every
    # signal's in-flight accumulator: the phase bookkeeping is kept, but the
    # next chunk re-integrates cleanly from the snapped phase to the next
    # boundary. Block-aligned starts have no residue, so this is a no-op there.
    new_signals = (new_driver, new_passengers...)
    just_synced = _any_signal_just_synced(sat.signals, new_signals)
    snapped_code_phase =
        just_synced ? _snap_code_phase_from_synced_signal(new_signals, sat.code_phase) :
        sat.code_phase
    final_signals =
        just_synced ? map(_reset_inflight_integration, new_signals) : new_signals

    TrackedSat(
        sat;
        code_phase = snapped_code_phase,
        carrier_doppler = new_carrier_doppler,
        code_doppler = new_code_doppler,
        signals = final_signals,
        doppler_estimator_state = new_doppler_estimator_state,
    )
end

# Drop an in-flight (partial) integration: zero the accumulator and its sample
# counter, leaving all other per-signal state intact. Used at the sync-transition
# phase snap, where the shared `code_phase` moves and any partial accumulated at
# the old phase must not be carried into the re-anchored window.
@inline _reset_inflight_integration(s::TrackedSignal) =
    TrackedSignal(s; correlator = zero(s.correlator), integrated_samples = 0)

# De-rotation applied to a component's bit-buffer prompt — and, when combining,
# to its correlator before the PLL discriminator — so its own energy is real
# again, given the loops lock the driver onto the real axis. The rotation is
# `cis(driver_carrier_phase_offset − get_carrier_phase_offset(signal))`. Both
# terms are the *nominal* per-signal carrier phase offset the ICD specifies
# (radians, relative to the band's in-phase reference), from
# `get_carrier_phase_offset` — a compile-time constant per signal, not the
# satellite's tracked `carrier_phase`, which cancels out of the difference
# because every signal of a group rides one carrier. Formed once per signal per
# chunk (see `RecordFoldContext`), never per record.
# For an in-phase component (co-phased with the driver, or the driver itself)
# the difference is 0 and `cis(0) === 1 + 0im`, a bit-identical no-op; a
# quadrature component (GPS L5 / Galileo E5a I-vs-Q) rotates by `±90°` onto the
# real axis, where the navigation decoder resolves the residual sign via its
# preamble.
@inline _carrier_phase_derotation(driver_carrier_phase_offset::Real, signal) =
    cis(driver_carrier_phase_offset - get_carrier_phase_offset(signal))

# The same de-rotation applied to every tap of a correlator, so the unchanged
# `pll_disc` reads a passenger's phase error in the driver's frame: a component
# in carrier quadrature with the driver (GPS L5I vs L5Q, Galileo E5aI vs E5aQ)
# would otherwise read `±π/2`. `dll_disc` and `fll_disc` are rotation-invariant
# and keep reading the un-rotated correlator. `nothing` is the driver's own
# record, which is not rotated at all rather than rotated by the identity.
# The call-site `@inline` is load-bearing — without it the closure escapes and
# the combining fold allocates per record.
@inline _derotate_correlator(correlator, derotation::Complex) =
    @inline apply(tap -> tap * derotation, correlator)
@inline _derotate_correlator(correlator, ::Nothing) = correlator

# Whether this record was correlated with a replica the sync detected *in this
# same chunk* had not yet corrected. `found_before_chunk` is the signal's sync
# state frozen at the chunk boundary (see `RecordFoldContext`), so the two
# together say "synced during this fold, and this record precedes the fix".
@inline _correlated_pre_sync(tracked_signal::TrackedSignal, found_before_chunk::Bool) =
    !found_before_chunk && has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer)

# Apply one completed `CorrelatorOutput` record to a signal — shared by the
# estimator-driver and passenger folds so they cannot drift (issue #133):
# normalize the record's (raw) correlator by its sample count, update/apply the
# post-corr filter, record the filtered prompt, advance the CN0 estimator and
# bit buffer, and rebuild the `TrackedSignal` with the record moved to
# `last_fully_integrated_*`. Returns the rebuilt signal plus the intermediate
# values the driver's loop-filter section needs (`filtered_correlator`,
# `integrated_code_blocks`).
#
# Unlike the old per-integration advance, this does NOT reset the live
# accumulator or `integrated_samples`: the correlate phase already reset them
# when it snapshotted this record and began (or is carrying) the next
# integration. It consumes only the record's stored correlator.
#
# The bit accumulator is credited with the blocks *actually* integrated
# (`calc_num_code_blocks_for_bit_buffer`), recovered from the record's sample
# count: post-sync the first integration is truncated to land on the data-bit
# boundary, so crediting the intended length would misalign the decoded bits
# (issue #125). The block count returned for the driver's `1/N` carrier-bandwidth
# scaling is the same actual count (floored at 1): the bandwidth must pair with
# the record's true integration time, so it only switches when the integration
# actually lengthened — not already on the fold where sync was detected but the
# records were still single-block (see `_advance_driver_record`).
# `correlated_pre_sync = true` marks a record that follows a bit/secondary sync
# detected earlier in the same fold, i.e. one that was correlated with a
# pre-sync replica. Its *prompt* is only unusable where sync changed the replica
# — the secondary-code wipe-off, whose absence would feed sign-corrupted prompts
# into the first post-sync bits — so it is dropped from the coherent bit
# accumulation for secondary-coded signals only; a signal without a secondary
# code (GPS L1 C/A) correlates identically either side of the sync instant and
# keeps its prompt. The code blocks such a record covers are real either way and
# are always credited to the accumulator's block count: dropping the count
# slides the bit window one block off the navigation-bit grid for the rest of
# the run, which costs ~0.9 dB of bit-decision SNR and makes every coherent
# window that follows the grid straddle a bit flip (issue #219).
@inline function _apply_correlator_output(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    prn::Integer,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
    derotation::Complex = one(ComplexF64);
    correlated_pre_sync::Bool = false,
)
    signal = tracked_signal.signal
    normalized_correlator =
        normalize(output.correlator, output.integrated_samples, get_code_amplitude(signal))
    post_corr_filter =
        update(tracked_signal.post_corr_filter, get_prompt(normalized_correlator))
    # The filter's weights, read once and used twice: to combine the antennas
    # here, and to reduce the shared noise covariance to *this* satellite's floor
    # below. Both sides of the C/N₀ ratio must go through the same `w` or they
    # describe different channels — which is precisely what went wrong when the
    # floor was measured on one fixed antenna and the prompt was beamformed.
    weights = get_weights(post_corr_filter, _num_ants_val(normalized_correlator))
    filtered_correlator = _combine_correlator(normalized_correlator, weights)
    prompt = get_prompt(filtered_correlator)
    push!(tracked_signal.filtered_prompts, prompt)
    bit_block_count = calc_num_code_blocks_for_bit_buffer(
        signal,
        output.integrated_samples,
        sampling_frequency,
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    )
    # Blocks this record actually covered, for the driver's `1/N` carrier-bandwidth
    # scaling. The floor at 1 covers the fractional-block record right after a
    # sync phase-snap accumulator reset, whose rounded block count can be 0.
    integrated_code_blocks = max(1, bit_block_count)
    # De-rotate the prompt onto the driver's (real) phase frame before both the
    # secondary/bit sync search and the coherent bit accumulation inside
    # `buffer`, so a quadrature component's data lands on the real axis it is
    # decided on. No-op for the driver and for co-phased pairs. `derotation` is
    # `_carrier_phase_derotation`'s value, formed once per signal per chunk by
    # the caller rather than per record: it depends only on the two nominal
    # carrier phase offsets, both compile-time constants per signal.
    bit_prompt = prompt * derotation
    # Keep a pre-sync-correlated record's prompt out of the coherent sum where
    # the sync changed the replica under it (secondary-code wipe-off), but
    # always let it advance the accumulator's block count — see above.
    drop_prompt = correlated_pre_sync && get_secondary_code_length(signal) > 1
    # The CN0 estimator is handed the navigation-bit state along with the prompt
    # (`CN0UpdateContext`, built from the bit buffer as it stands *before* this
    # record): `NWPRCN0Estimator` needs to know where the data-bit
    # boundaries are to sum prompts coherently over exactly one bit — nothing a
    # downstream consumer of `get_filtered_prompts` could reconstruct. The bit
    # grid is trustworthy for exactly the records whose prompt is: this record's
    # blocks are always credited to the accumulator, so a pre-sync-correlated
    # record is only unusable where the sync changed the replica under it, and
    # `drop_prompt` is that condition. Where it holds the context reports "no bit
    # grid", which keeps the sign-corrupted prompt out of any coherent window and
    # drops the window that was open — exactly what happens before sync.
    #
    # It also carries this signal's own noise density and this record's own
    # integration time, for an estimator that divides by a *measured* floor
    # (`NoiseRefCN0Estimator`). `noise_density_ready == false` means a source is
    # configured but its window is still empty, and then the update is skipped —
    # but only for the estimators that would actually read the density. Skipping
    # the whole record would corrupt a co-resident `NWPRCN0Estimator`: a record
    # missing from the bit grid makes `_update_nwpr` drop its open narrowband
    # window, silently demoting NWPR to its fallback for a whole `num_records`.
    # `requires_noise_density` is a compile-time constant on the estimator's
    # type, so the gate costs nothing at run time.
    #
    # The density arrives as whatever the signal's window measures — a scalar for
    # a single antenna, a spatial covariance for an array — and is reduced here,
    # through the same weights that produced `prompt`, to the one scalar floor
    # this satellite's combiner actually sees. `nothing` (no estimator configured)
    # passes straight through, so the context's type parameter stays `Nothing` and
    # the wiring mistake still surfaces at the first record.
    scalar_noise_density = _reduce_noise_density(noise_density, weights)
    cn0_estimator = _update_cn0_estimator(
        get_cn0_estimator(tracked_signal),
        prompt,
        signal,
        tracked_signal.bit_buffer,
        bit_block_count,
        !drop_prompt,
        scalar_noise_density,
        noise_density_ready,
        output.integrated_samples / sampling_frequency,
    )
    bit_buffer = buffer(
        signal,
        prn,
        tracked_signal.bit_buffer,
        bit_block_count,
        drop_prompt ? zero(bit_prompt) : bit_prompt,
    )
    # Such a record also moves the secondary-code anchor: the code-phase snap
    # runs after this fold and aligns the *upcoming* integration to
    # `bit_buffer.secondary_phase`, which the detector reported for the block
    # right after the syncing record.
    if correlated_pre_sync
        bit_buffer = _advance_secondary_phase(signal, bit_buffer, bit_block_count)
    end
    new_signal = TrackedSignal(
        tracked_signal;
        last_fully_integrated_filtered_prompt = prompt,
        bit_buffer,
        cn0_estimator,
        post_corr_filter,
        last_fully_integrated_correlator = output.correlator,
        # What the CN0 estimator's newest prompt was integrated over.
        last_fully_integrated_num_code_blocks = integrated_code_blocks,
    )
    return new_signal, filtered_correlator, integrated_code_blocks
end

# Build the context and fold the record into one signal's CN0 estimator, or skip
# it where the estimator needs a density that signal has not measured yet. Both
# branches return the same concrete estimator type, so this stays type-stable
# and allocation-free; the `requires_noise_density` half of the condition folds
# away at compile time.
@inline function _update_cn0_estimator(
    estimator::AbstractCN0Estimator,
    prompt,
    signal::AbstractGNSSSignal,
    bit_buffer::BitBuffer,
    num_code_blocks::Integer,
    bit_sync_usable::Bool,
    noise_density,
    noise_density_ready::Bool,
    integration_time,
)
    if !noise_density_ready && requires_noise_density(estimator)
        return estimator
    end
    # Positional, not keyword: this is the per-record site, and a keyword call
    # allocates per call on Julia 1.10 (see `CN0UpdateContext`).
    update(
        estimator,
        prompt,
        CN0UpdateContext(
            signal,
            bit_buffer,
            num_code_blocks,
            bit_sync_usable,
            noise_density,
            integration_time,
        ),
    )
end

# ---------------------------------------------------------------------------
# Per-signal chunk constants
# ---------------------------------------------------------------------------

"""
Everything **one signal's** fold needs that is fixed for the whole chunk, so the
recursive passenger walk threads one value per passenger instead of several
parallel tuples, and `signals`/`cursors` are the only state it rewrites. One of
these per signal, built once per chunk — the driver's included.

  - `signal_index` is this signal's position in `sat.signals` — `1` for the
    driver — for a measurement collector that keeps a per-signal quantity.
  - `code_doppler` is chunk-fixed: the Doppler that actually generated this
    chunk's replicas, which is what every `dll_disc` in the chunk must be
    evaluated at, driver and passengers alike.
  - `derotation` is `_carrier_phase_derotation`'s value for this signal against
    the driver — a function of two compile-time-constant carrier phase offsets
    and nothing else, so its `sincos` runs once per signal per chunk instead of
    once per record. The identity for the driver, which de-rotates against
    itself.
  - `group_delay_chips` is how far this signal's code phase sits from the
    driver's, in chips: `signals[1]`'s [`get_group_delay`](@ref) minus this
    signal's, at the code frequency in effect. `nothing` where either end is
    unknown, which withholds this signal from the combined **code** loop and
    leaves its carrier contributions untouched. A structural `0.0` for the
    driver, whose code phase *is* the reference — never `nothing`, so the driver
    stays in its own code loop even where its own datum is unknown.
  - `found_before_chunk` is the signal's sync state at the *chunk* boundary,
    frozen here rather than re-read during the walk: `correlated_pre_sync` must
    stay true for every record of this chunk that follows a mid-chunk sync
    detection, and re-reading the (now synced) bit buffer later would silently
    clear it.
  - `noise_density` / `noise_density_ready` are this signal's own measured floor
    and whether it is usable (see `_signal_noise_densities`) — per signal, never
    per band.
"""
struct RecordFoldContext{FS,D}
    prn::Int
    signal_index::Int
    sampling_frequency::FS
    code_doppler::typeof(1.0Hz)
    derotation::ComplexF64
    group_delay_chips::Maybe{Float64}
    found_before_chunk::Bool
    noise_density::D
    noise_density_ready::Bool
end

# Build one passenger's context from the passenger itself, so the values derived
# from the signal (`derotation`, `group_delay_chips`, `found_before_chunk`) are
# derived in exactly one place.
@inline function RecordFoldContext(
    tracked_signal::TrackedSignal,
    prn::Integer,
    signal_index::Integer,
    sampling_frequency,
    code_doppler,
    driver_carrier_phase_offset::Real,
    group_delay_datum::Maybe{typeof(1.0s)},
    noise::Tuple{Any,Bool},
)
    signal = tracked_signal.signal
    delay = _driver_relative_group_delay(get_group_delay(tracked_signal), group_delay_datum)
    RecordFoldContext(
        Int(prn),
        Int(signal_index),
        sampling_frequency,
        code_doppler,
        _carrier_phase_derotation(driver_carrier_phase_offset, signal),
        isnothing(delay) ? nothing :
        _group_delay_to_chips(delay, code_doppler + get_code_frequency(signal)),
        has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
        noise[1],
        noise[2],
    )
end

# The driver's own context. The values a passenger derives from its signal are
# structural here: the identity rotation, because the driver de-rotates against
# itself, and a `0.0` chip offset, because its code phase is the datum the others
# are referred to.
@inline _driver_fold_context(
    tracked_signal::TrackedSignal,
    prn::Integer,
    sampling_frequency,
    code_doppler,
    noise::Tuple{Any,Bool},
) = RecordFoldContext(
    Int(prn),
    1,
    sampling_frequency,
    code_doppler,
    one(ComplexF64),
    0.0,
    has_bit_or_secondary_code_been_found(tracked_signal.bit_buffer),
    noise[1],
    noise[2],
)

# One context per passenger, in `sat.signals[2:end]` order — paired positionally
# with `Base.tail(noise)`, which `_signal_noise_densities` builds in that same
# order. Heterogeneous signal types fold at inference time, so this stays
# type-stable and allocation-free. `signal_index` counts in `sat.signals`
# coordinates, so the first passenger is 2.
@inline _passenger_fold_contexts(
    ::Tuple{},
    ::Tuple{},
    ::Integer,
    sampling_frequency,
    code_doppler,
    ::Real,
    ::Maybe{typeof(1.0s)},
    signal_index::Int = 2,
) = ()
@inline function _passenger_fold_contexts(
    passengers::Tuple,
    noise::Tuple,
    prn::Integer,
    sampling_frequency,
    code_doppler,
    driver_carrier_phase_offset::Real,
    group_delay_datum::Maybe{typeof(1.0s)},
    signal_index::Int = 2,
)
    (
        RecordFoldContext(
            first(passengers),
            prn,
            signal_index,
            sampling_frequency,
            code_doppler,
            driver_carrier_phase_offset,
            group_delay_datum,
            first(noise),
        ),
        _passenger_fold_contexts(
            Base.tail(passengers),
            Base.tail(noise),
            prn,
            sampling_frequency,
            code_doppler,
            driver_carrier_phase_offset,
            group_delay_datum,
            signal_index + 1,
        )...,
    )
end

# ---------------------------------------------------------------------------
# Combining arithmetic
# ---------------------------------------------------------------------------

"""
Weighted sums of the discriminators contributed by the passenger records that
coincide with **one** driver record — built fresh for that record and consumed
by it, never stored: combining keeps no state across driver records or chunks.

One `(weight, sum)` pair per loop, because the three loops weight a record
differently (see `_weighted_record_contribution`).
"""
struct DiscriminatorAccumulator
    dll_weight::Float64
    dll_sum::Float64                # chips × weight
    pll_weight::Float64
    pll_sum::Float64                # radians × weight
    fll_weight::Float64
    fll_sum::typeof(1.0Hz)          # Hz × weight
end

DiscriminatorAccumulator() = DiscriminatorAccumulator(0.0, 0.0, 0.0, 0.0, 0.0, 0.0Hz)

# Fold one record's `_weighted_record_contribution` in, taking it by name so no
# caller can pair a value with another loop's weight.
@inline _accumulate_all_loops(acc::DiscriminatorAccumulator, contribution::NamedTuple) =
    DiscriminatorAccumulator(
        acc.dll_weight + contribution.dll_weight,
        acc.dll_sum + contribution.dll_weight * contribution.dll,
        acc.pll_weight + contribution.pll_weight,
        acc.pll_sum + contribution.pll_weight * contribution.pll,
        acc.fll_weight + contribution.fll_weight,
        acc.fll_sum + contribution.fll_weight * contribution.fll,
    )

# Two passengers' sums, for the tuple walk.
@inline _accumulate(a::DiscriminatorAccumulator, b::DiscriminatorAccumulator) =
    DiscriminatorAccumulator(
        a.dll_weight + b.dll_weight,
        a.dll_sum + b.dll_sum,
        a.pll_weight + b.pll_weight,
        a.pll_sum + b.pll_sum,
        a.fll_weight + b.fll_weight,
        a.fll_sum + b.fll_sum,
    )

"""
$(SIGNATURES)

Weight-normalized mean of one loop's accumulated `(weight, sum)`. A mean rather
than a sum keeps the loop gain independent of how many signals reported at this
epoch.

Zero weight answers a zero of the sum's own type rather than `0/0`. That is not a
degenerate case: the FLL's weight is zero on a signal's very first record, where
there is no previous prompt to difference against, and the code loop's is zero for
a record whose group delay could not be referred.
"""
@inline _normalized_mean(weighted_sum, weight) =
    iszero(weight) ? zero(weighted_sum) : weighted_sum / weight

# One loop's combined discriminator: the driver's own value and weight, plus
# whatever the coincident passengers contributed to that loop. With no passenger
# weight this is the driver's value *as is* — not `(w·d)/w`, which is not `d`
# bit for bit — so a driver record that no passenger record coincides with closes
# the loops exactly as the driver-only path does.
@inline _combined_discriminator(own, own_weight, passenger_sum, passenger_weight) =
    iszero(passenger_weight) ? own :
    _normalized_mean(own_weight * own + passenger_sum, own_weight + passenger_weight)

"""
$(SIGNATURES)

Everything one record contributes to the three combined loops, as a named
`(; dll, dll_weight, pll, pll_weight, fll, fll_weight)` — formed for the driver's
record and for every coincident passenger record alike, so they cannot end up on
different scales. Named rather than positional because two of the three values
are bare `Float64`s: a swapped pair would type-check, compile and quietly weight
one loop's error with another's.

Weights are nominal post-integration SNR (`_relative_record_snr`, the power
share `P` times the sample count `N`) over the loop's discriminator noise gain:

  - code loop: `P · N / dll_disc_noise_gain(…)`;
  - carrier phase loop: `P · N` — every `pll_disc` is the same Costas
    discriminator with a signal-independent noise gain, which cancels;
  - carrier frequency loop: `P · N³` — `fll_disc` divides an inter-prompt phase
    difference by `2π·T`, which scales its noise by `1/T`.

Only *ratios* within one satellite matter, so every factor shared across a group
cancels and is never formed — including the sampling frequency, which is why the
FLL's `T²` appears as `N²`.

`derotation` rotates the correlator onto the driver's phase frame before the PLL
discriminator (`_derotate_correlator`); `nothing` for the driver itself.
`group_delay_chips` is how far this signal's code phase sits from the driver's
(`RecordFoldContext.group_delay_chips`): a structural `0.0` for the driver, and
`nothing` for a signal whose datum or own value the caller has not supplied —
which zeroes the **code** weight only. A zeroed previous prompt (a signal's
first record, or the first after [`reset_loop_filters!`](@ref)) is "no frequency
measurement yet" and zeroes the FLL weight.
"""
@inline function _weighted_record_contribution(
    signal::AbstractGNSSSignal,
    filtered_correlator,
    integrated_samples::Integer,
    previous_prompt,
    derotation::Maybe{Complex},
    group_delay_chips::Maybe{Float64},
    code_doppler,
    sampling_frequency,
)
    relative_snr = _relative_record_snr(signal, integrated_samples)

    # --- Code loop -------------------------------------------------------
    # A passenger's code phase differs from the driver's by the satellite's
    # payload group delay difference, so its raw discriminator measures the
    # driver's error *plus* that offset. Subtract it to refer the measurement to
    # the driver's code phase, which is what the shared `sat.code_phase` (and
    # everything downstream of it) means. Unknown is not zero: see
    # `set_group_delay!` for why guessing it is the unsafe direction.
    dll, dll_weight = if !isnothing(group_delay_chips)
        raw = dll_disc(signal, filtered_correlator, code_doppler, sampling_frequency)
        gain = dll_disc_noise_gain(
            signal,
            filtered_correlator,
            code_doppler,
            sampling_frequency,
        )
        (raw - group_delay_chips, relative_snr / gain)
    else
        (0.0, 0.0)
    end

    # --- Carrier phase loop ----------------------------------------------
    # No group-delay analogue: the carrier carries no inter-signal group delay,
    # so a passenger contributes from its very first integration. The Costas
    # discriminator is blind to the ±1 data sign, which is what lets a data
    # component's value be averaged with a pilot's.
    pll = pll_disc(signal, _derotate_correlator(filtered_correlator, derotation))
    pll_weight = relative_snr

    # --- Carrier frequency (FLL) loop ------------------------------------
    fll, fll_weight = if iszero(previous_prompt)
        (0.0Hz, 0.0)
    else
        integration_time = integrated_samples / sampling_frequency
        (
            fll_disc(signal, filtered_correlator, previous_prompt, integration_time),
            relative_snr * integrated_samples^2,
        )
    end

    (; dll, dll_weight, pll, pll_weight, fll, fll_weight)
end

# A passenger record `_apply_correlator_output` marked `correlated_pre_sync` on a
# secondary-coded signal is excluded from combining outright rather than
# weighted down: it was correlated without the secondary-code wipe-off the sync
# had just established, so a multi-block record's coherent sum partially cancels
# and even a single-block one carries an unresolved secondary-chip sign. Exactly
# the condition under which that function drops the prompt from the bit
# accumulation, and for the same reason. The driver's own record is never
# excluded — it is the record whose loop update is being formed.
@inline _record_invalidated_by_sync(signal::AbstractGNSSSignal, correlated_pre_sync::Bool) =
    correlated_pre_sync && get_secondary_code_length(signal) > 1

# ---------------------------------------------------------------------------
# Estimator hooks
# ---------------------------------------------------------------------------
#
# The traversal below is shared by both shipped estimators. What differs between
# them is dispatched on the per-satellite state (or on the value a hook itself
# returned), so each estimator's rules stay explicit next to its state type:
#
#   - `_with_measurements(f, state)` — calls `f` with whatever the estimator
#     collects per record besides the loop update, `nothing` for this one. A
#     continuation rather than a return value, so a runtime choice between two
#     collector types becomes two specialisations instead of a `Union`.
#   - `_collect_measurements(measurements, signal, output, filtered_correlator,
#     previous_prompt, ctx)` — fold one record, of any signal, into it.
#   - `_close_loops(state, …)` — the per-record loop closure, returning the two
#     NCO updates plus the two threaded filters.
#   - `_park_state(state, carrier_loop_filter, code_loop_filter, measurements)` —
#     the per-satellite state the chunk leaves behind.

@inline _with_measurements(f::F, ::SatConventionalPLLAndDLL) where {F} = f(nothing)

@inline _collect_measurements(::Nothing, _signal, _output, _filtered, _previous, _ctx) =
    nothing

@inline _park_state(
    state::SatConventionalPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    ::Nothing,
) = SatConventionalPLLAndDLL(state; carrier_loop_filter, code_loop_filter)

# The plain FLL-assisted scalar PLL/DLL: both loops local, both discriminators
# read. Written without the per-satellite state, because it needs nothing from
# it — which is what lets an estimator whose own closure only *overrides* some
# records fall back to this one by calling it.
@inline function _close_scalar_loops(
    carrier_loop_filter,
    code_loop_filter,
    discriminators,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
)
    carrier_freq_update, carrier_loop_filter = _filter_carrier_loop(
        carrier_loop_filter,
        discriminators.pll,
        discriminators.fll,
        integration_time,
        carrier_bandwidth,
    )
    code_freq_update, code_loop_filter =
        filter_loop(code_loop_filter, discriminators.dll, integration_time, code_bandwidth)
    (; carrier_freq_update, code_freq_update, carrier_loop_filter, code_loop_filter)
end

# The default loop closure, the conventional estimator's.
@inline _close_loops(
    _state,
    carrier_loop_filter,
    code_loop_filter,
    discriminators,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
) = _close_scalar_loops(
    carrier_loop_filter,
    code_loop_filter,
    discriminators,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
)

# Push already-formed carrier discriminators through the carrier loop filter.
# This is the single place that knows an FLL-assisted filter takes a
# `(pll, fll)` pair while every other filter takes the PLL discriminator alone.
# The second argument of the pair is whatever the caller's closure feeds the FLL
# branch: the FLL discriminator in scalar tracking, an externally supplied
# frequency correction where another component closes that loop.
@inline _filter_carrier_loop(
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    pll_discriminator,
    fll_input,
    integration_time,
    loop_bandwidth,
) = filter_loop(
    carrier_loop_filter,
    (pll_discriminator, fll_input),
    integration_time,
    loop_bandwidth,
)

@inline _filter_carrier_loop(
    carrier_loop_filter::AbstractLoopFilter,
    pll_discriminator,
    _fll_input,
    integration_time,
    loop_bandwidth,
) = filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)

# ---------------------------------------------------------------------------
# The chunk traversal
# ---------------------------------------------------------------------------

# Everything every estimator does identically to one of the driver's
# `CorrelatorOutput`s before the loops are closed on it. Returns the rebuilt
# signal plus what loop closure needs, as a NamedTuple so positional results do
# not have to be unpacked in the right order.
@inline function _advance_driver_record(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    ctx::RecordFoldContext,
    carrier_loop_filter_bandwidth,
    code_loop_filter_bandwidth,
)
    sampling_frequency = ctx.sampling_frequency
    # FLL needs the previous record's filtered prompt; the first record of the
    # chunk chains from the sat's carried-over
    # `last_fully_integrated_filtered_prompt` (the previous chunk's last). Read
    # it off the signal BEFORE the advance overwrites it.
    previous_prompt = get_last_fully_integrated_filtered_prompt(tracked_signal)
    # Per-record integration time — the block time, NOT the chunk time.
    integration_time = output.integrated_samples / sampling_frequency
    # A record that follows a sync detected earlier in THIS fold was correlated
    # with pre-sync replicas — its blocks still count towards the bit, only its
    # prompt may have to be dropped (see `_apply_correlator_output`).
    correlated_pre_sync = _correlated_pre_sync(tracked_signal, ctx.found_before_chunk)
    tsig, filtered_correlator, integrated_code_blocks = _apply_correlator_output(
        tracked_signal,
        output,
        ctx.prn,
        sampling_frequency,
        ctx.noise_density,
        ctx.noise_density_ready,
        ctx.derotation;
        correlated_pre_sync,
    )

    # The configured CARRIER bandwidth is referenced to a
    # one-primary-code-period integration. Coherently integrating N periods
    # grows the loop update interval by that factor, so scale the effective
    # bandwidth by 1/N to hold the loop's BL·Δt stability product at its
    # single-period value. N (`integrated_code_blocks`) is the number of blocks
    # this record ACTUALLY covered — recovered from its sample count — not the
    # intended integration length: the bandwidth pairs with the record's true
    # `integration_time`, so it only switches once the integrations really
    # lengthen. Scaling by the intended length instead would under-gain the loop
    # for single-block records folded after a mid-fold sync detection
    # (correlated pre-sync, but the live bit buffer already reports the post-sync
    # length) and for the first post-sync integration, which is truncated to land
    # on the data-bit boundary. For the N=1 path this divides by 1 and is
    # bit-identical to before.
    carrier_bandwidth = carrier_loop_filter_bandwidth / integrated_code_blocks
    # The DLL's is an absolute bandwidth, so it is capped by its own stability
    # product against this record's integration time instead of scaled by N —
    # see `effective_code_loop_filter_bandwidth`.
    code_bandwidth =
        effective_code_loop_filter_bandwidth(code_loop_filter_bandwidth, integration_time)
    (;
        tracked_signal = tsig,
        filtered_correlator,
        previous_prompt,
        integration_time,
        carrier_bandwidth,
        code_bandwidth,
    )
end

# Apply one passenger record to its signal (prompt / CN0 / bit buffer) and hand
# it to the estimator's measurement collector. Returns the rebuilt signal, what
# a combining caller needs to form the record's contribution, and whether the
# pre-sync gate excludes the record — from combining and from the collector
# alike.
@inline function _apply_passenger_record(
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    ctx::RecordFoldContext,
    measurements,
)
    previous_prompt = get_last_fully_integrated_filtered_prompt(tracked_signal)
    correlated_pre_sync = _correlated_pre_sync(tracked_signal, ctx.found_before_chunk)
    tsig, filtered_correlator, _ = _apply_correlator_output(
        tracked_signal,
        output,
        ctx.prn,
        ctx.sampling_frequency,
        ctx.noise_density,
        ctx.noise_density_ready,
        ctx.derotation;
        correlated_pre_sync,
    )
    excluded = _record_invalidated_by_sync(tsig.signal, correlated_pre_sync)
    measurements =
        excluded ? measurements :
        _collect_measurements(
            measurements,
            tsig.signal,
            output,
            filtered_correlator,
            previous_prompt,
            ctx,
        )
    (tsig, filtered_correlator, previous_prompt, excluded, measurements)
end

# Advance one passenger up to the driver record ending at `sample_index`: every
# record that ended *before* it is applied to the signal and nothing more; a
# record ending on exactly that sample is applied and its contribution returned
# for the driver's loop update. A record ending later is left for the next
# driver record, or for the drain. Returns the rebuilt signal, the new cursor,
# the contribution (zero when nothing coincided or the record was excluded) and
# the grown measurements.
@inline function _advance_one_passenger_to(
    tracked_signal::TrackedSignal,
    cursor::Int,
    ctx::RecordFoldContext,
    sample_index::Int,
    measurements,
)
    tsig = tracked_signal
    outputs = tsig.correlator_outputs
    @inbounds while cursor <= length(outputs) && outputs[cursor].sample_index < sample_index
        tsig, _, _, _, measurements =
            _apply_passenger_record(tsig, outputs[cursor], ctx, measurements)
        cursor += 1
    end
    contribution = DiscriminatorAccumulator()
    @inbounds if cursor <= length(outputs) && outputs[cursor].sample_index == sample_index
        output = outputs[cursor]
        tsig, filtered_correlator, previous_prompt, excluded, measurements =
            _apply_passenger_record(tsig, output, ctx, measurements)
        cursor += 1
        if !excluded
            contribution = _accumulate_all_loops(
                contribution,
                _weighted_record_contribution(
                    tsig.signal,
                    filtered_correlator,
                    output.integrated_samples,
                    previous_prompt,
                    ctx.derotation,
                    ctx.group_delay_chips,
                    ctx.code_doppler,
                    ctx.sampling_frequency,
                ),
            )
        end
    end
    (tsig, cursor, contribution, measurements)
end

# Tuple-recursive walk of `_advance_one_passenger_to` over all passengers, in
# lockstep with their cursors and their `RecordFoldContext`s, summing the
# contributions. Heterogeneous signal types fold at inference time, so this stays
# type-stable and allocation-free.
#
# `Val(false)` — the group does not combine — walks nothing and contributes
# nothing: the passengers are left entirely to the drain, which applies them
# after the driver exactly as a driver-only fold does.
@inline _advance_passengers_to(
    ::Val{false},
    passengers::Tuple,
    cursors::Tuple,
    ::Tuple,
    ::Int,
    measurements,
) = (passengers, cursors, nothing, measurements)
@inline _advance_passengers_to(
    ::Val{true},
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Int,
    measurements,
) = ((), (), DiscriminatorAccumulator(), measurements)
@inline function _advance_passengers_to(
    combine::Val{true},
    passengers::Tuple,
    cursors::Tuple,
    contexts::Tuple,
    sample_index::Int,
    measurements,
)
    tsig, cursor, contribution, measurements = _advance_one_passenger_to(
        first(passengers),
        first(cursors),
        first(contexts),
        sample_index,
        measurements,
    )
    rest, rest_cursors, acc, measurements = _advance_passengers_to(
        combine,
        Base.tail(passengers),
        Base.tail(cursors),
        Base.tail(contexts),
        sample_index,
        measurements,
    )
    (
        (tsig, rest...),
        (cursor, rest_cursors...),
        _accumulate(contribution, acc),
        measurements,
    )
end

# Apply every record still unconsumed — the ones that ended after the chunk's
# last driver record, or all of them where the group does not combine — and
# empty each passenger's buffer. Their discriminators reach no loop.
@inline _drain_passengers(::Tuple{}, ::Tuple{}, ::Tuple{}, measurements) =
    ((), measurements)
@inline function _drain_passengers(
    passengers::Tuple,
    cursors::Tuple,
    contexts::Tuple,
    measurements,
)
    tsig = first(passengers)
    outputs = tsig.correlator_outputs
    ctx = first(contexts)
    @inbounds for k = first(cursors):length(outputs)
        tsig, _, _, _, measurements =
            _apply_passenger_record(tsig, outputs[k], ctx, measurements)
    end
    empty!(outputs)
    rest, measurements = _drain_passengers(
        Base.tail(passengers),
        Base.tail(cursors),
        Base.tail(contexts),
        measurements,
    )
    ((tsig, rest...), measurements)
end

# The loop discriminators for one driver record. `nothing` — the group does not
# combine, or the satellite has no passengers — is the driver's own three,
# exactly what the driver-only path has always computed, hence bit-identical
# Dopplers.
@inline _loop_discriminators(
    ::Nothing,
    signal::AbstractGNSSSignal,
    carrier_loop_filter::AbstractLoopFilter,
    folded,
    ::Integer,
    ctx::RecordFoldContext,
) = (;
    pll = pll_disc(signal, folded.filtered_correlator),
    fll = _driver_fll_discriminator(carrier_loop_filter, signal, folded),
    dll = dll_disc(
        signal,
        folded.filtered_correlator,
        ctx.code_doppler,
        ctx.sampling_frequency,
    ),
)

# With nothing combined the FLL discriminator's one consumer is the carrier loop
# filter's FLL branch, and only the FLL-assisted filter has one — so on every
# other filter this answers `nothing`, which `_filter_carrier_loop` ignores, and
# the `atan` is never run.
@inline _driver_fll_discriminator(
    ::ThirdOrderAssistedBilinearLF,
    signal::AbstractGNSSSignal,
    folded,
) = fll_disc(
    signal,
    folded.filtered_correlator,
    folded.previous_prompt,
    folded.integration_time,
)
@inline _driver_fll_discriminator(::AbstractLoopFilter, ::AbstractGNSSSignal, _) = nothing

# Combining: the weighted mean, per loop, of the driver's own record and the
# coincident passenger records summed in `passengers`. A loop no passenger
# contributed to gets the driver's own value unchanged (see
# `_combined_discriminator`).
@inline function _loop_discriminators(
    passengers::DiscriminatorAccumulator,
    signal::AbstractGNSSSignal,
    ::AbstractLoopFilter,
    folded,
    integrated_samples::Integer,
    ctx::RecordFoldContext,
)
    own = _weighted_record_contribution(
        signal,
        folded.filtered_correlator,
        integrated_samples,
        folded.previous_prompt,
        nothing,
        0.0,
        ctx.code_doppler,
        ctx.sampling_frequency,
    )
    (;
        pll = _combined_discriminator(
            own.pll,
            own.pll_weight,
            passengers.pll_sum,
            passengers.pll_weight,
        ),
        fll = _combined_discriminator(
            own.fll,
            own.fll_weight,
            passengers.fll_sum,
            passengers.fll_weight,
        ),
        dll = _combined_discriminator(
            own.dll,
            own.dll_weight,
            passengers.dll_sum,
            passengers.dll_weight,
        ),
    )
end

# ---------------------------------------------------------------------------
# Signal processing entry point
# ---------------------------------------------------------------------------
#
# `_seed_and_fold` is the seam `_update_tracked_sat_doppler` calls. All it decides
# is whether passenger records may combine (`Val(true)`) or not (`Val(false)`),
# and what the estimator collects on the side (`_with_measurements`).
#
# The empty/non-empty passenger split is at the *type* level, so a single-signal
# satellite compiles to exactly the driver-only code — no combining arithmetic in
# the N=1 hot path, and bit-identical Dopplers whatever the flag says. The flag
# itself is a runtime branch, so that flipping it changes no container's
# concrete type; the traversal specialises on each side of it.

@inline _seed_and_fold(
    driver::TrackedSignal,
    passengers::Tuple{},
    sat::TrackedSat,
    state,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    ::Bool,
) = _with_measurements(state) do measurements
    _fold_satellite_chunk(
        driver,
        passengers,
        sat,
        state,
        sampling_frequency,
        noise,
        driver_carrier_phase_offset,
        Val(false),
        measurements,
    )
end

@inline _seed_and_fold(
    driver::TrackedSignal,
    passengers::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    sat::TrackedSat,
    state,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    discriminator_combining::Bool,
) = _with_measurements(state) do measurements
    fold(combine) = _fold_satellite_chunk(
        driver,
        passengers,
        sat,
        state,
        sampling_frequency,
        noise,
        driver_carrier_phase_offset,
        combine,
        measurements,
    )
    discriminator_combining ? fold(Val(true)) : fold(Val(false))
end

# Walk the driver's records in order and, before closing the loops on each one,
# advance every passenger up to it: passenger records that ended earlier are
# applied to their signal only, and a passenger record that ended on the very
# same sample (`sample_index`) is applied *and* combined into that driver
# record's loop update. After the last driver record the passengers are drained
# — applied, never combined. Nothing is carried to the next chunk: a passenger
# record that coincides with no driver record simply does not reach the loops.
#
# This is the only chunk traversal in the package. A single-signal satellite
# reaches it with an empty passenger tuple and `Val(false)`, which erases the
# passenger walk and all the combining arithmetic at compile time, leaving the
# plain driver fold it has always run.
#
# Returns `(new_driver, new_passengers, new_state, carrier_doppler, code_doppler)`
# — the chunk's last Dopplers, since the NCO is written once per chunk. With no
# driver record the Dopplers hold.
@inline function _fold_satellite_chunk(
    driver::TrackedSignal,
    passengers::Tuple,
    sat::TrackedSat,
    state,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    combine::Val,
    measurements,
)
    signal = driver.signal
    driver_outputs = driver.correlator_outputs
    contexts = _passenger_fold_contexts(
        passengers,
        Base.tail(noise),
        sat.prn,
        sampling_frequency,
        sat.code_doppler,
        driver_carrier_phase_offset,
        get_group_delay(driver),
    )
    driver_context = _driver_fold_context(
        driver,
        sat.prn,
        sampling_frequency,
        sat.code_doppler,
        first(noise),
    )

    tsig = driver
    cursors = map(_ -> 1, passengers)
    carrier_loop_filter = state.carrier_loop_filter
    code_loop_filter = state.code_loop_filter
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler

    @inbounds for k in eachindex(driver_outputs)
        output = driver_outputs[k]
        passengers, cursors, coincident, measurements = _advance_passengers_to(
            combine,
            passengers,
            cursors,
            contexts,
            output.sample_index,
            measurements,
        )

        folded = _advance_driver_record(
            tsig,
            output,
            driver_context,
            state.carrier_loop_filter_bandwidth,
            state.code_loop_filter_bandwidth,
        )
        tsig = folded.tracked_signal
        measurements = _collect_measurements(
            measurements,
            signal,
            output,
            folded.filtered_correlator,
            folded.previous_prompt,
            driver_context,
        )

        discriminators = _loop_discriminators(
            coincident,
            signal,
            carrier_loop_filter,
            folded,
            output.integrated_samples,
            driver_context,
        )
        closed = _close_loops(
            state,
            carrier_loop_filter,
            code_loop_filter,
            discriminators,
            folded.integration_time,
            folded.carrier_bandwidth,
            folded.code_bandwidth,
        )
        carrier_loop_filter = closed.carrier_loop_filter
        code_loop_filter = closed.code_loop_filter
        carrier_doppler, code_doppler = aid_dopplers(
            signal,
            state.init_carrier_doppler,
            state.init_code_doppler,
            closed.carrier_freq_update,
            closed.code_freq_update,
        )
    end
    empty!(driver_outputs)

    passengers, measurements =
        _drain_passengers(passengers, cursors, contexts, measurements)

    (
        tsig,
        passengers,
        _park_state(state, carrier_loop_filter, code_loop_filter, measurements),
        carrier_doppler,
        code_doppler,
    )
end

"""
$(SIGNATURES)

Estimate Dopplers and filter prompts for all satellites where the correlation has reached
the end of the code or multiples of that. This function uses the
conventional PLL and DLL implementation to estimate Dopplers for
carrier and code. Those Doppler estimations will be used to create the next
replicas to downconvert and decode the incoming signal. In addition to the
Doppler estimation it will also filter the prompt with the configured
post correlation filter.
In the case that the that the correlation hasn't reached the end, e.g. in the case
the incoming signal did not provide enough samples, it will return struct with
zeroed values.

The sampling-frequency argument may be either a [`BandMeasurements`](@ref)
NamedTuple (the `track` path, from which the per-band rate is read) or a bare
per-band sampling-frequency source — a `NamedTuple`/`Dict` keyed by
`get_band_id` mapping each band to its sampling frequency. The latter is
the entry point for an **external correlator producer** (e.g. an FPGA): it needs
no sample buffer, only the per-signal `correlator_outputs` and the rate that
maps `integrated_samples` to an integration time. See
[External correlator producers](@ref).
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Detach the slot *values* from the input (sharing the key set), then
    # delegate to the in-place form. This step never changes the key set, so
    # sharing the `Indices` is safe and avoids copying the hash table every
    # `track` loop iteration; the key set is detached once at the `track`
    # boundary (`reset_start_sample_and_bit_buffer`, #123). The per-sat
    # doppler update is identical between the two forms — only the storage
    # ownership differs.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, sampling_frequencies)
end

# Per-band sampling frequency for a group, from either a `BandMeasurements`
# NamedTuple (read the rate off the band's `BandMeasurement`) or a bare
# per-band rate source keyed by `get_band_id` (a `NamedTuple`/`Dict`). Both are
# looked up by the group's band id, so the estimator stays per-band and is
# never handed a scalar (groups may sit on different bands).
@inline _band_sampling_frequency(m::BandMeasurements, key) = m[key].sampling_frequency
@inline _band_sampling_frequency(fs::NamedTuple, key) = fs[key]
@inline _band_sampling_frequency(fs::AbstractDict, key) = fs[key]

# Per-signal noise density, as `(density, ready)` — the union-free pair the fold
# threads down to `_apply_correlator_output`. Looked up by signal id, because the
# floor a record divides by is the *post-correlation* one and that is a property
# of the despreading modulation, not of the RF band (see
# [`AbstractNoiseEstimator`](@ref)).
#
# Three cases, and only the first two are ever seen by a shipped estimator:
#
#   * no entry for the signal — no `AbstractNoiseEstimator` is configured, a
#     *static* property of the setup. The density is `nothing`, which makes the
#     context's `N` type parameter `Nothing` and `NoiseRefCN0Estimator.update`
#     throw. `ready` is `true`, because there is nothing to wait for: the wiring
#     mistake must surface at the first record, not be silently skipped forever.
#   * an entry whose window is still empty — a *runtime* condition. The density
#     is a dummy scalar of the right type and `ready` is `false`, so the fold
#     skips the update for the estimators that would read it, without ever
#     letting a `Union{Nothing,D}` into the context.
#   * an entry with a density — the normal case.
@inline _signal_noise_density(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? _noise_density_and_ready(noise_estimators[K]) :
    (nothing, true)

# One `(density, ready)` pair per signal of a group's slot type, in the slot
# type's own order — which is the order of every satellite's `signals` tuple, so
# the fold can pair them off positionally with `first`/`Base.tail` and never
# needs a lookup per satellite. The whole tuple folds to a constant shape: the
# signal ids come out of `TrackedSignal`'s type parameters and `haskey` on a
# NamedTuple is compile-time.
@inline _signal_noise_densities(
    noise_estimators::NamedTuple,
    ::Type{<:TrackedSat{Signals}},
) where {Signals} = _signal_noise_densities(noise_estimators, Signals)
@inline _signal_noise_densities(::NamedTuple, ::Type{Tuple{}}) = ()
@inline _signal_noise_densities(noise_estimators::NamedTuple, ::Type{T}) where {T<:Tuple} =
    (
        _signal_noise_density(
            noise_estimators,
            Val(get_signal_id(_signal_type(Base.tuple_type_head(T)))),
        ),
        _signal_noise_densities(noise_estimators, Base.tuple_type_tail(T))...,
    )

# Per-group body for the doppler estimator. Pulled out so
# `_foreach_group!` can call it without boxing when the groups tuple
# is heterogeneous (e.g. GPS L1 + Galileo E1B). The per-signal type
# is recovered from each signal inside `_update_tracked_sat_doppler`, and the
# group's `discriminator_combining` flag is handed to every satellite of it.
# Routes to this group's band's sampling frequency (see
# `_band_sampling_frequency`), and to each signal's own noise density (see
# `_signal_noise_densities`) — the sampling rate is a band property, the noise
# floor is not.
@inline function _est_one_group!(
    g::SignalGroup,
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
    noise_estimators::NamedTuple,
)
    vals = g.satellites.values
    isempty(vals) && return nothing
    combining = g.discriminator_combining
    sampling_frequency = _band_sampling_frequency(sampling_frequencies, get_band_id(g.band))
    noise = _signal_noise_densities(noise_estimators, eltype(g.satellites))
    _warn_noise_density_missing(eltype(g.satellites), noise, noise_estimators)
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_doppler(vals[i], sampling_frequency, noise, combining)
    end
    return nothing
end

# A signal with a configured noise estimator that has no *usable* density is a
# loud *symptom* — every satellite reports `-Inf dB-Hz` on it — but not a loud
# diagnosis. There are two causes, and the message names both because the fold
# cannot tell them apart from here:
#
#   - The window is still empty. The likeliest hardware-integration mistake: a
#     `CorrelatorNoiseEstimator` configured per the docs whose
#     `append_noise_observation!` is never called, so the static "no source
#     configured" check does not fire and the runtime skip repeats forever.
#   - The window holds a floor of zero, i.e. the input carries no power at all —
#     a front-end dropout or a buffer underrun. `_noise_density_and_ready`
#     reports that as not-ready rather than dividing by it (see there), so it
#     arrives here as the same flag. This is the one cause that *can* reach a
#     software-path caller, and the one that can appear after a signal has been
#     reporting a real C/N₀ for a while.
#
# This is the only point that knows a fold actually ran *and* the density was
# unusable, so the warning belongs here. `maxlog` is
# keyed per callsite rather than per signal, so the `_id` is made
# signal-specific — otherwise a second misconfigured signal would be silenced by
# the first, and the message names the signal precisely because a multi-signal
# setup is where the mistake is most likely.
#
# Walks the slot type and the density tuple in lockstep: `requires_noise_density`
# is a compile-time constant per signal, so all that survives is one branch on
# each signal's runtime `ready` flag.
#
# Warn rather than throw, and the asymmetry with the static case is the point:
# "no source configured" is unambiguously a mistake, whereas "no usable density at
# this instant" has legitimate transient readings (a producer that folds before it
# appends, a buffer shorter than one sub-integration, or a momentary dead input),
# so making it fatal would break a caller streaming short buffers.
#
# Takes the estimators as well as the `(density, ready)` pairs, for one reason: a
# multi-antenna window that is merely *filling* to its own dimension count is
# not-ready too, and that is a normal startup transient rather than a
# misconfiguration. Reading `_noise_window_filling` off the estimator keeps the
# fold's pair two-state — nothing downstream has to learn a third case — while
# still telling the two apart here, which is the only place that cares.
@inline _warn_noise_density_missing(
    ::Type{<:TrackedSat{Signals}},
    noise::Tuple,
    noise_estimators::NamedTuple,
) where {Signals} = _warn_noise_density_missing(Signals, noise, noise_estimators)
@inline _warn_noise_density_missing(::Type{Tuple{}}, ::Tuple{}, ::NamedTuple) = nothing
@inline function _warn_noise_density_missing(
    ::Type{T},
    noise::Tuple,
    noise_estimators::NamedTuple,
) where {T<:Tuple}
    head = Base.tuple_type_head(T)
    if !last(first(noise)) && requires_noise_density(_cn0_estimator_type(head))
        signal_id = get_signal_id(_signal_type(head))
        _noise_window_still_filling(noise_estimators, Val(signal_id)) ||
            _emit_noise_density_warning(signal_id)
    end
    _warn_noise_density_missing(Base.tuple_type_tail(T), Base.tail(noise), noise_estimators)
end

# `haskey` on a NamedTuple is compile-time, so the no-estimator case (a static
# wiring mistake, which must warn) folds away rather than being tested per fold.
@inline _noise_window_still_filling(noise_estimators::NamedTuple, ::Val{K}) where {K} =
    haskey(noise_estimators, K) ? _noise_window_filling(noise_estimators[K]) : false

@noinline function _emit_noise_density_warning(signal_id::Symbol)
    @warn """
          Signal `:$signal_id` has a noise estimator but no usable noise density, so \
          `NoiseRefCN0Estimator` will not update and C/N₀ stays at `-Inf dB-Hz`. \
          Either its window is still empty — call `append_noise_observation!` before \
          `estimate_dopplers_and_filter_prompt!`, or use a noise estimator that \
          measures from the band's samples — or the measured floor is zero, which \
          means the samples carry no power at all (a dead input or an underrun).""" _id =
        Symbol(:no_noise_density_, signal_id) maxlog = 1
    nothing
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref). Walks each
group's `Vector{TrackedSat}` backing storage and overwrites slots with the
new immutable `TrackedSat` value. Returns the same `track_state` object —
allocation-free in steady state when [`track!`](@ref)'s preconditions are met.

As with the immutable form, the second argument is either a
[`BandMeasurements`](@ref) NamedTuple or a bare per-band sampling-frequency
source keyed by `get_band_id`. The estimator only reads the per-band
sampling frequency (to turn each output's `integrated_samples` into an
integration time and to normalize the DLL discriminator); it consumes each
signal's `correlator_outputs` and **clears** them, whether they were produced
by the software correlate phase or appended by an external producer via
[`append_correlator_output!`](@ref).
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:SignalGroups,<:ConventionalPLLAndDLL},
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    _foreach_group!(
        _est_one_group!,
        track_state.groups,
        sampling_frequencies,
        track_state.noise_estimators,
    )
    return track_state
end

# Per-record loop updates on one correlator's own discriminators. The shared fold
# no longer calls these — it forms the discriminators first so they can be
# combined across signals (see `_loop_discriminators` / `_close_loops`) — but
# they are kept for callers that close a single signal's loops by hand.
function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    fll_discriminator = fll_disc(signal, correlator, previous_prompt, integration_time)
    filter_loop(
        carrier_loop_filter,
        (pll_discriminator, fll_discriminator),
        integration_time,
        loop_bandwidth,
    )
end

function calculate_carrier_frequency_update(
    signal::AbstractGNSSSignal,
    carrier_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    loop_bandwidth,
)
    pll_discriminator = pll_disc(signal, correlator)
    filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)
end

function calculate_code_frequency_update(
    signal::AbstractGNSSSignal,
    code_loop_filter::AbstractLoopFilter,
    correlator::AbstractCorrelator,
    code_doppler,
    sampling_frequency,
    integration_time,
    loop_bandwidth,
)
    dll_discriminator = dll_disc(signal, correlator, code_doppler, sampling_frequency)
    filter_loop(code_loop_filter, dll_discriminator, integration_time, loop_bandwidth)
end
