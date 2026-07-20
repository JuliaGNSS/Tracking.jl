"""
Per-satellite state for the vector PLL and DLL Doppler estimator
([`VectorPLLAndDLL`](@ref)).

On top of the conventional loop-filter state it carries the
vector-tracking (VT) interface to an external navigation filter
(e.g. GNSSReceiver.jl's VDFLL):

  - `code_discr_acc` / `carrier_discr_acc`: `(count, sum)` accumulators of
    the DLL discriminator (chips) and the FLL discriminator (Hz) since the
    navigation filter last read and reset them
    ([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
    Only accumulated while `vt_on`.
  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.
  - `vt_on`: whether the navigation filter controls this satellite's NCOs.
    While `false` the satellite runs a conventional (scalar) PLL/DLL as a
    fallback and nothing is accumulated. Promoted by [`update_vt_states!`](@ref)
    / [`set_vt_on!`](@ref).
"""
@kwdef struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    code_discr_acc::Tuple{Int,Float64} = (0, 0.0)
    code_freq_update::typeof(0.0Hz) = 0.0Hz
    carrier_discr_acc::Tuple{Int,typeof(0.0Hz)} = (0, 0.0Hz)
    carrier_freq_update::typeof(0.0Hz) = 0.0Hz
    vt_on::Bool = false
end

function SatVectorPLLAndDLL(
    sat::TrackedSat,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

function SatVectorPLLAndDLL(
    sat_vector_pll_and_dll::SatVectorPLLAndDLL{CA,CO};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_discr_acc::Maybe{Tuple{Int,Float64}} = nothing,
    code_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    carrier_discr_acc::Maybe{Tuple{Int,typeof(0.0Hz)}} = nothing,
    carrier_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    vt_on::Maybe{Bool} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatVectorPLLAndDLL{CA,CO}(
        sat_vector_pll_and_dll.init_carrier_doppler,
        sat_vector_pll_and_dll.init_code_doppler,
        isnothing(carrier_loop_filter) ? sat_vector_pll_and_dll.carrier_loop_filter :
        carrier_loop_filter,
        isnothing(code_loop_filter) ? sat_vector_pll_and_dll.code_loop_filter :
        code_loop_filter,
        isnothing(carrier_loop_filter_bandwidth) ?
        sat_vector_pll_and_dll.carrier_loop_filter_bandwidth :
        carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ?
        sat_vector_pll_and_dll.code_loop_filter_bandwidth : code_loop_filter_bandwidth,
        isnothing(code_discr_acc) ? sat_vector_pll_and_dll.code_discr_acc : code_discr_acc,
        isnothing(code_freq_update) ? sat_vector_pll_and_dll.code_freq_update :
        code_freq_update,
        isnothing(carrier_discr_acc) ? sat_vector_pll_and_dll.carrier_discr_acc :
        carrier_discr_acc,
        isnothing(carrier_freq_update) ? sat_vector_pll_and_dll.carrier_freq_update :
        carrier_freq_update,
        isnothing(vt_on) ? sat_vector_pll_and_dll.vt_on : vt_on,
    )
end

"""
$(SIGNATURES)

Vector-tracking Phase-Locked Loop (PLL) and Delay-Locked Loop (DLL) Doppler
estimator. Configuration-only — per-satellite state lives in each
[`TrackedSat`](@ref) wrapper as a [`SatVectorPLLAndDLL`](@ref), produced via
[`init_estimator_state`](@ref).

In vector tracking, the per-satellite tracking loops are closed centrally by
a navigation filter (living outside this package, e.g. in GNSSReceiver.jl)
instead of by per-satellite loop filters. The division of labor per
integration:

  - This estimator accumulates each satellite's DLL / FLL discriminator
    outputs for the navigation filter to consume (and reset via
    [`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
  - The navigation filter feeds NCO corrections back via
    [`set_code_freq_updates!`](@ref) / [`set_carrier_freq_updates!`](@ref).
    While a satellite's `vt_on` flag is set, its code Doppler follows the
    navigation filter's `code_freq_update` directly and the FLL branch of
    the FLL-assisted carrier loop filter is driven by the navigation
    filter's `carrier_freq_update` (a vector delay / frequency lock loop)
    while the PLL branch still runs on the satellite's own discriminator.
  - Satellites with `vt_on` unset (fresh from acquisition, before
    [`update_vt_states!`](@ref) promotes them) run a conventional scalar
    PLL/DLL as a fallback.

Type parameters `CA` and `CO` select the carrier and code loop filter types.
The FLL-assisted `ThirdOrderAssistedBilinearLF` carrier filter is the
default — with any non-assisted filter the navigation filter's
`carrier_freq_update` has no input path into the carrier loop.

Each bandwidth field is `Maybe{typeof(1.0Hz)}`: a `nothing` field (the
default) means **auto** — [`init_estimator_state`](@ref) sizes the bandwidth
per satellite from that sat's estimator-driver signal (`signals[1]`) via
[`default_carrier_loop_filter_bandwidth`](@ref) /
[`default_code_loop_filter_bandwidth`](@ref) — the same sizing as the
conventional estimator, which the scalar fallback loop is. Like the
conventional estimator, the effective bandwidth is scaled by `1/N` at filter
time when a signal coherently integrates `N` primary code blocks.
"""
struct VectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
end

function VectorPLLAndDLL(
    ::Type{CA} = ThirdOrderAssistedBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(carrier_loop_filter_bandwidth, code_loop_filter_bandwidth)
end

# Kwarg-update constructor for tweaking bandwidths in place.
function VectorPLLAndDLL(
    pll_and_dll::VectorPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
    )
end

"""
$(SIGNATURES)

Build the per-satellite estimator state stored in a [`TrackedSat`](@ref) for a
satellite tracked under [`VectorPLLAndDLL`](@ref). Auto bandwidths (`nothing`
on the estimator) are resolved here, per satellite, from the sat's
estimator-driver signal (`signals[1]`).

New satellites start with `vt_on = false` — the scalar fallback loop —
until [`update_vt_states!`](@ref) (or [`set_vt_on!`](@ref)) promotes them
into the vector loop.
"""
function init_estimator_state(
    estimator::VectorPLLAndDLL{CA,CO},
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
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Per-sat bandwidth overrides and the
# `vt_on` flag survive the reset.
function _reset_estimator_state(
    ::VectorPLLAndDLL,
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatVectorPLLAndDLL},
)
    state = sat.doppler_estimator_state
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter = constructorof(typeof(state.carrier_loop_filter))(),
        code_loop_filter = constructorof(typeof(state.code_loop_filter))(),
        carrier_loop_filter_bandwidth = state.carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth = state.code_loop_filter_bandwidth,
        vt_on = state.vt_on,
    )
end

# Carrier loop filtering with an explicit FLL-branch input: the raw FLL
# discriminator in the scalar fallback, or the navigation filter's
# `carrier_freq_update` under vector closure. Only the FLL-assisted filter
# has an input path for it; any other filter runs on the PLL discriminator
# alone (and the vector carrier closure degrades to PLL-only).
function _filter_vector_carrier_loop(
    carrier_loop_filter::ThirdOrderAssistedBilinearLF,
    pll_discriminator,
    fll_input,
    integration_time,
    loop_bandwidth,
)
    filter_loop(
        carrier_loop_filter,
        (pll_discriminator, fll_input),
        integration_time,
        loop_bandwidth,
    )
end

function _filter_vector_carrier_loop(
    carrier_loop_filter::AbstractLoopFilter,
    pll_discriminator,
    fll_input,
    integration_time,
    loop_bandwidth,
)
    filter_loop(carrier_loop_filter, pll_discriminator, integration_time, loop_bandwidth)
end

# Process the estimator-driver signal (signals[1]) under vector tracking.
# Same structure as the conventional method (dispatching on the per-sat
# state type plugs this into the shared `_update_tracked_sat_doppler`), but
# with the vector closure: while `vt_on`, the navigation filter's NCO
# corrections replace the scalar DLL output and the FLL branch input, and
# the discriminator outputs are accumulated for the navigation filter.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    pll_and_dll_state::SatVectorPLLAndDLL,
    sampling_frequency,
    driver_carrier_phase::Real = 0.0,
)
    if !tracked_signal.is_integration_completed || tracked_signal.integrated_samples == 0
        return tracked_signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    # The previous fully-integrated prompt must be read off the *old*
    # signal — the shared advance overwrites it with this integration's.
    previous_prompt = get_last_fully_integrated_filtered_prompt(tracked_signal)
    integration_time = tracked_signal.integrated_samples / sampling_frequency
    new_signal, filtered_correlator, integrated_code_blocks =
        _advance_signal_after_integration(
            tracked_signal,
            sat.prn,
            sampling_frequency,
            driver_carrier_phase,
        )

    # Same 1/N effective-bandwidth scaling as the conventional estimator —
    # holds the loop's BL·Δt stability product at its single-period value
    # when a signal coherently integrates N primary code blocks.
    carrier_bandwidth =
        pll_and_dll_state.carrier_loop_filter_bandwidth / integrated_code_blocks
    code_bandwidth = pll_and_dll_state.code_loop_filter_bandwidth / integrated_code_blocks

    pll_discriminator = pll_disc(signal, filtered_correlator)
    fll_discriminator =
        fll_disc(signal, filtered_correlator, previous_prompt, integration_time)
    dll_discriminator =
        dll_disc(signal, filtered_correlator, sat.code_doppler, sampling_frequency)

    if pll_and_dll_state.vt_on
        # Vector loop closure: accumulate the discriminator outputs for the
        # navigation filter and apply its NCO corrections — `code_freq_update`
        # directly, `carrier_freq_update` through the FLL branch of the
        # carrier loop filter (the PLL branch still runs on this satellite's
        # own discriminator). The code loop filter is bypassed and keeps its
        # state.
        code_discr_acc = pll_and_dll_state.code_discr_acc .+ (1, dll_discriminator)
        carrier_discr_acc = pll_and_dll_state.carrier_discr_acc .+ (1, fll_discriminator)
        code_freq_update = pll_and_dll_state.code_freq_update
        code_loop_filter = pll_and_dll_state.code_loop_filter
        fll_input = pll_and_dll_state.carrier_freq_update
    else
        # Scalar tracking fallback — a conventional FLL-assisted PLL and DLL.
        code_discr_acc = pll_and_dll_state.code_discr_acc
        carrier_discr_acc = pll_and_dll_state.carrier_discr_acc
        code_freq_update, code_loop_filter = filter_loop(
            pll_and_dll_state.code_loop_filter,
            dll_discriminator,
            integration_time,
            code_bandwidth,
        )
        fll_input = fll_discriminator
    end
    carrier_freq_update, carrier_loop_filter = _filter_vector_carrier_loop(
        pll_and_dll_state.carrier_loop_filter,
        pll_discriminator,
        fll_input,
        integration_time,
        carrier_bandwidth,
    )
    carrier_doppler, code_doppler = aid_dopplers(
        signal,
        pll_and_dll_state.init_carrier_doppler,
        pll_and_dll_state.init_code_doppler,
        carrier_freq_update,
        code_freq_update,
    )
    # `carrier_freq_update` is deliberately NOT stored back into the state:
    # the stored field carries the navigation filter's correction, which the
    # carrier-loop-filter output would otherwise overwrite between two
    # `set_carrier_freq_updates!` calls.
    new_doppler_estimator_state = SatVectorPLLAndDLL(
        pll_and_dll_state;
        carrier_loop_filter,
        code_loop_filter,
        code_discr_acc,
        code_freq_update,
        carrier_discr_acc,
    )
    return new_signal, new_doppler_estimator_state, carrier_doppler, code_doppler
end

"""
$(SIGNATURES)

Estimate Dopplers and filter prompts for all satellites where the correlation
has reached the end of the code or multiples of that, using the vector
PLL and DLL implementation — see [`VectorPLLAndDLL`](@ref) for how the
per-satellite loops are closed. In the case that the correlation hasn't
reached the end, e.g. in the case the incoming signal did not provide enough
samples, the state is passed through unchanged.
"""
function estimate_dopplers_and_filter_prompt(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    measurements::BandMeasurements,
)
    # Detach the slot *values* (sharing the key set) then delegate to the
    # in-place form — same key-sharing copy the conventional estimator uses;
    # see its `estimate_dopplers_and_filter_prompt` for the rationale.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, measurements)
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref) for the
vector PLL and DLL estimator.
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    measurements::BandMeasurements,
)
    _foreach_group!(_est_one_group!, track_state.groups, measurements)
    return track_state
end

# Shared in-place walker for the vector-tracking state managers below:
# overwrite each satellite's `doppler_estimator_state` with the per-sat rule
# `f(sat, sat.doppler_estimator_state, args...)`, threading each manager's own
# arguments (`prns`, a freq-update table, …) through unchanged — the same
# named-worker-plus-trailing-args pattern `_foreach_group!` uses. `f` must
# return the same concrete `SatVectorPLLAndDLL` type (guaranteed by the
# kwarg-update constructor).
#
# Two entry points give the managers their two addressing forms:
# `_map_vt_states!` walks every group; `_map_vt_states_in_group!` walks only the
# addressed group (Symbol / Integer / Val). Keeping them as distinct names —
# rather than one function overloaded on the group position — is what lets a
# threaded `Symbol`/`Integer` argument not be mistaken for a group selector;
# `Bool <: Integer` makes that a real hazard for `set_vt_on!`. The group-scoped
# form is what makes the managers usable in a multi-constellation receiver,
# where PRNs alone are ambiguous (GPS PRN 5 and Galileo PRN 5 are different
# satellites in different groups).
@inline function _map_vt_states!(
    f::F,
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    args::Vararg{Any,N},
) where {F,N}
    _foreach_group!(_map_vt_states_one_group!, track_state.groups, f, args...)
    return track_state
end

@inline function _map_vt_states_in_group!(
    f::F,
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    args::Vararg{Any,N},
) where {F,N}
    _map_vt_states_one_group!(_index_group(track_state.groups, group), f, args...)
    return track_state
end

@inline function _map_vt_states_one_group!(
    g::SignalGroup,
    f::F,
    args::Vararg{Any,N},
) where {F,N}
    vals = g.satellites.values
    @inbounds for i in eachindex(vals)
        sat = vals[i]
        vals[i] = TrackedSat(
            sat;
            doppler_estimator_state = f(sat, sat.doppler_estimator_state, args...),
        )
    end
    return nothing
end

# The per-satellite promotion rule of `update_vt_states!`. Promotion is
# one-directional: a satellite whose PRN is in lock and is not yet in the
# vector loop joins it; nothing is ever demoted here (use `set_vt_on!` to drop
# a satellite).
_update_vt_membership(sat, state, prns_in_lock) =
    !state.vt_on && sat.prn in prns_in_lock ? SatVectorPLLAndDLL(state; vt_on = true) :
    state

"""
$(SIGNATURES)

Promote every satellite whose PRN is in `prns_in_lock` (the set the receiver
currently considers usable, e.g. from a carrier-to-noise-density threshold)
into the vector loop by setting its `vt_on` flag; a satellite keeps running
the scalar fallback loop until it appears here. Promotion is one-directional —
a satellite already in the vector loop is left untouched, and none is ever
demoted (use [`set_vt_on!`](@ref) to drop one).

The two-argument form walks every group; pass a `group` (Symbol or index) to
address a single group — required in a multi-constellation receiver, where a
PRN alone is ambiguous across groups.

Mutates `track_state` in place and returns it.
"""
function update_vt_states!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    prns_in_lock,
)
    _map_vt_states!(_update_vt_membership, track_state, prns_in_lock)
end

function update_vt_states!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    prns_in_lock,
)
    _map_vt_states_in_group!(_update_vt_membership, track_state, group, prns_in_lock)
end

"""
$(SIGNATURES)

Set the `vt_on` flag to `vt_on` for every satellite whose PRN is in `prns` —
the manual override for [`update_vt_states!`](@ref)'s automatic promotion.
Walks every group, or only `group` when one is given (the unambiguous form
for multi-constellation receivers). Mutates `track_state` in place and
returns it.
"""
function set_vt_on!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    vt_on::Bool,
    prns,
)
    _map_vt_states!(_set_sat_vt_on, track_state, vt_on, prns)
end

function set_vt_on!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    vt_on::Bool,
    prns,
)
    _map_vt_states_in_group!(_set_sat_vt_on, track_state, group, vt_on, prns)
end

_set_sat_vt_on(sat, state, vt_on, prns) =
    sat.prn in prns ? SatVectorPLLAndDLL(state; vt_on) : state

"""
$(SIGNATURES)

Reset the code (DLL) discriminator accumulator of every satellite in the
vector loop — called by the navigation filter after it has consumed the
accumulated values. Mutates `track_state` in place and returns it.
"""
function reset_code_discr_acc!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL})
    _map_vt_states!(_reset_sat_code_discr_acc, track_state)
end

_reset_sat_code_discr_acc(_sat, state) =
    state.vt_on ? SatVectorPLLAndDLL(state; code_discr_acc = (0, 0.0)) : state

"""
$(SIGNATURES)

Reset the carrier (FLL) discriminator accumulator of every satellite in the
vector loop — called by the navigation filter after it has consumed the
accumulated values. Mutates `track_state` in place and returns it.
"""
function reset_carrier_discr_acc!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL})
    _map_vt_states!(_reset_sat_carrier_discr_acc, track_state)
end

_reset_sat_carrier_discr_acc(_sat, state) =
    state.vt_on ? SatVectorPLLAndDLL(state; carrier_discr_acc = (0, 0.0Hz)) : state

"""
$(SIGNATURES)

Mean DLL (code) discriminator accumulated on `state` since the last
[`reset_code_discr_acc!`](@ref), in chips, or `nothing` if nothing has been
accumulated yet (the `(count, sum)` accumulator's `count` is 0). This is the
single place the accumulator's averaging convention lives — read it here
rather than dividing `code_discr_acc` by hand, so a change to how the
accumulator is stored can't silently diverge between consumers.
"""
function mean_code_discriminator(state::SatVectorPLLAndDLL)
    count, discr_sum = state.code_discr_acc
    count == 0 ? nothing : discr_sum / count
end

"""
$(SIGNATURES)

Mean FLL (carrier) discriminator accumulated on `state` since the last
[`reset_carrier_discr_acc!`](@ref), in Hz, or `nothing` if nothing has been
accumulated yet (`count == 0`). The carrier counterpart to
[`mean_code_discriminator`](@ref).
"""
function mean_carrier_discriminator(state::SatVectorPLLAndDLL)
    count, discr_sum = state.carrier_discr_acc
    count == 0 ? nothing : discr_sum / count
end

"""
$(SIGNATURES)

Feed the navigation filter's code-frequency NCO corrections back into the
tracking loops. `code_freq_updates` maps PRN to the correction (anything
indexable by PRN, e.g. a `Dictionary`) and must have an entry for every
satellite in the vector loop; satellites with `vt_on` unset are skipped.
Walks every group, or only `group` when one is given — the group-scoped form
is required in a multi-constellation receiver, where each group carries its
own corrections and PRNs collide across groups. Mutates `track_state` in
place and returns it.
"""
function set_code_freq_updates!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    code_freq_updates,
)
    _map_vt_states!(_set_sat_code_freq_update, track_state, code_freq_updates)
end

function set_code_freq_updates!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    code_freq_updates,
)
    _map_vt_states_in_group!(
        _set_sat_code_freq_update,
        track_state,
        group,
        code_freq_updates,
    )
end

_set_sat_code_freq_update(sat, state, code_freq_updates) =
    state.vt_on ? SatVectorPLLAndDLL(state; code_freq_update = code_freq_updates[sat.prn]) :
    state

"""
$(SIGNATURES)

Feed the navigation filter's carrier-frequency NCO corrections back into the
tracking loops. `carrier_freq_updates` maps PRN to the correction (anything
indexable by PRN, e.g. a `Dictionary`) and must have an entry for every
satellite in the vector loop; satellites with `vt_on` unset are skipped.
Walks every group, or only `group` when one is given — the group-scoped form
is required in a multi-constellation receiver, where each group carries its
own corrections and PRNs collide across groups. Mutates `track_state` in
place and returns it.
"""
function set_carrier_freq_updates!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    carrier_freq_updates,
)
    _map_vt_states!(_set_sat_carrier_freq_update, track_state, carrier_freq_updates)
end

function set_carrier_freq_updates!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    carrier_freq_updates,
)
    _map_vt_states_in_group!(
        _set_sat_carrier_freq_update,
        track_state,
        group,
        carrier_freq_updates,
    )
end

_set_sat_carrier_freq_update(sat, state, carrier_freq_updates) =
    state.vt_on ?
    SatVectorPLLAndDLL(state; carrier_freq_update = carrier_freq_updates[sat.prn]) : state
