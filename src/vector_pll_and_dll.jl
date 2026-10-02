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
    **One pair per signal**, in `sat.signals` order (see
    [Multi-signal satellites](@ref)); read them with [`mean_code_discr`](@ref) /
    [`mean_carrier_discr`](@ref). Only accumulated while `vt_on`.
  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.
  - `vt_on`: whether the navigation filter controls this satellite's NCOs.
    While `false` the satellite runs a conventional (scalar) PLL/DLL as a
    fallback and nothing is accumulated. Set by [`enable_vt!`](@ref) /
    [`disable_vt!`](@ref).
  - `group_delays`: the per-signal group delays (`NaN` for unknown), as for
    the conventional estimator's state; see [`set_group_delay!`](@ref).

Type parameter `N` is the satellite's signal count and `C` whether the
passengers' discriminators are combined into this satellite's loops
(`combine_discriminators`, see [`VectorPLLAndDLL`](@ref)).

The per-signal accumulators have no default: their length is the satellite's
signal count. Build the state from the satellite (`SatVectorPLLAndDLL(sat, …)`),
through [`init_estimator_state`](@ref), or by keyword with both accumulators;
`group_delays` then defaults to every delay unknown and takes `nothing` for an
unknown entry.
"""
struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N,C}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA
    code_loop_filter::CO
    carrier_loop_filter_bandwidth::typeof(1.0Hz)
    code_loop_filter_bandwidth::typeof(1.0Hz)
    code_discr_acc::NTuple{N,Tuple{Int,Float64}}
    code_freq_update::typeof(0.0Hz)
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}}
    carrier_freq_update::typeof(0.0Hz)
    vt_on::Bool
    group_delays::NTuple{N,typeof(1.0s)}
end

function SatVectorPLLAndDLL(;
    init_carrier_doppler,
    init_code_doppler,
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF(),
    code_loop_filter::CO = SecondOrderBilinearLF(),
    carrier_loop_filter_bandwidth = 18.0Hz,
    code_loop_filter_bandwidth = 1.0Hz,
    code_discr_acc::NTuple{N,Tuple{Int,Float64}},
    code_freq_update = 0.0Hz,
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}},
    carrier_freq_update = 0.0Hz,
    vt_on::Bool = false,
    group_delays::NTuple{N,Any} = map(_ -> nothing, code_discr_acc),
    combine_discriminators::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    SatVectorPLLAndDLL{CA,CO,N,combine_discriminators}(
        init_carrier_doppler,
        init_code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        code_discr_acc,
        code_freq_update,
        carrier_discr_acc,
        carrier_freq_update,
        vt_on,
        map(_stored_group_delay, group_delays),
    )
end

# One zeroed accumulator per signal. Mapping over the signal tuple keeps the
# count a compile-time constant, so every satellite of a group has the same
# state type.
@inline _zero_code_discr_acc(sat::TrackedSat) = map(_ -> (0, 0.0), sat.signals)
@inline _zero_carrier_discr_acc(sat::TrackedSat) = map(_ -> (0, 0.0Hz), sat.signals)

# The same accumulators emptied, each `(count, sum)` zeroed in its own types.
@inline _zeroed(accs::Tuple) = map(acc -> zero.(acc), accs)

function SatVectorPLLAndDLL(
    sat::TrackedSat,
    carrier_loop_filter::AbstractLoopFilter,
    code_loop_filter::AbstractLoopFilter;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
    combine_discriminators::Bool = false,
)
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        code_discr_acc = _zero_code_discr_acc(sat),
        carrier_discr_acc = _zero_carrier_discr_acc(sat),
        combine_discriminators,
    )
end

function SatVectorPLLAndDLL(
    sat_vector_pll_and_dll::SatVectorPLLAndDLL{CA,CO,N,C};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_discr_acc::Maybe{NTuple{N,Tuple{Int,Float64}}} = nothing,
    code_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    carrier_discr_acc::Maybe{NTuple{N,Tuple{Int,typeof(0.0Hz)}}} = nothing,
    carrier_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    vt_on::Maybe{Bool} = nothing,
    group_delays::Maybe{NTuple{N,typeof(1.0s)}} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N,C}
    SatVectorPLLAndDLL{CA,CO,N,C}(
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
        isnothing(group_delays) ? sat_vector_pll_and_dll.group_delays : group_delays,
    )
end

_combines_discriminators(::SatVectorPLLAndDLL{<:Any,<:Any,<:Any,C}) where {C} = C
_group_delays(state::SatVectorPLLAndDLL) = state.group_delays
_with_group_delays(state::SatVectorPLLAndDLL, group_delays) =
    SatVectorPLLAndDLL(state; group_delays)
get_group_delay(state::SatVectorPLLAndDLL, signal_index::Integer) =
    _group_delay(state, signal_index)

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

  - This estimator accumulates every signal's DLL / FLL discriminator outputs
    for the navigation filter to consume (and reset via
    [`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
  - The navigation filter feeds NCO corrections back via
    [`set_code_freq_updates!`](@ref) / [`set_carrier_freq_updates!`](@ref).
    While a satellite's `vt_on` flag is set, its code Doppler follows the
    navigation filter's `code_freq_update` directly and the FLL branch of
    the FLL-assisted carrier loop filter is driven by the navigation
    filter's `carrier_freq_update` (a vector delay / frequency lock loop)
    while the PLL branch still runs on the satellite's own discriminator.
  - Satellites with `vt_on` unset (fresh from acquisition, before
    [`enable_vt!`](@ref) puts them in the loop) run a conventional scalar
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

`combine_discriminators = true` combines a multi-signal satellite's passengers
into the loops this package still closes, by the rules of
[`ConventionalPLLAndDLL`](@ref)'s `combine_discriminators` (see
[Discriminator combining](@ref)): all three while `vt_on` is unset, so that the
scalar fallback is [`ConventionalAssistedPLLAndDLL`](@ref) combining with the
same loop filters, and only the carrier phase loop under vector closure. Off by
default. The per-signal accumulators hold every signal's raw readings either
way. The flag is the type parameter `C`, as for the conventional estimator.
"""
struct VectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,C} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
end

function VectorPLLAndDLL(
    ::Type{CA} = ThirdOrderAssistedBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    combine_discriminators::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO,combine_discriminators}(
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Kwarg-update constructor for tweaking bandwidths in place.
function VectorPLLAndDLL(
    pll_and_dll::VectorPLLAndDLL{CA,CO,C};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,C}
    VectorPLLAndDLL{CA,CO,C}(
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
until [`enable_vt!`](@ref) puts them into the vector loop.
"""
function init_estimator_state(
    estimator::VectorPLLAndDLL{CA,CO,C},
    sat::TrackedSat,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,C}
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
    # Every signal starts with an unknown group delay.
    SatVectorPLLAndDLL{
        typeof(carrier_loop_filter),
        typeof(code_loop_filter),
        length(sat.signals),
        C,
    }(
        sat.carrier_doppler,
        sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        _zero_code_discr_acc(sat),
        0.0Hz,
        _zero_carrier_discr_acc(sat),
        0.0Hz,
        false,
        map(_ -> _UNKNOWN_GROUP_DELAY, sat.signals),
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Per-sat bandwidth overrides, the `vt_on`
# flag, the group delays and the combining flag survive the reset.
function _reset_estimator_state(
    ::VectorPLLAndDLL,
    sat::TrackedSat{<:Tuple{Vararg{TrackedSignal}},<:SatVectorPLLAndDLL},
)
    state = sat.doppler_estimator_state
    typeof(state)(
        sat.carrier_doppler,
        sat.code_doppler,
        constructorof(typeof(state.carrier_loop_filter))(),
        constructorof(typeof(state.code_loop_filter))(),
        state.carrier_loop_filter_bandwidth,
        state.code_loop_filter_bandwidth,
        _zeroed(state.code_discr_acc),
        0.0Hz,
        _zeroed(state.carrier_discr_acc),
        0.0Hz,
        state.vt_on,
        state.group_delays,
    )
end

# The vector closure of one record, plugged into the shared driver fold by
# dispatch on the per-sat state. While `vt_on`, the navigation filter's NCO
# corrections drive the loops — `code_freq_update` directly (code loop filter
# bypassed) and `carrier_freq_update` through the FLL branch of the carrier loop
# filter (the PLL branch still runs on this satellite's own discriminator) — and
# the raw DLL/FLL discriminator outputs are accumulated for the navigation
# filter. With `vt_on` unset it is the conventional scalar closure. Only the
# FLL-assisted carrier filter has an input path for `carrier_freq_update`; with
# any other filter the vector carrier closure degrades to PLL-only.
@inline function _close_loops(
    state::SatVectorPLLAndDLL,
    own,
    mixed,
    contributions,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
)
    code_loop_filter = state.code_loop_filter
    code_discr_acc = state.code_discr_acc
    carrier_discr_acc = state.carrier_discr_acc
    if state.vt_on
        # The driver's own measurements in its slot, and those of the passengers
        # combined in this step in theirs; the others follow in `_fold_driver`.
        code_discr_acc = _with_combined_measurements(
            Base.setindex(code_discr_acc, _accumulated(code_discr_acc[1], own.dll), 1),
            contributions,
            c -> c.dll,
        )
        carrier_discr_acc = _with_combined_measurements(
            Base.setindex(
                carrier_discr_acc,
                _accumulated(carrier_discr_acc[1], own.fll),
                1,
            ),
            contributions,
            c -> c.fll,
        )
        code_freq_update = state.code_freq_update
        fll_input = state.carrier_freq_update
    else
        code_freq_update, code_loop_filter = calculate_code_frequency_update(
            code_loop_filter,
            mixed.dll,
            integration_time,
            code_bandwidth,
        )
        fll_input = mixed.fll
    end
    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        state.carrier_loop_filter,
        mixed.pll,
        fll_input,
        integration_time,
        carrier_bandwidth,
    )
    # Neither NCO-correction field (`code_freq_update` / `carrier_freq_update`)
    # is written back: both are owned by the navigation filter (set via
    # `set_code_freq_updates!` / `set_carrier_freq_updates!`) and read as this
    # closure's inputs, so overwriting them with a loop output would clobber the
    # navigation filter's value between two of its calls.
    carrier_freq_update,
    code_freq_update,
    SatVectorPLLAndDLL(
        state;
        carrier_loop_filter,
        code_loop_filter,
        code_discr_acc,
        carrier_discr_acc,
    )
end

# One record's measurement added to a `(count, sum)` accumulator: the single
# accumulation rule, for the driver's slot and every passenger's.
@inline _accumulated(acc::Tuple, measurement) = acc .+ (1, measurement)

# The passengers' slots with the raw readings of the records combined in this
# driver step added (`reading` picks the DLL or FLL one); a passenger that did
# not take part, and a fold without combining (`nothing`), adds nothing.
@inline _with_combined_measurements(accs::Tuple, ::Nothing, _) = accs
@inline _with_combined_measurements(
    accs::Tuple,
    contributions::Tuple,
    reading::F,
) where {F} = (
    first(accs),
    map(Base.tail(accs), contributions) do acc, contribution
        contribution.coincident ? _accumulated(acc, reading(contribution)) : acc
    end...,
)

# The driver's step of the per-sat update under vector tracking: the shared
# driver fold (closing the loops through the vector closure above), with the
# coincident passengers combined when `combine_discriminators` is set, then
# every other passenger through the shared passenger fold, its raw readings added
# to its own accumulators while `vt_on`. The passengers come back with their
# records consumed, leaving nothing for the shared passenger fold.
@inline function _fold_driver(
    tracked_signal::TrackedSignal,
    passengers::Tuple,
    sat::TrackedSat,
    state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase::Real,
)
    # A runtime flag, but both calls return the same types.
    new_driver, driver_state, carrier_doppler, code_doppler, passengers =
        if _combines_discriminators(state) && !isempty(tracked_signal.correlator_outputs)
            _fold_driver_records(
                tracked_signal,
                passengers,
                sat,
                state,
                sampling_frequency,
                noise,
                driver_carrier_phase,
                _passenger_contexts(
                    tracked_signal,
                    passengers,
                    sat,
                    state,
                    noise,
                    driver_carrier_phase,
                ),
            )
        else
            _fold_driver_records(
                tracked_signal,
                passengers,
                sat,
                state,
                sampling_frequency,
                noise,
                driver_carrier_phase,
                nothing,
            )
        end
    vt_on = driver_state.vt_on
    measure(accs, signal, output, filtered_correlator, previous_prompt) =
        vt_on ?
        _accumulated_records(
            accs,
            _record_discriminators(
                signal.signal,
                filtered_correlator,
                previous_prompt,
                output.integrated_samples / sampling_frequency,
                sat.code_doppler,
                sampling_frequency,
            ),
        ) : accs
    folded = map(
        passengers,
        Base.tail(noise),
        Base.tail(driver_state.code_discr_acc),
        Base.tail(driver_state.carrier_discr_acc),
    ) do passenger, (noise_density, noise_density_ready), code_acc, carrier_acc
        _fold_passenger_records(
            passenger,
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase,
            (code_acc, carrier_acc),
            measure,
        )
    end
    new_state = SatVectorPLLAndDLL(
        driver_state;
        code_discr_acc = (
            first(driver_state.code_discr_acc),
            map(f -> last(f)[1], folded)...,
        ),
        carrier_discr_acc = (
            first(driver_state.carrier_discr_acc),
            map(f -> last(f)[2], folded)...,
        ),
    )
    new_driver, new_state, carrier_doppler, code_doppler, map(first, folded)
end

# A passenger's `(code, carrier)` accumulator pair with one record's readings added.
@inline _accumulated_records((code_acc, carrier_acc), discriminators) = (
    _accumulated(code_acc, discriminators.dll),
    _accumulated(carrier_acc, discriminators.fll),
)

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
    sampling_frequencies::Union{BandMeasurements,NamedTuple,AbstractDict},
)
    # Detach the slot *values* (sharing the key set) then delegate to the
    # in-place form — same key-sharing copy the conventional estimator uses;
    # see its `estimate_dopplers_and_filter_prompt` for the rationale.
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    estimate_dopplers_and_filter_prompt!(new_track_state, sampling_frequencies)
end

"""
$(SIGNATURES)

In-place version of [`estimate_dopplers_and_filter_prompt`](@ref) for the
vector PLL and DLL estimator. Like the conventional form, the second argument
is a [`BandMeasurements`](@ref) NamedTuple or a bare per-band
sampling-frequency source keyed by `get_band_id`; each signal's
`correlator_outputs` are consumed and cleared.
"""
function estimate_dopplers_and_filter_prompt!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
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
# threaded `Symbol`/`Integer` argument not be mistaken for a group selector.
# `Bool <: Integer` makes that a real hazard for the `vt_on` flag `enable_vt!` /
# `disable_vt!` thread through, which is precisely why that flag never appears
# in a public signature. The group-scoped form is what makes the managers usable
# in a multi-constellation receiver, where PRNs alone are ambiguous (GPS PRN 5
# and Galileo PRN 5 are different satellites in different groups).
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

# The per-satellite membership rule shared by `enable_vt!` / `disable_vt!`:
# write `vt_on` to the addressed PRNs, leave every other satellite as it is.
# The flag is a plain field write, so re-issuing the same membership every
# cycle is a no-op for satellites already in that state.
_set_sat_vt_on(sat, state, vt_on, prns) =
    sat.prn in prns ? SatVectorPLLAndDLL(state; vt_on) : state

"""
$(SIGNATURES)

Put every satellite whose PRN is in `prns` into the vector loop by setting its
`vt_on` flag: from the next integration on, the navigation filter's NCO
corrections drive its loops ([`set_code_freq_updates!`](@ref) /
[`set_carrier_freq_updates!`](@ref)) and its DLL/FLL discriminators are
accumulated for the filter to read. Satellites outside `prns` are left
untouched, and re-enabling one already in the loop changes nothing — pass the
currently usable set (e.g. the satellites in lock) every cycle rather than
only the newly eligible ones.

The two-argument form walks every group; pass a `group` (Symbol or index) to
address a single group — required in a multi-constellation receiver, where a
PRN alone is ambiguous across groups.

Mutates `track_state` in place and returns it. [`disable_vt!`](@ref) is the
inverse.
"""
function enable_vt!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL}, prns)
    _map_vt_states!(_set_sat_vt_on, track_state, true, prns)
end

function enable_vt!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    prns,
)
    _map_vt_states_in_group!(_set_sat_vt_on, track_state, group, true, prns)
end

"""
$(SIGNATURES)

Take every satellite whose PRN is in `prns` out of the vector loop by clearing
its `vt_on` flag — the inverse of [`enable_vt!`](@ref), for satellites the
navigation filter drops (unhealthy, diverged, no longer visible). Each falls
back to its own scalar PLL/DLL, whose loop filters still carry their
pre-vector state and whose Doppler now contains the navigation filter's last
NCO correction; follow up with [`reset_loop_filters!`](@ref) to re-seed the
scalar loop from the current Doppler and make the handoff transient-free.

Addressed exactly like [`enable_vt!`](@ref): every group, or only `group` when
one is given. Mutates `track_state` in place and returns it.
"""
function disable_vt!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL}, prns)
    _map_vt_states!(_set_sat_vt_on, track_state, false, prns)
end

function disable_vt!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
    prns,
)
    _map_vt_states_in_group!(_set_sat_vt_on, track_state, group, false, prns)
end

"""
$(SIGNATURES)

Reset the code (DLL) discriminator accumulator of every satellite in the
vector loop — called by the navigation filter after it has consumed the
accumulated values via [`mean_code_discr`](@ref). Satellites outside
the vector loop are left untouched (they accumulate nothing).

The one-argument form walks every group; pass a `group` (Symbol or index) to
reset a single group's satellites. Mutates `track_state` in place and returns
it.
"""
function reset_code_discr_acc!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL})
    _map_vt_states!(_reset_sat_code_discr_acc, track_state)
end

function reset_code_discr_acc!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
)
    _map_vt_states_in_group!(_reset_sat_code_discr_acc, track_state, group)
end

_reset_sat_code_discr_acc(_sat, state) =
    state.vt_on ?
    SatVectorPLLAndDLL(state; code_discr_acc = _zeroed(state.code_discr_acc)) : state

"""
$(SIGNATURES)

Reset the carrier (FLL) discriminator accumulator of every satellite in the
vector loop — called by the navigation filter after it has consumed the
accumulated values via [`mean_carrier_discr`](@ref). Addressed exactly
like [`reset_code_discr_acc!`](@ref): every group, or only `group` when
one is given. Mutates `track_state` in place and returns it.
"""
function reset_carrier_discr_acc!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL})
    _map_vt_states!(_reset_sat_carrier_discr_acc, track_state)
end

function reset_carrier_discr_acc!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
)
    _map_vt_states_in_group!(_reset_sat_carrier_discr_acc, track_state, group)
end

_reset_sat_carrier_discr_acc(_sat, state) =
    state.vt_on ?
    SatVectorPLLAndDLL(state; carrier_discr_acc = _zeroed(state.carrier_discr_acc)) : state

"""
$(SIGNATURES)

Mean DLL (code) discriminator of one signal, accumulated since the last
[`reset_code_discr_acc!`](@ref), in chips, or `nothing` if nothing has been
accumulated yet (the `(count, sum)` accumulator's `count` is 0). This is the
single place the accumulator's averaging convention lives — read it here
rather than dividing `code_discr_acc` by hand.

Each signal has its own accumulator; see [Multi-signal satellites](@ref) for what
it holds.

Addressed like the other per-signal accessors: from a [`TrackState`](@ref) or a
[`TrackedSat`](@ref) with a trailing signal selector — an index or a signal
type — which may be omitted only for a single-signal satellite. The
`SatVectorPLLAndDLL` form takes an index.

```julia
mean_code_discr(track_state, :galileo_e1, 11, GalileoE1B)  # by signal type
mean_code_discr(sat, 2)                                     # by index
mean_code_discr(get_doppler_estimator_state(sat))           # single-signal satellite
```
"""
function mean_code_discr(state::SatVectorPLLAndDLL, signal_index::Integer)
    count, discr_sum = state.code_discr_acc[signal_index]
    count == 0 ? nothing : discr_sum / count
end

"""
$(SIGNATURES)

Mean FLL (carrier) discriminator of one signal, accumulated since the last
[`reset_carrier_discr_acc!`](@ref), in Hz, or `nothing` if nothing has been
accumulated yet (`count == 0`). The carrier counterpart to
[`mean_code_discr`](@ref), addressed the same way.
"""
function mean_carrier_discr(state::SatVectorPLLAndDLL, signal_index::Integer)
    count, discr_sum = state.carrier_discr_acc[signal_index]
    count == 0 ? nothing : discr_sum / count
end

# A single-signal satellite needs no selector; a multi-signal one is refused, as
# an unqualified per-signal read is.
mean_code_discr(state::SatVectorPLLAndDLL) =
    mean_code_discr(state, _sole_signal_index(state))
mean_carrier_discr(state::SatVectorPLLAndDLL) =
    mean_carrier_discr(state, _sole_signal_index(state))
_sole_signal_index(::SatVectorPLLAndDLL{<:Any,<:Any,N}) where {N} =
    N == 1 ? 1 : _throw_needs_signal_selector()

# The satellite rung turns a signal type into a slot, since the estimator state
# holds the accumulators but not the signals; the `TrackState` rungs forward to
# it like the other per-signal accessors.
for fn in (:mean_code_discr, :mean_carrier_discr)
    @eval begin
        $fn(sat::TrackedSat, sel...) =
            $fn(get_doppler_estimator_state(sat), _signal_index(sat.signals, sel...))
        $fn(s::TrackState{<:SignalGroups,<:VectorPLLAndDLL}, id...) =
            $fn(get_sat_state(s, id...))
        $fn(
            s::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
            group::Union{Symbol,Integer,Val},
            sat_id,
            sig::_SignalSelector,
        ) = $fn(get_sat_state(s, group, sat_id), sig)
    end
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
