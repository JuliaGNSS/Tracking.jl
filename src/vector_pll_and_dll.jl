"""
Per-satellite state for the vector PLL and DLL Doppler estimator
([`VectorPLLAndDLL`](@ref)).

On top of the conventional loop-filter state it carries the
vector-tracking (VT) interface to an external navigation filter
(e.g. GNSSReceiver.jl's VDFLL):

  - `code_discr_accs` / `carrier_discr_accs`: `(count, sum)` accumulators of
    the DLL discriminator (chips) and the FLL discriminator (Hz) since the
    navigation filter last read and reset them
    ([`reset_code_discr_accs!`](@ref) / [`reset_carrier_discr_accs!`](@ref)),
    **one per signal**, in `sat.signals` order; read them with
    [`mean_code_discr`](@ref) / [`mean_carrier_discr`](@ref). Only
    accumulated while `vt_on`. They have no default, as their number is the
    satellite's signal count: build the state from the satellite
    (`SatVectorPLLAndDLL(sat, …)`) or through [`init_estimator_state`](@ref).
  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.
  - `vt_on`: whether the navigation filter controls this satellite's NCOs.
    While `false` the satellite runs a conventional (scalar) PLL/DLL as a
    fallback and nothing is accumulated. Set by [`enable_vt!`](@ref) /
    [`disable_vt!`](@ref).
  - `frequency_lock`: the carrier loop's [`FrequencyLockIndicator`](@ref), run
    as in the conventional estimator while `vt_on` is unset.
  - `signal_combining_sums`: the passengers' pending
    [`SignalCombiningSums`](@ref).
"""
@kwdef struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    code_discr_accs::NTuple{N,Tuple{Int,Float64}}
    code_freq_update::typeof(0.0Hz) = 0.0Hz
    carrier_discr_accs::NTuple{N,Tuple{Int,typeof(0.0Hz)}}
    carrier_freq_update::typeof(0.0Hz) = 0.0Hz
    vt_on::Bool = false
    frequency_lock::FrequencyLockIndicator = FrequencyLockIndicator()
    signal_combining_sums::SignalCombiningSums = SignalCombiningSums()
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
        code_discr_accs = map(_ -> (0, 0.0), sat.signals),
        carrier_discr_accs = map(_ -> (0, 0.0Hz), sat.signals),
    )
end

# The same accumulators emptied.
@inline _zeroed(accs::Tuple) = map(acc -> zero.(acc), accs)

# One reading added to a `(count, sum)` accumulator.
@inline _accumulated(acc::Tuple, reading) = acc .+ (1, reading)

function SatVectorPLLAndDLL(
    sat_vector_pll_and_dll::SatVectorPLLAndDLL{CA,CO,N};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_discr_accs::Maybe{NTuple{N,Tuple{Int,Float64}}} = nothing,
    code_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    carrier_discr_accs::Maybe{NTuple{N,Tuple{Int,typeof(0.0Hz)}}} = nothing,
    carrier_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    vt_on::Maybe{Bool} = nothing,
    frequency_lock::Maybe{FrequencyLockIndicator} = nothing,
    signal_combining_sums::Maybe{SignalCombiningSums} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    SatVectorPLLAndDLL{CA,CO,N}(
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
        isnothing(code_discr_accs) ? sat_vector_pll_and_dll.code_discr_accs :
        code_discr_accs,
        isnothing(code_freq_update) ? sat_vector_pll_and_dll.code_freq_update :
        code_freq_update,
        isnothing(carrier_discr_accs) ? sat_vector_pll_and_dll.carrier_discr_accs :
        carrier_discr_accs,
        isnothing(carrier_freq_update) ? sat_vector_pll_and_dll.carrier_freq_update :
        carrier_freq_update,
        isnothing(vt_on) ? sat_vector_pll_and_dll.vt_on : vt_on,
        isnothing(frequency_lock) ? sat_vector_pll_and_dll.frequency_lock : frequency_lock,
        isnothing(signal_combining_sums) ? sat_vector_pll_and_dll.signal_combining_sums :
        signal_combining_sums,
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

  - This estimator accumulates every signal's DLL / FLL discriminator
    outputs, one accumulator pair per signal, for the navigation filter to
    consume (and reset via
    [`reset_code_discr_accs!`](@ref) / [`reset_carrier_discr_accs!`](@ref)).
    The FLL readings are four-quadrant on a wiped-off prompt, see
    [The per-integration contract](@ref).
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

Bandwidths (auto by default) are sized and capped exactly as for
[`ConventionalPLLAndDLL`](@ref).

`combine_signals = true` combines the passengers' discriminators into the
driver's loops as for [`ConventionalPLLAndDLL`](@ref); under vector closure
into the PLL only.
"""
struct VectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter} <:
       AbstractDopplerEstimator
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)}
    combine_signals::Bool
end

function VectorPLLAndDLL(
    ::Type{CA} = ThirdOrderAssistedBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    combine_signals::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        combine_signals,
    )
end

# Kwarg-update constructor for tweaking bandwidths in place.
function VectorPLLAndDLL(
    pll_and_dll::VectorPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    combine_signals::Maybe{Bool} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
        isnothing(combine_signals) ? pll_and_dll.combine_signals : combine_signals,
    )
end

"""
$(SIGNATURES)

Build the per-satellite estimator state stored in a [`TrackedSat`](@ref) for a
satellite tracked under [`VectorPLLAndDLL`](@ref). Auto bandwidths are resolved
as for [`ConventionalPLLAndDLL`](@ref). New satellites start with
`vt_on = false` (the scalar fallback loop) until [`enable_vt!`](@ref) puts them
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
    SatVectorPLLAndDLL(
        sat,
        carrier_loop_filter,
        code_loop_filter;
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
    )
end

# Under vector closure the navigation filter owns the code loop and the FLL
# branch, so passengers are combined into the PLL only.
@inline _locally_closed_loops(state::SatVectorPLLAndDLL) =
    state.vt_on ? (pll = true, fll = false, dll = false) :
    (pll = true, fll = true, dll = true)

# Re-seed hook used by `reset_loop_filters!`, as for the conventional estimator
# (the carrier loop's staging restarts too), plus zeroing the discriminator
# accumulators and the NCO corrections: the current Dopplers already contain the
# last correction, so keeping it would apply it twice. Bandwidths and `vt_on`
# survive the reset.
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
        code_discr_accs = _zeroed(state.code_discr_accs),
        carrier_discr_accs = _zeroed(state.carrier_discr_accs),
        vt_on = state.vt_on,
    )
end

# The vector closure of one record, plugged into the shared driver fold by
# dispatch on the per-sat state; with `vt_on` unset it is the scalar closure. See
# "The per-integration contract" in docs/src/vector_tracking.md.
#
# The vector closure is not staged: it stays FLL-assisted and does not run the
# frequency lock indicator. Its carrier discriminators are picked, and `pending`
# combined, as in the scalar closure.
@inline function _close_loops(
    state::SatVectorPLLAndDLL,
    signal::AbstractGNSSSignal,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    code_doppler,
    sampling_frequency,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
    wiped_off::Bool,
    polarity::Int8,
    pending::SignalCombiningSums = SignalCombiningSums(),
)
    state.vt_on || return _close_scalar_loops(
        state,
        signal,
        correlator,
        previous_prompt,
        code_doppler,
        sampling_frequency,
        integration_time,
        carrier_bandwidth,
        code_bandwidth,
        wiped_off,
        polarity,
        pending,
    )
    pll_discriminator = _gated_mean(
        pll_disc(signal, correlator; polarity),
        _discriminator_weight(signal, integration_time),
        pending.pll,
        _TWO_QUADRANT_PLL_RANGE,
    )
    fll_discriminator = fll_disc(
        signal,
        correlator,
        previous_prompt,
        integration_time;
        four_quadrant = wiped_off,
    )
    dll_discriminator = dll_disc(signal, correlator, code_doppler, sampling_frequency)
    carrier_freq_update, carrier_loop_filter = calculate_carrier_frequency_update(
        state.carrier_loop_filter,
        pll_discriminator,
        state.carrier_freq_update,
        integration_time,
        carrier_bandwidth,
    )
    # The NCO-correction fields are not written back: they are owned by the
    # navigation filter, and overwriting them would clobber its value.
    carrier_freq_update,
    state.code_freq_update,
    SatVectorPLLAndDLL(
        state;
        carrier_loop_filter,
        # The driver's own readings, in its slot; the passengers' follow in the
        # driver fold.
        code_discr_accs = Base.setindex(
            state.code_discr_accs,
            _accumulated(first(state.code_discr_accs), dll_discriminator),
            1,
        ),
        carrier_discr_accs = Base.setindex(
            state.carrier_discr_accs,
            _accumulated_fll(
                first(state.carrier_discr_accs),
                fll_discriminator,
                previous_prompt,
            ),
            1,
        ),
    )
end

# A record without a previous prompt has no FLL reading (`fll_disc` reads 0 Hz)
# and is left out of the accumulator, as the scalar loop leaves it out.
@inline _accumulated_fll(acc::Tuple, fll_discriminator, previous_prompt::Complex) =
    iszero(previous_prompt) ? acc : _accumulated(acc, fll_discriminator)

@inline _with_loop_state(state::SatVectorPLLAndDLL; kwargs...) =
    SatVectorPLLAndDLL(state; kwargs...)

_carrier_phase_polarity(::SatVectorPLLAndDLL, sat::TrackedSat) =
    _sync_polarity(first(sat.signals).signal, first(sat.signals).bit_buffer, sat.prn)

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
    # See the conventional `estimate_dopplers_and_filter_prompt`.
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
        track_state.doppler_estimator.combine_signals,
    )
    return track_state
end

# Shared in-place walker for the vector-tracking state managers below: replace
# each satellite's `doppler_estimator_state` with
# `f(sat, sat.doppler_estimator_state, args...)`, which must return the same
# concrete `SatVectorPLLAndDLL` type.
#
# `_map_vt_states!` walks every group, `_map_vt_states_in_group!` only the
# addressed one. They are distinct names, not overloads on the group position,
# so a threaded `Symbol`/`Integer` argument (including the `vt_on` Bool, since
# `Bool <: Integer`) is never mistaken for a group selector.
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

# Shared by `enable_vt!` / `disable_vt!`: write `vt_on` to the addressed PRNs.
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

Reset the code (DLL) discriminator accumulators, one per signal, of every
satellite in the vector loop — called by the navigation filter after it has consumed the
accumulated values via [`mean_code_discr`](@ref). Satellites outside
the vector loop are left untouched (they accumulate nothing).

The one-argument form walks every group; pass a `group` (Symbol or index) to
reset a single group's satellites. Mutates `track_state` in place and returns
it.
"""
function reset_code_discr_accs!(track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL})
    _map_vt_states!(_reset_sat_code_discr_accs, track_state)
end

function reset_code_discr_accs!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
)
    _map_vt_states_in_group!(_reset_sat_code_discr_accs, track_state, group)
end

_reset_sat_code_discr_accs(_sat, state) =
    state.vt_on ?
    SatVectorPLLAndDLL(state; code_discr_accs = _zeroed(state.code_discr_accs)) : state

"""
$(SIGNATURES)

Reset the carrier (FLL) discriminator accumulators, one per signal, of every
satellite in the vector loop — called by the navigation filter after it has consumed the
accumulated values via [`mean_carrier_discr`](@ref). Addressed exactly
like [`reset_code_discr_accs!`](@ref): every group, or only `group` when
one is given. Mutates `track_state` in place and returns it.
"""
function reset_carrier_discr_accs!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
)
    _map_vt_states!(_reset_sat_carrier_discr_accs, track_state)
end

function reset_carrier_discr_accs!(
    track_state::TrackState{<:SignalGroups,<:VectorPLLAndDLL},
    group::Union{Symbol,Integer,Val},
)
    _map_vt_states_in_group!(_reset_sat_carrier_discr_accs, track_state, group)
end

_reset_sat_carrier_discr_accs(_sat, state) =
    state.vt_on ?
    SatVectorPLLAndDLL(state; carrier_discr_accs = _zeroed(state.carrier_discr_accs)) :
    state

"""
$(SIGNATURES)

Mean DLL (code) discriminator of one signal accumulated since the last
[`reset_code_discr_accs!`](@ref), in chips, or `nothing` if nothing has been
accumulated yet (`count == 0`). Use this rather than dividing `code_discr_accs` by
hand, so consumers do not depend on how the accumulators are stored.

Every signal of a satellite has its own accumulator, filled with that signal's
raw readings. A passenger's DLL reads its received code against the satellite's
one shared replica code phase, so its reading includes the inter-signal bias
(see [`set_group_delay!`](@ref)); removing it and fusing the signals is the
navigation filter's job.

Addressed like the other per-signal accessors: from a [`TrackState`](@ref) or a
[`TrackedSat`](@ref) with a trailing signal selector, an index or a signal type,
which may be omitted only for a single-signal satellite. The
`SatVectorPLLAndDLL` form takes an index.

```julia
mean_code_discr(track_state, :e1, 11, GalileoE1B)        # by signal type
mean_code_discr(sat, 2)                                  # by index
mean_code_discr(get_doppler_estimator_state(sat))        # single-signal satellite
```
"""
function mean_code_discr(state::SatVectorPLLAndDLL, signal_index::Integer)
    count, discr_sum = state.code_discr_accs[signal_index]
    count == 0 ? nothing : discr_sum / count
end

"""
$(SIGNATURES)

Mean FLL (carrier) discriminator of one signal accumulated since the last
[`reset_carrier_discr_accs!`](@ref), in Hz, or `nothing` if nothing has been
accumulated yet (`count == 0`). The carrier counterpart to
[`mean_code_discr`](@ref), addressed the same way. A signal's reading is
four-quadrant where its prompt is wiped off, see [The per-integration
contract](@ref). A record without a previous prompt (the first after
[`add_satellite!`](@ref), [`reset_loop_filters!`](@ref) or a pilot's
secondary-code sync, or one whose length differs from the previous record's) has
no FLL reading and is not counted.
"""
function mean_carrier_discr(state::SatVectorPLLAndDLL, signal_index::Integer)
    count, discr_sum = state.carrier_discr_accs[signal_index]
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

# The satellite form turns a signal type into an index, since the estimator
# state holds the accumulators but not the signals; the `TrackState` forms
# forward to it like the other per-signal accessors.
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
            sig::Union{Integer,Type{<:AbstractGNSSSignal}},
        ) = $fn(get_sat_state(s, group, sat_id), sig)
    end
end

# The passengers' records come through the driver fold while `vt_on`, each adding
# its raw DLL and FLL readings to that passenger's own accumulators.
@inline _measures_passengers(state::SatVectorPLLAndDLL) = state.vt_on
@inline _passenger_measurement_accs(state::SatVectorPLLAndDLL, passengers::Tuple) = map(
    (_, code_acc, carrier_acc) -> (code_acc, carrier_acc),
    passengers,
    Base.tail(state.code_discr_accs),
    Base.tail(state.carrier_discr_accs),
)
@inline _with_passenger_measurements(state::SatVectorPLLAndDLL, accs::Tuple) =
    SatVectorPLLAndDLL(
        state;
        code_discr_accs = (first(state.code_discr_accs), map(first, accs)...),
        carrier_discr_accs = (first(state.carrier_discr_accs), map(last, accs)...),
    )

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
Addressed like [`set_code_freq_updates!`](@ref). Mutates `track_state` in place
and returns it.
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
