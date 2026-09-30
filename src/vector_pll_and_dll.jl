"""
Per-satellite state for the vector PLL and DLL Doppler estimator
([`VectorPLLAndDLL`](@ref)).

On top of the conventional loop-filter state it carries the
vector-tracking (VT) interface to an external navigation filter
(e.g. GNSSReceiver.jl's VDFLL):

  - `code_discr_acc` / `carrier_discr_acc`: `(count, sum)` accumulators of the
    raw DLL discriminator (chips) and the raw FLL discriminator (Hz) since the
    navigation filter last read and reset them
    ([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
    **One pair per signal and per loop**, in `sat.signals` order, so a
    multi-signal satellite hands the filter every component's measurement
    rather than the driver's alone — read them with [`mean_code_discr`](@ref) /
    [`mean_carrier_discr`](@ref), which take the same signal selector as every
    other per-signal accessor. Collected while `vt_on`, whatever the group's
    `discriminator_combining` flag says.

    A count and a sum, with no per-record timing alongside them: every record
    folded in adds one to the count and its discriminator output to the sum.
    Which interval the mean covers is therefore the consumer's to know —
    `count` records of the coherent integration time in force — and that holds
    while the length does not change under it. See [`VectorPLLAndDLL`](@ref).

  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.

  - `vt_on`: whether the navigation filter controls this satellite's NCOs.
    While `false` the satellite runs a conventional (scalar) PLL/DLL as a
    fallback and nothing is accumulated. Set by [`enable_vt!`](@ref) /
    [`disable_vt!`](@ref).
"""
@kwdef struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    # No defaults, deliberately: a default would have to pick an `N`, and a
    # one-slot default silently truncates a multi-signal satellite — slot 2's
    # measurements would have nowhere to go and `mean_code_discr(state, 2)` would
    # throw. `N` is the satellite's signal count, so it comes from the
    # satellite: use the `SatVectorPLLAndDLL(sat, …)` constructor below, or pass
    # `_zero_code_discr_acc(sat)` / `_zero_carrier_discr_acc(sat)` here.
    code_discr_acc::NTuple{N,Tuple{Int,Float64}}
    code_freq_update::typeof(0.0Hz) = 0.0Hz
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}}
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
        code_discr_acc = _zero_code_discr_acc(sat),
        carrier_discr_acc = _zero_carrier_discr_acc(sat),
    )
end

# One zeroed accumulator per signal of the satellite. The `Val` keeps `N` a
# compile-time constant, so the state's type is fixed per group — every
# satellite of a group shares its signal-tuple shape.
@inline _zero_code_discr_acc(sat::TrackedSat) = ntuple(_ -> (0, 0.0), _num_signals_val(sat))
@inline _zero_carrier_discr_acc(sat::TrackedSat) =
    ntuple(_ -> (0, 0.0Hz), _num_signals_val(sat))

function SatVectorPLLAndDLL(
    sat_vector_pll_and_dll::SatVectorPLLAndDLL{CA,CO,N};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_discr_acc::Maybe{NTuple{N,Tuple{Int,Float64}}} = nothing,
    code_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    carrier_discr_acc::Maybe{NTuple{N,Tuple{Int,typeof(0.0Hz)}}} = nothing,
    carrier_freq_update::Maybe{typeof(0.0Hz)} = nothing,
    vt_on::Maybe{Bool} = nothing,
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
    [`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)) —
    one accumulator per *signal*, so a multi-signal satellite hands the filter
    every component's raw measurement and fuses none of them itself.
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

Collecting those measurements follows `vt_on` alone: a satellite in the vector
loop hands the filter every signal's measurement whether or not its group
combines anything for its own loops.

!!! note "A `(count, sum)` pair maps to an interval only at one length"

    `count × integration_time` is the interval the mean covers only while every
    record in it ran at the same coherent integration time. Calling
    [`set_preferred_num_code_blocks_to_integrate!`](@ref) while `vt_on` is set
    leaves the mean spanning two lengths with nothing in the pair to say so, so
    reset the accumulators across such a change or make it with the satellite
    out of the vector loop. The vector-tracking manual works this through.

With the group's `discriminator_combining = true` a multi-signal satellite folds
its passengers' coincident records into the driver's loop update exactly as
[Multi-signal discriminator combining](@ref Multi-signal-discriminator-combining)
describes for the conventional estimator — same coincidence rule, same weighting,
same de-rotation onto the driver's phase frame, same equal-length assumption.
Which loops it reaches follows `vt_on`, because that is what decides which loops
this package still closes:

  - **`vt_on = false`.** Every loop is local, so this is scalar tracking and all
    three discriminators combine — bit-for-bit the conventional estimator's
    behaviour, [`set_group_delay!`](@ref) included.
  - **`vt_on = true`.** Only the carrier phase loop is still closed here, and
    only it combines. The code and carrier frequency loops are the navigation
    filter's, and every signal's discriminators reach it as that signal's
    **own, raw** accumulator ([`mean_code_discr`](@ref) /
    [`mean_carrier_discr`](@ref)), which is where they should be fused: the
    filter holds each signal's measured C/N₀, tap spacing and integration
    length, so it can weigh them better than a nominal power split can, and a
    mean formed here would reach it carrying one signal's variance. The group
    delay goes with the loop that reads it, so it is **not** applied to these
    accumulators; whatever inter-signal bias the filter needs it applies itself,
    from [`get_group_delay`](@ref) or from its own tables.

Ranging is on the shared `code_phase` of `signals[1]` in either mode.

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
until [`enable_vt!`](@ref) puts them into the vector loop.
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
        code_discr_acc = _zero_code_discr_acc(sat),
        carrier_discr_acc = _zero_carrier_discr_acc(sat),
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Every signal's accumulator slot is zeroed;
# per-sat bandwidth overrides and the `vt_on` flag survive the reset.
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
        code_discr_acc = map(_ -> (0, 0.0), state.code_discr_acc),
        carrier_discr_acc = map(_ -> (0, 0.0Hz), state.carrier_discr_acc),
        vt_on = state.vt_on,
    )
end

# What this estimator does with one record's discriminators — the one
# `_close_loops` method that is not the default.
#
# With `vt_on` unset this satellite is running its scalar fallback, and the
# fallback *is* the conventional loop closure, so it defers to the same
# `_close_scalar_loops` rather than restating its two filter calls.
#
# While `vt_on`, the navigation filter's NCO corrections drive the loops:
# `code_freq_update` directly, the code loop filter bypassed and holding its
# state, and `carrier_freq_update` through the FLL branch of the carrier loop
# filter, the PLL branch still running on this satellite's own phase
# discriminator — combined with the coincident passenger records' where the
# group combines, since the carrier phase loop is the one loop still closed
# here. The fold hands in combined `dll` and `fll` values too, and they are
# ignored: those two loops are the navigation filter's, fed by every signal's
# raw measurement instead (see `PerSignalMeasurements`). Only the FLL-assisted filter has an input path for that branch;
# on any other `_filter_carrier_loop` runs on the PLL discriminator alone, so the
# vector carrier closure degrades to PLL-only.
@inline function _close_loops(
    state::SatVectorPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    discriminators,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
)
    state.vt_on || return _close_scalar_loops(
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
        state.carrier_freq_update,
        integration_time,
        carrier_bandwidth,
    )
    (;
        carrier_freq_update,
        code_freq_update = state.code_freq_update,
        carrier_loop_filter,
        code_loop_filter,
    )
end

"""
The per-signal navigation measurements a satellite in the vector loop (`vt_on`)
collects during a chunk: one `(count, sum)` pair per signal for the raw DLL
discriminator and one for the raw FLL discriminator, in `sat.signals` order,
growing until the navigation filter reads and resets them
([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).

Every signal accumulates its **own** here and never a combination across
signals — fusing them is the navigation filter's job.
"""
struct PerSignalMeasurements{N}
    code_discr_acc::NTuple{N,Tuple{Int,Float64}}
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}}
end

# Measurements are collected only while the navigation filter has this
# satellite, and then for every signal whatever the group's combining flag says:
# `vt_on` decides whether the filter is handed the measurements, the flag only
# what the loops this package still closes are fed. Outside the vector loop
# there is nothing to collect, and the fold runs with `nothing` exactly as the
# conventional estimator's does.
@inline _with_measurements(f::F, state::SatVectorPLLAndDLL) where {F} =
    state.vt_on ? f(PerSignalMeasurements(state.code_discr_acc, state.carrier_discr_acc)) :
    f(nothing)

# Fold one signal's discriminator into its own `(count, sum)`. `Val(N)` keeps
# the rebuilt tuple's length a compile-time constant, so this is
# allocation-free and the slot write costs a predicted branch per signal.
@inline _accumulate_one_discr(accs::NTuple{N,Any}, slot::Integer, value) where {N} =
    ntuple(k -> k == slot ? (accs[k][1] + 1, accs[k][2] + value) : accs[k], Val(N))

# One record — any signal's, the driver's included — into its own slot. The raw
# discriminators, with no group-delay referral and no de-rotation: the group delay
# goes with the combined code loop, and under `vt_on` the code loop is the
# navigation filter's (see `VectorPLLAndDLL`); `dll_disc` and `fll_disc` are
# rotation-invariant anyway. Nothing here withholds a value, so `count` is the
# record count. The one record the fold skips is a passenger's correlated with a
# pre-sync replica on a secondary-coded signal, and it skips both slots
# together.
@inline function _collect_measurements(
    m::PerSignalMeasurements,
    signal::AbstractGNSSSignal,
    output::CorrelatorOutput,
    filtered_correlator,
    previous_prompt,
    ctx::RecordFoldContext,
)
    integration_time = output.integrated_samples / ctx.sampling_frequency
    # `dll_disc` is fed the chunk-fixed `sat.code_doppler` — the code Doppler
    # that generated this chunk's replicas — for every record.
    dll = dll_disc(signal, filtered_correlator, ctx.code_doppler, ctx.sampling_frequency)
    fll = fll_disc(signal, filtered_correlator, previous_prompt, integration_time)
    PerSignalMeasurements(
        _accumulate_one_discr(m.code_discr_acc, ctx.signal_index, dll),
        _accumulate_one_discr(m.carrier_discr_acc, ctx.signal_index, fll),
    )
end

# Neither NCO-correction field (`code_freq_update` / `carrier_freq_update`) is
# written back: both are owned by the navigation filter (set via
# `set_code_freq_updates!` / `set_carrier_freq_updates!`) and read as this fold's
# inputs, so overwriting them with a per-record loop output would clobber the
# navigation filter's value between two of its calls.
@inline _park_state(
    state::SatVectorPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    ::Nothing,
) = SatVectorPLLAndDLL(state; carrier_loop_filter, code_loop_filter)

@inline _park_state(
    state::SatVectorPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    m::PerSignalMeasurements,
) = SatVectorPLLAndDLL(
    state;
    carrier_loop_filter,
    code_loop_filter,
    code_discr_acc = m.code_discr_acc,
    carrier_discr_acc = m.carrier_discr_acc,
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

Reset the code (DLL) discriminator accumulators of every satellite in the
vector loop — every signal's slot — called by the navigation filter after it
has consumed the accumulated values via [`mean_code_discr`](@ref). Satellites
outside the vector loop are left untouched (they accumulate nothing).

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
    SatVectorPLLAndDLL(state; code_discr_acc = map(_ -> (0, 0.0), state.code_discr_acc)) :
    state

"""
$(SIGNATURES)

Reset the carrier (FLL) discriminator accumulators of every satellite in the
vector loop — every signal's slot — called by the navigation filter after it has consumed the
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
    SatVectorPLLAndDLL(
        state;
        carrier_discr_acc = map(_ -> (0, 0.0Hz), state.carrier_discr_acc),
    ) : state

"""
$(SIGNATURES)

Mean DLL (code) discriminator of one signal, accumulated since the last
[`reset_code_discr_acc!`](@ref), in chips — or `nothing` if nothing has been
accumulated yet (the `(count, sum)` accumulator's `count` is 0). This is the
single place the accumulator's averaging convention lives — read it here rather
than dividing `code_discr_acc` by hand, so a change to how the accumulator is
stored can't silently diverge between consumers.

Every signal of a satellite accumulates its **own** discriminator, never a
combination across signals: fusing them is the navigation filter's job, and
[`VectorPLLAndDLL`](@ref) says why.

The value is the signal's **raw** code discriminator, on that signal's own code
phase: no [`get_group_delay`](@ref) difference is subtracted. This package
applies that difference only where it closes the combined code loop itself, and
under `vt_on` that loop is the consumer's — so the inter-signal bias is the
consumer's too, to apply with whatever it already applies downstream when it
forms a pseudorange. [`get_group_delay`](@ref) hands it the per-signal value if
it wants Tracking's.

Nothing is withheld for want of a group delay, so `count` is the number of
records folded since the reset. The one record a fold does skip is a
passenger's correlated with a replica a mid-chunk bit/secondary-code sync had
not yet corrected — it carries no usable delay or phase — and that can only
happen on the fold that *detects* sync, which is long before a navigation filter
admits the satellite to the vector loop. Both loops' counts skip it together,
so they never diverge.

Addressed like every other per-signal accessor ([`estimate_cn0`](@ref),
[`get_group_delay`](@ref)): from a [`TrackState`](@ref) or a
[`TrackedSat`](@ref) with a trailing signal selector — an index or a signal type
— which may be omitted only when the satellite tracks one signal. The
`SatVectorPLLAndDLL` form takes an index, since the estimator state alone cannot
resolve a signal type.

```julia
mean_code_discr(track_state, :gps_l1, 7, GPSL1C_D)   # by signal type
mean_code_discr(sat, 2)                              # by index
mean_code_discr(get_doppler_estimator_state(sat))    # single-signal satellite
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
[`mean_code_discr`](@ref), addressed the same way and per signal for the same
reason — except that the signals of one satellite share a carrier, so these
measure one Doppler and need no inter-signal correction to be fused.
"""
function mean_carrier_discr(state::SatVectorPLLAndDLL, signal_index::Integer)
    count, discr_sum = state.carrier_discr_acc[signal_index]
    count == 0 ? nothing : discr_sum / count
end

# A satellite tracking one signal has one slot, so its accumulators need no
# selector; one tracking several is refused exactly as `_find_signal` refuses an
# unqualified per-signal read, because silently answering with the driver's is
# the mistake per-signal accumulation exists to remove.
@inline _sole_signal_index(::NTuple{1,Any}) = 1
@inline _sole_signal_index(::Tuple) = _throw_needs_signal_selector()

mean_code_discr(state::SatVectorPLLAndDLL) =
    mean_code_discr(state, _sole_signal_index(state.code_discr_acc))
mean_carrier_discr(state::SatVectorPLLAndDLL) =
    mean_carrier_discr(state, _sole_signal_index(state.carrier_discr_acc))

# The satellite-level rung of the accessor ladder, which is what turns a signal
# *type* into a slot: the estimator state holds the accumulators but not the
# signals, so only the satellite can resolve one. The `TrackState` rungs are
# generated in tracking_state.jl with the other per-signal accessors.
for fn in (:mean_code_discr, :mean_carrier_discr)
    @eval $fn(sat::TrackedSat, sel...) =
        $fn(get_doppler_estimator_state(sat), _signal_index(sat.signals, sel...))
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
