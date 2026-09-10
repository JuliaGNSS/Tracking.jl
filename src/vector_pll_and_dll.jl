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
    Only accumulated while `vt_on`. **One pair per signal**, in `sat.signals`
    order, so a multi-signal satellite hands the filter every component's
    measurement rather than the driver's alone — read them with
    [`mean_code_discr`](@ref) / [`mean_carrier_discr`](@ref), which take the
    same signal selector as every other per-signal accessor.
  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.
  - `vt_on`: whether the navigation filter controls this satellite's NCOs.
    While `false` the satellite runs a conventional (scalar) PLL/DLL as a
    fallback and nothing is accumulated. Set by [`enable_vt!`](@ref) /
    [`disable_vt!`](@ref).

`signal_combining` and `pending_discriminators` carry multi-signal combining,
whose reach follows `vt_on`: every loop while this satellite runs its own scalar
fallback, the carrier phase loop alone once the navigation filter has the other
two (see [`VectorPLLAndDLL`](@ref)). The flag is seeded from the shared
[`VectorPLLAndDLL`](@ref) by [`init_estimator_state`](@ref) and owned here
afterwards, like the bandwidths, and survives [`reset_loop_filters!`](@ref); the
accumulator does not, and is dropped on every `vt_on` transition as well — sums
gathered under one mode's reach must not reach a loop update under the other's.
"""
@kwdef struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    signal_combining::Bool = false
    pending_discriminators::DiscriminatorAccumulator = zero(DiscriminatorAccumulator)
    code_discr_acc::NTuple{N,Tuple{Int,Float64}} = ((0, 0.0),)
    code_freq_update::typeof(0.0Hz) = 0.0Hz
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}} = ((0, 0.0Hz),)
    carrier_freq_update::typeof(0.0Hz) = 0.0Hz
    vt_on::Bool = false
end

function SatVectorPLLAndDLL(
    sat::TrackedSat,
    carrier_loop_filter::CA,
    code_loop_filter::CO;
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz,
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz,
    signal_combining::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        signal_combining,
        code_discr_acc = _zero_code_discr_accs(sat),
        carrier_discr_acc = _zero_carrier_discr_accs(sat),
    )
end

# One zeroed accumulator per signal of the satellite. The `Val` keeps `N` a
# compile-time constant, so the state's type is fixed per group — every
# satellite of a group shares its signal-tuple shape.
@inline _zero_code_discr_accs(sat::TrackedSat) =
    ntuple(_ -> (0, 0.0), _num_signals_val(sat))
@inline _zero_carrier_discr_accs(sat::TrackedSat) =
    ntuple(_ -> (0, 0.0Hz), _num_signals_val(sat))

function SatVectorPLLAndDLL(
    sat_vector_pll_and_dll::SatVectorPLLAndDLL{CA,CO,N};
    carrier_loop_filter::Maybe{CA} = nothing,
    code_loop_filter::Maybe{CO} = nothing,
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    signal_combining::Maybe{Bool} = nothing,
    pending_discriminators::Maybe{DiscriminatorAccumulator} = nothing,
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
        isnothing(signal_combining) ? sat_vector_pll_and_dll.signal_combining :
        signal_combining,
        isnothing(pending_discriminators) ? sat_vector_pll_and_dll.pending_discriminators :
        pending_discriminators,
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
    every component's measurement and fuses none of them itself.
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

With `signal_combining = true` a multi-signal satellite folds its signals'
discriminators into one minimum-variance weighted mean before the loop filters
see it, exactly as
[Multi-signal discriminator combining](@ref Multi-signal-discriminator-combining)
describes for the conventional estimator — same weighting, same de-rotation onto
the driver's phase frame, same precondition that `signals[1]` must be the group's
longest-integrating signal.

Which loops it reaches follows `vt_on`, because that is what decides which loops
this package still closes:

  - **`vt_on = false`.** Every loop is local, so this is scalar tracking and all
    three discriminators combine — bit-for-bit the conventional estimator's
    behaviour, [`set_differential_group_delay!`](@ref) included, since the code
    loop needs each passenger referred to the driver's code phase.
  - **`vt_on = true`.** Only the carrier phase loop is still closed here, and
    only it combines. The code and carrier frequency loops are the navigation
    filter's, and every signal's discriminators reach it as that signal's
    **own** accumulator ([`mean_code_discr`](@ref) / [`mean_carrier_discr`](@ref)),
    which is where they should be fused: the filter holds each signal's measured
    C/N₀, tap spacing and integration length, so it can weigh them better than
    a nominal power split can, and a mean formed here would reach it carrying
    one signal's variance. For the same reason it, and not this package, applies
    the differential group delay to the code measurements — referring them to
    one ranging datum belongs with the consumer that ranges on them.

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
    signal_combining::Bool
end

function VectorPLLAndDLL(
    ::Type{CA} = ThirdOrderAssistedBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    signal_combining::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        signal_combining,
    )
end

# Kwarg-update constructor for tweaking bandwidths / combining in place.
function VectorPLLAndDLL(
    pll_and_dll::VectorPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    signal_combining::Maybe{Bool} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
        isnothing(signal_combining) ? pll_and_dll.signal_combining : signal_combining,
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
        signal_combining = estimator.signal_combining,
        code_discr_acc = _zero_code_discr_accs(sat),
        carrier_discr_acc = _zero_carrier_discr_accs(sat),
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Per-sat bandwidth overrides, the
# `signal_combining` flag and the `vt_on` flag survive the reset; the pending
# passenger carrier phase discriminators do not — they are pre-reset history,
# like the filter integrators and each signal's FLL prompt.
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
        signal_combining = state.signal_combining,
        code_discr_acc = map(_ -> (0, 0.0), state.code_discr_acc),
        carrier_discr_acc = map(_ -> (0, 0.0Hz), state.carrier_discr_acc),
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

# What this estimator does with one record's three discriminators, shared by the
# driver-only fold and the combining one so the two cannot disagree
# about what `vt_on` changes.
#
# While `vt_on`, the navigation filter's NCO corrections drive the loops —
# `code_freq_update` directly (the code loop filter bypassed, holding its state)
# and `carrier_freq_update` through the FLL branch of the carrier loop filter,
# the PLL branch still running on this satellite's own phase discriminator. With
# `vt_on` unset it is a conventional FLL-assisted scalar PLL/DLL.
#
# The discriminators arrive already formed because only `pll_discriminator` is
# ever a combination of several signals.
@inline function _close_vector_loops(
    state::SatVectorPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    pll_discriminator,
    fll_discriminator,
    dll_discriminator,
    integration_time,
    carrier_bandwidth,
    code_bandwidth,
)
    if state.vt_on
        code_freq_update = state.code_freq_update
        fll_input = state.carrier_freq_update
    else
        code_freq_update, code_loop_filter = filter_loop(
            code_loop_filter,
            dll_discriminator,
            integration_time,
            code_bandwidth,
        )
        fll_input = fll_discriminator
    end
    carrier_freq_update, carrier_loop_filter = _filter_vector_carrier_loop(
        carrier_loop_filter,
        pll_discriminator,
        fll_input,
        integration_time,
        carrier_bandwidth,
    )
    (; carrier_freq_update, code_freq_update, carrier_loop_filter, code_loop_filter)
end

# Fold one signal's discriminator into its own `(count, sum)`. `Val(N)` keeps
# the rebuilt tuple's length a compile-time constant, so this is
# allocation-free and the slot write costs a predicted branch per signal.
@inline _accumulate_one_discr(accs::NTuple{N,Any}, slot::Integer, value) where {N} =
    ntuple(k -> k == slot ? (accs[k][1] + 1, accs[k][2] + value) : accs[k], Val(N))

# Process the estimator-driver signal (signals[1]) under vector tracking.
# Mirrors the conventional driver fold (dispatching on the per-sat state type
# plugs this into the shared `_update_tracked_sat_doppler`): fold over every
# `CorrelatorOutput` collected during this chunk, threading the loop filters,
# the FLL `previous_prompt`, and the discriminator accumulators across records,
# then return the *last* record's carrier/code Doppler (the NCO is written once
# per chunk). With no outputs the Doppler holds.
#
# Every discriminator here is the driver's own; `_close_vector_loops` does the
# rest. This is the path a satellite takes with combining off, and the one a
# single-signal satellite takes either way.
@inline function _process_estimator_driver_signal(
    tracked_signal::TrackedSignal,
    sat::TrackedSat,
    pll_and_dll_state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise_density,
    noise_density_ready::Bool,
)
    outputs = tracked_signal.correlator_outputs
    if isempty(outputs)
        return tracked_signal, pll_and_dll_state, sat.carrier_doppler, sat.code_doppler
    end
    signal = tracked_signal.signal
    tsig = tracked_signal
    carrier_loop_filter = pll_and_dll_state.carrier_loop_filter
    code_loop_filter = pll_and_dll_state.code_loop_filter
    code_discr_acc = pll_and_dll_state.code_discr_acc
    carrier_discr_acc = pll_and_dll_state.carrier_discr_acc
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    found_before_fold = has_bit_or_secondary_code_been_found(tsig.bit_buffer)
    @inbounds for k in eachindex(outputs)
        # Same per-record prologue as the conventional folds, bandwidth handling
        # included: the carrier's per-primary-period reference scaled by 1/N,
        # the DLL's absolute bandwidth capped by the same stability product.
        folded = _advance_driver_record(
            tsig,
            outputs[k],
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            found_before_fold,
            pll_and_dll_state.carrier_loop_filter_bandwidth,
            pll_and_dll_state.code_loop_filter_bandwidth,
        )
        tsig = folded.tracked_signal
        filtered_correlator = folded.filtered_correlator
        integration_time = folded.integration_time

        pll_discriminator = pll_disc(signal, filtered_correlator)
        fll_discriminator =
            fll_disc(signal, filtered_correlator, folded.previous_prompt, integration_time)
        # `dll_disc` is fed the chunk-fixed `sat.code_doppler` — the code Doppler
        # that generated this chunk's replicas — for every record.
        dll_discriminator =
            dll_disc(signal, filtered_correlator, sat.code_doppler, sampling_frequency)

        if pll_and_dll_state.vt_on
            # This signal is the driver, hence slot 1 — the only slot a
            # single-signal satellite has, and the only one this fold fills:
            # the passengers' are the combining fold's business.
            code_discr_acc = _accumulate_one_discr(code_discr_acc, 1, dll_discriminator)
            carrier_discr_acc =
                _accumulate_one_discr(carrier_discr_acc, 1, fll_discriminator)
        end
        closed = _close_vector_loops(
            pll_and_dll_state,
            carrier_loop_filter,
            code_loop_filter,
            pll_discriminator,
            fll_discriminator,
            dll_discriminator,
            integration_time,
            folded.carrier_bandwidth,
            folded.code_bandwidth,
        )
        carrier_loop_filter = closed.carrier_loop_filter
        code_loop_filter = closed.code_loop_filter

        carrier_doppler, code_doppler = aid_dopplers(
            signal,
            pll_and_dll_state.init_carrier_doppler,
            pll_and_dll_state.init_code_doppler,
            closed.carrier_freq_update,
            closed.code_freq_update,
        )
    end
    empty!(outputs)
    # Neither NCO-correction field (`code_freq_update` / `carrier_freq_update`)
    # is written back: both are owned by the navigation filter (set via
    # `set_code_freq_updates!` / `set_carrier_freq_updates!`) and read as this
    # fold's inputs, so overwriting them with a per-record loop output would
    # clobber the navigation filter's value between two of its calls.
    new_doppler_estimator_state = SatVectorPLLAndDLL(
        pll_and_dll_state;
        carrier_loop_filter,
        code_loop_filter,
        code_discr_acc,
        carrier_discr_acc,
    )
    return tsig, new_doppler_estimator_state, carrier_doppler, code_doppler
end

# ---------------------------------------------------------------------------
# Multi-signal discriminator combining
# ---------------------------------------------------------------------------

"""
The passenger accumulator of a **vector-closed** satellite: the carrier phase
loop combines, and every signal's code and carrier frequency discriminators
accumulate in its own `(count, sum)` for the navigation filter to fuse — see
[`VectorPLLAndDLL`](@ref) for why the split falls there.

A marker wrapping the same [`DiscriminatorAccumulator`](@ref) rather than a
narrower struct, because which loops a satellite combines follows `vt_on` and
therefore changes at run time, while the field the sums park in cannot change
type. The fold is written once against `AbstractDiscriminatorAccumulator` and
specialises on this, so a vector-closed satellite forms no weight and no noise
gain for the two loops it does not combine.

They only ride along for the duration of a chunk's fold before parking back on
the `SatVectorPLLAndDLL`, and a driver record leaves them alone: it closes a
*loop* window, while they answer to the navigation filter's reading cadence
([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
"""
struct VectorClosureAccumulator{N} <: AbstractDiscriminatorAccumulator
    combined::DiscriminatorAccumulator
    code_discr_accs::NTuple{N,Tuple{Int,Float64}}
    carrier_discr_accs::NTuple{N,Tuple{Int,typeof(1.0Hz)}}
end

@inline _sums(acc::VectorClosureAccumulator) = acc.combined

# A driver record closes the combining window, not the per-signal accumulators.
@inline _zero_window(acc::VectorClosureAccumulator) = VectorClosureAccumulator(
    zero(DiscriminatorAccumulator),
    acc.code_discr_accs,
    acc.carrier_discr_accs,
)

# Where the per-signal accumulators end up when the chunk is done. In the scalar
# fallback nothing touched them, so the state's own tuples pass through.
@inline _code_discr_accs(acc::VectorClosureAccumulator, ::SatVectorPLLAndDLL) =
    acc.code_discr_accs
@inline _carrier_discr_accs(acc::VectorClosureAccumulator, ::SatVectorPLLAndDLL) =
    acc.carrier_discr_accs
@inline _code_discr_accs(::DiscriminatorAccumulator, state::SatVectorPLLAndDLL) =
    state.code_discr_acc
@inline _carrier_discr_accs(::DiscriminatorAccumulator, state::SatVectorPLLAndDLL) =
    state.carrier_discr_acc

# Fold one signal's own two discriminators into its slot. The no-op method is
# the scalar fallback, where the navigation filter reads nothing.
@inline _accumulate_own_discrs(acc::DiscriminatorAccumulator, ::Integer, _, _) = acc
@inline _accumulate_own_discrs(
    acc::VectorClosureAccumulator,
    slot::Integer,
    dll_discriminator,
    fll_discriminator,
) = VectorClosureAccumulator(
    acc.combined,
    _accumulate_one_discr(acc.code_discr_accs, slot, dll_discriminator),
    _accumulate_one_discr(acc.carrier_discr_accs, slot, fll_discriminator),
)

# One passenger record: its carrier phase measurement into the combination, its
# other two into its own accumulator. Nothing reads the differential group delay —
# referring a code measurement to one ranging datum is the navigation filter's
# job here, and it is the party that knows which datum it ranges on.
#
# The pre-sync gate excludes an invalidated record from *both*, and for one
# reason: a replica the sync had not yet corrected carries no usable phase or
# delay information, whether the value is about to be averaged or shipped.
@inline function _accumulate_passenger_discriminators(
    acc::VectorClosureAccumulator,
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    filtered_correlator,
    previous_prompt,
    ctx::PassengerFoldContext,
    correlated_pre_sync::Bool,
)
    signal = tracked_signal.signal
    _passenger_record_invalidated_by_sync(signal, correlated_pre_sync) && return acc
    combined = _accumulate_pll(
        acc.combined,
        _record_pll_contribution(
            signal,
            filtered_correlator,
            output.integrated_samples,
            ctx.rotation,
        )...,
    )
    integration_time = output.integrated_samples / ctx.sampling_frequency
    _accumulate_own_discrs(
        VectorClosureAccumulator(combined, acc.code_discr_accs, acc.carrier_discr_accs),
        ctx.signal_index,
        dll_disc(signal, filtered_correlator, ctx.code_doppler, ctx.sampling_frequency),
        fll_disc(signal, filtered_correlator, previous_prompt, integration_time),
    )
end

# Carrier phase combined, the other two this signal's own — and formed by the
# plain discriminators rather than by `_record_contribution`, so no weight, no
# noise gain and no group delay is computed for a value nothing will weigh. The
# conventional estimator's method, which combines all three, is in
# conventional_pll_and_dll.jl.
@inline function _driver_record_discriminators(
    marked::VectorClosureAccumulator,
    signal::AbstractGNSSSignal,
    folded,
    integrated_samples::Integer,
    code_doppler,
    sampling_frequency,
)
    acc = marked.combined
    own_pll, own_pll_weight = _record_pll_contribution(
        signal,
        folded.filtered_correlator,
        integrated_samples,
        one(ComplexF64),
    )
    (
        _combine_with_own(acc.pll_sum, acc.pll_weight, own_pll, own_pll_weight),
        fll_disc(
            signal,
            folded.filtered_correlator,
            folded.previous_prompt,
            folded.integration_time,
        ),
        dll_disc(signal, folded.filtered_correlator, code_doppler, sampling_frequency),
    )
end

# Passengers present under vector tracking: fold their discriminators into the
# loop update if this satellite combines, else run the pre-combining path (the
# driver fold above, plus the passengers' prompt / CN0 / bit-buffer advance). The
# flag is read off the per-satellite state, not the shared estimator, and the
# branch is a runtime one rather than a type-level one, for the reasons given at
# the `SatConventionalPLLAndDLL` method.
#
# *Which* loops combine follows `vt_on`, and the second branch is what says so:
# every loop in the scalar fallback, the carrier phase loop alone under vector
# closure, where the other two are the navigation filter's (see
# `VectorPLLAndDLL`). The `VectorClosureAccumulator` marker carries that choice
# into a fold written once — a runtime branch over two concrete accumulator
# types, hence the barrier: the fold proper is its own function, called with
# whichever seed the mode selects, so its body specialises on each.
@inline function _process_signals(
    driver::TrackedSignal,
    passengers::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    sat::TrackedSat,
    state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    ::VectorPLLAndDLL,
)
    state.signal_combining || return _process_signals_separately(
        driver,
        passengers,
        sat,
        state,
        sampling_frequency,
        noise,
        driver_carrier_phase_offset,
    )
    pending = state.pending_discriminators
    if state.vt_on
        _process_signals_combined_vector(
            driver,
            passengers,
            sat,
            state,
            sampling_frequency,
            noise,
            driver_carrier_phase_offset,
            VectorClosureAccumulator(
                pending,
                state.code_discr_acc,
                state.carrier_discr_acc,
            ),
        )
    else
        _process_signals_combined_vector(
            driver,
            passengers,
            sat,
            state,
            sampling_frequency,
            noise,
            driver_carrier_phase_offset,
            pending,
        )
    end
end

# The combining fold: the same interleaved walk as the conventional estimator's
# `_process_signals_combined` — before each driver record closes the loops,
# advance every passenger through the records that completed inside that
# record's window, so a loop update sees exactly the measurements of its own
# interval, and records completing after the chunk's last driver record stay
# pending for the next chunk (see `DiscriminatorAccumulator`).
#
# `acc` arrives seeded and marked by the caller above:
# `_driver_record_discriminators` reads the marker for which loops combine, and
# `_close_vector_loops` closes them.
@inline function _process_signals_combined_vector(
    driver::TrackedSignal,
    passengers::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    sat::TrackedSat,
    state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
    acc::AbstractDiscriminatorAccumulator,
)
    signal = driver.signal
    driver_outputs = driver.correlator_outputs
    driver_noise_density, driver_noise_density_ready = first(noise)
    contexts = _passenger_fold_contexts(
        passengers,
        Base.tail(noise),
        sat.prn,
        sampling_frequency,
        sat.code_doppler,
        driver_carrier_phase_offset,
    )

    tsig = driver
    cursors = map(_ -> 1, passengers)

    carrier_loop_filter = state.carrier_loop_filter
    code_loop_filter = state.code_loop_filter
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    driver_found_before_chunk = has_bit_or_secondary_code_been_found(tsig.bit_buffer)

    @inbounds for k in eachindex(driver_outputs)
        output = driver_outputs[k]
        passengers, cursors, acc =
            _advance_passengers_to(passengers, cursors, contexts, acc, output.sample_index)

        folded = _advance_driver_record(
            tsig,
            output,
            sat.prn,
            sampling_frequency,
            driver_noise_density,
            driver_noise_density_ready,
            driver_found_before_chunk,
            state.carrier_loop_filter_bandwidth,
            state.code_loop_filter_bandwidth,
        )
        tsig = folded.tracked_signal

        pll_discriminator, fll_discriminator, dll_discriminator =
            _driver_record_discriminators(
                acc,
                signal,
                folded,
                output.integrated_samples,
                sat.code_doppler,
                sampling_frequency,
            )

        # The driver's own, slot 1. Its passengers' are folded by the
        # accumulator as the walk reaches them, each into its own slot, so every
        # signal accumulates its own discriminator and never a combination.
        acc = _accumulate_own_discrs(acc, 1, dll_discriminator, fll_discriminator)

        closed = _close_vector_loops(
            state,
            carrier_loop_filter,
            code_loop_filter,
            pll_discriminator,
            fll_discriminator,
            dll_discriminator,
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
        # This window is closed; the next driver record starts a fresh one.
        acc = _zero_window(acc)
    end

    # Drain whatever completed after the last driver record (all of it, when the
    # driver had no record this chunk): the signals must still see those records,
    # and their discriminators stay pending for the next chunk's driver update.
    passengers, _, acc =
        _advance_passengers_to(passengers, cursors, contexts, acc, typemax(Int))
    empty!(driver_outputs)
    _empty_correlator_outputs(passengers)

    # As in the driver-only fold, neither NCO-correction field is written back:
    # both are the navigation filter's to set. The marker is dropped here — what
    # parks on the state is the sums, and the next chunk re-marks them from the
    # `vt_on` in force then.
    new_state = SatVectorPLLAndDLL(
        state;
        carrier_loop_filter,
        code_loop_filter,
        code_discr_acc = _code_discr_accs(acc, state),
        carrier_discr_acc = _carrier_discr_accs(acc, state),
        pending_discriminators = _sums(acc),
    )
    (tsig, passengers, new_state, carrier_doppler, code_doppler)
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
        track_state.doppler_estimator,
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
# cycle is a no-op for satellites already in that state — including for the
# pending combining sums, which a *change* of state must drop: `vt_on` decides
# which loops those sums belong to (see `SatVectorPLLAndDLL`), so passenger
# measurements gathered for three loops must not arrive at a driver record that
# now combines one, or the other way round.
function _set_sat_vt_on(sat, state, vt_on, prns)
    (sat.prn in prns && vt_on != state.vt_on) || return state
    SatVectorPLLAndDLL(
        state;
        vt_on,
        pending_discriminators = zero(DiscriminatorAccumulator),
    )
end

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
    SatVectorPLLAndDLL(state; code_discr_acc = map(_ -> (0, 0.0), state.code_discr_acc)) :
    state

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

!!! warning "A passenger's value is unreferred — read it with its group delay"

    Every signal measures against the satellite's one shared `code_phase`, which
    is the *driver's*, so a passenger's value carries the satellite's
    differential payload group delay on top of the error being measured. Nothing
    in the returned number says so. Fusing it as if it were the driver's puts
    that bias straight into the ranging solution, where — like a wrong group
    delay sign — it is invisible in every tracking metric.

    Reading a passenger is therefore two calls, and the second decides whether
    the first is usable:

    ```julia
    discr = mean_code_discr(sat, i)
    delay = get_differential_group_delay(sat, i)
    isnothing(delay) && continue          # unreferable: withhold, do not assume 0
    referred = discr - delay * code_frequency   # chips
    ```

    `nothing` means the caller has supplied no delay for that signal, and
    withholding is what this package does in the mode where it fuses: the scalar
    path gives such a record zero *code* weight while still taking its carrier
    contribution. Substituting `0.0` is the one unsafe direction — see
    [`set_differential_group_delay!`](@ref).

    Slot 1 needs none of this: the driver's code phase *is* the reference, and
    its delay is `0.0` by definition.

Addressed like every other per-signal accessor
([`estimate_cn0`](@ref), [`get_differential_group_delay`](@ref)): from a
[`TrackState`](@ref) or a [`TrackedSat`](@ref) with a trailing signal selector —
an index or a signal type — which may be omitted only when the satellite tracks
one signal. The `SatVectorPLLAndDLL` form takes an index, since the estimator
state alone cannot resolve a signal type.

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
