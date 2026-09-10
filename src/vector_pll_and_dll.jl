"""
Per-satellite state for the vector PLL and DLL Doppler estimator
([`VectorPLLAndDLL`](@ref)).

On top of the conventional loop-filter state it carries the
vector-tracking (VT) interface to an external navigation filter
(e.g. GNSSReceiver.jl's VDFLL):

  - `code_discr_acc` / `carrier_discr_acc`: `(count, sum)` accumulators of the
    DLL discriminator (chips) and the FLL discriminator (Hz) since the
    navigation filter last read and reset them
    ([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).
    Collected whenever `vt_on`, whatever `discriminator_combining` says. **One pair per
    signal and per loop**, in `sat.signals` order, so a multi-signal satellite
    hands the filter every component's measurement rather than the driver's
    alone — read them with [`mean_code_discr`](@ref) /
    [`mean_carrier_discr`](@ref), which take the same signal selector as every
    other per-signal accessor.

    A count and a sum, with no per-record timing alongside them: every record
    folded in adds one to the count and its discriminator output to the sum,
    and `mean_code_discr` / `mean_carrier_discr` divide the two. Which interval
    the mean covers is therefore the consumer's to know — `count` records of
    the coherent integration time in force — and that holds while the length
    does not change under it. See [`VectorPLLAndDLL`](@ref).

  - `code_freq_update` / `carrier_freq_update`: the NCO corrections the
    navigation filter feeds back ([`set_code_freq_updates!`](@ref),
    [`set_carrier_freq_updates!`](@ref)). While `vt_on`, they replace the
    scalar DLL loop-filter output and the FLL branch of the carrier loop
    filter respectively.

  - `vt_on`: whether the navigation filter controls this satellite's NCOs, and
    therefore whether measurements are collected for it at all. While `false`
    the satellite runs a conventional (scalar) PLL/DLL as a fallback and
    nothing is accumulated. Set by [`enable_vt!`](@ref) /
    [`disable_vt!`](@ref).

`discriminator_combining` and `pending_discriminators` carry multi-signal combining,
which is an independent question from collecting measurements: `vt_on` decides
whether the navigation filter is handed every signal's measurement, while this
flag decides only whether the loops *this* package still closes are fed a
combination. Which loops those are does follow `vt_on` in turn: every loop while
this satellite runs its own scalar fallback, the carrier phase loop alone once
the navigation filter has the other two (see [`VectorPLLAndDLL`](@ref)). The flag is seeded from the shared
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
    discriminator_combining::Bool = false
    pending_discriminators::DiscriminatorAccumulator = DiscriminatorAccumulator()
    # No defaults, deliberately: a default would have to pick an `N`, and a
    # one-slot default silently truncates a multi-signal satellite — slot 2's
    # measurements would be dropped by `_accumulate_one_discr` and
    # `mean_code_discr(state, 2)` would throw. `N` is the satellite's signal
    # count, so it comes from the satellite: use the `SatVectorPLLAndDLL(sat, …)`
    # constructor below, or pass `_zero_code_discr_acc(sat)` /
    # `_zero_carrier_discr_acc(sat)` here.
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
    discriminator_combining::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    SatVectorPLLAndDLL(;
        init_carrier_doppler = sat.carrier_doppler,
        init_code_doppler = sat.code_doppler,
        carrier_loop_filter,
        code_loop_filter,
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        discriminator_combining,
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
    discriminator_combining::Maybe{Bool} = nothing,
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
        isnothing(discriminator_combining) ?
        sat_vector_pll_and_dll.discriminator_combining : discriminator_combining,
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

Collecting those measurements follows `vt_on` alone: a satellite in the vector
loop hands the filter every signal's measurement whether or not it combines
anything for its own loops.

!!! note "A `(count, sum)` pair maps to an interval only at one length"

    `count × integration_time` is the interval the mean covers only while every
    record in it ran at the same coherent integration time. Calling
    [`set_preferred_num_code_blocks_to_integrate!`](@ref) while `vt_on` is set
    leaves the mean spanning two lengths with nothing in the pair to say so, so
    reset the accumulators across such a change or make it with the satellite
    out of the vector loop. The vector-tracking manual works this through.

With `discriminator_combining = true` a multi-signal satellite folds its signals'
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
    discriminator_combining::Bool
end

function VectorPLLAndDLL(
    ::Type{CA} = ThirdOrderAssistedBilinearLF,
    ::Type{CO} = SecondOrderBilinearLF;
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    discriminator_combining::Bool = false,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        discriminator_combining,
    )
end

# Kwarg-update constructor for tweaking bandwidths / combining in place.
function VectorPLLAndDLL(
    pll_and_dll::VectorPLLAndDLL{CA,CO};
    carrier_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    code_loop_filter_bandwidth::Maybe{typeof(1.0Hz)} = nothing,
    discriminator_combining::Maybe{Bool} = nothing,
) where {CA<:AbstractLoopFilter,CO<:AbstractLoopFilter}
    VectorPLLAndDLL{CA,CO}(
        isnothing(carrier_loop_filter_bandwidth) ?
        pll_and_dll.carrier_loop_filter_bandwidth : carrier_loop_filter_bandwidth,
        isnothing(code_loop_filter_bandwidth) ? pll_and_dll.code_loop_filter_bandwidth :
        code_loop_filter_bandwidth,
        isnothing(discriminator_combining) ? pll_and_dll.discriminator_combining :
        discriminator_combining,
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
        discriminator_combining = estimator.discriminator_combining,
        code_discr_acc = _zero_code_discr_acc(sat),
        carrier_discr_acc = _zero_carrier_discr_acc(sat),
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Per-sat bandwidth overrides, the
# `discriminator_combining` flag and the `vt_on` flag survive the reset; the pending
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
        discriminator_combining = state.discriminator_combining,
        code_discr_acc = map(_ -> (0, 0.0), state.code_discr_acc),
        carrier_discr_acc = map(_ -> (0, 0.0Hz), state.carrier_discr_acc),
        vt_on = state.vt_on,
    )
end

# What this estimator does with one record's three discriminators. The shared
# traversal takes it as its loop-closure callback (`_loop_closure`), which is
# the only thing that differs between this estimator's fold and the
# conventional one's.
#
# While `vt_on`, the navigation filter's NCO corrections drive the loops —
# `code_freq_update` directly (the code loop filter bypassed, holding its state)
# and `carrier_freq_update` through the FLL branch of the carrier loop filter,
# the PLL branch still running on this satellite's own phase discriminator. With
# `vt_on` unset it is a conventional FLL-assisted scalar PLL/DLL.
#
# The discriminators arrive already formed because only `discriminators.pll` is
# ever a combination of several signals.
@inline function _close_vector_loops(
    state::SatVectorPLLAndDLL,
    carrier_loop_filter,
    code_loop_filter,
    discriminators,
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
            discriminators.dll,
            integration_time,
            code_bandwidth,
        )
        fll_input = discriminators.fll
    end
    # The shared `_filter_carrier_loop`: only the FLL-assisted filter has an
    # input path for `fll_input`, and any other runs on the PLL discriminator
    # alone — so the vector carrier closure degrades to PLL-only there, exactly
    # as the scalar one does.
    carrier_freq_update, carrier_loop_filter = _filter_carrier_loop(
        carrier_loop_filter,
        discriminators.pll,
        fll_input,
        integration_time,
        carrier_bandwidth,
    )
    (; carrier_freq_update, code_freq_update, carrier_loop_filter, code_loop_filter)
end

@inline _loop_closure(::SatVectorPLLAndDLL) = _close_vector_loops

# Fold one signal's discriminator into its own `(count, sum)`. `Val(N)` keeps
# the rebuilt tuple's length a compile-time constant, so this is
# allocation-free and the slot write costs a predicted branch per signal.
@inline _accumulate_one_discr(accs::NTuple{N,Any}, slot::Integer, value) where {N} =
    ntuple(k -> k == slot ? (accs[k][1] + 1, accs[k][2] + value) : accs[k], Val(N))

# ---------------------------------------------------------------------------
# Multi-signal discriminator combining
# ---------------------------------------------------------------------------

"""
The per-signal navigation measurements of a satellite whose `vt_on` flag is set:
one `(count, sum)` pair per signal for the DLL and one for the FLL, in
`sat.signals` order, growing until the navigation filter reads and resets them
([`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)).

This is the accumulator of a satellite that collects measurements and combines
nothing; [`VectorLoopAccumulator`](@ref) wraps it when it does both (see
[`VectorPLLAndDLL`](@ref) for why the two are separate questions).

Every signal accumulates its **own** discriminators here and never a
combination across signals: fusing them is the navigation filter's job.
"""
struct PerSignalAccumulator{N} <: AbstractDiscriminatorAccumulator
    code_discr_acc::NTuple{N,Tuple{Int,Float64}}
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(1.0Hz)}}
end

@inline PerSignalAccumulator(state::SatVectorPLLAndDLL) =
    PerSignalAccumulator(state.code_discr_acc, state.carrier_discr_acc)

# Nothing combines here, so there is no window to close and no pending sum to
# park: a driver record closes a *loop* window, while these answer to the
# navigation filter's reading cadence.
@inline _reset_combining_window(acc::PerSignalAccumulator) = acc

"""
The accumulator of a satellite that both combines and measures: the carrier
phase loop combines across signals, and every signal's code and carrier
frequency discriminators go to the navigation filter as that signal's own — see
[`VectorPLLAndDLL`](@ref) for why the split falls there.

A composition of the two halves rather than one flat struct, so the measurement
half is the very same [`PerSignalAccumulator`](@ref) a non-combining
vector-closed satellite carries, folded by the same code.

It is a marker as much as a container: which loops a satellite combines follows
`vt_on` and therefore changes at run time, while the state field the sums park
in cannot change type. The traversal is written once against
`AbstractDiscriminatorAccumulator` and specialises on this, so a vector-closed
satellite forms no weight and no noise gain for the two loops it does not
combine.
"""
struct VectorLoopAccumulator{N} <: AbstractDiscriminatorAccumulator
    combined::DiscriminatorAccumulator
    measurements::PerSignalAccumulator{N}
end

@inline VectorLoopAccumulator(state::SatVectorPLLAndDLL) =
    VectorLoopAccumulator(state.pending_discriminators, PerSignalAccumulator(state))

# A driver record closes the combining window, not the per-signal measurements.
@inline _reset_combining_window(acc::VectorLoopAccumulator) =
    VectorLoopAccumulator(DiscriminatorAccumulator(), acc.measurements)

# Where a chunk's accumulator ends up on the per-satellite state. One method per
# accumulator type, so each mode writes exactly the fields it owns and leaves
# every other one as it found it: a mode that does not combine must not clear
# the pending combining sums, and one that does not measure must not clear the
# navigation filter's accumulators.
@inline _park_accumulator(
    state::SatVectorPLLAndDLL,
    acc::VectorLoopAccumulator,
    carrier_loop_filter,
    code_loop_filter,
) = SatVectorPLLAndDLL(
    state;
    carrier_loop_filter,
    code_loop_filter,
    code_discr_acc = acc.measurements.code_discr_acc,
    carrier_discr_acc = acc.measurements.carrier_discr_acc,
    pending_discriminators = acc.combined,
)

@inline _park_accumulator(
    state::SatVectorPLLAndDLL,
    acc::PerSignalAccumulator,
    carrier_loop_filter,
    code_loop_filter,
) = SatVectorPLLAndDLL(
    state;
    carrier_loop_filter,
    code_loop_filter,
    code_discr_acc = acc.code_discr_acc,
    carrier_discr_acc = acc.carrier_discr_acc,
)

@inline _park_accumulator(
    state::SatVectorPLLAndDLL,
    acc::DiscriminatorAccumulator,
    carrier_loop_filter,
    code_loop_filter,
) = SatVectorPLLAndDLL(
    state;
    carrier_loop_filter,
    code_loop_filter,
    pending_discriminators = acc,
)

# Neither combining nor measuring: only the two threaded loop filters are this
# chunk's output. The pending combining sums stay as they were — a satellite
# that stopped combining mid-flight must not have them cleared behind its back.
@inline _park_accumulator(
    state::SatVectorPLLAndDLL,
    ::Nothing,
    carrier_loop_filter,
    code_loop_filter,
) = SatVectorPLLAndDLL(state; carrier_loop_filter, code_loop_filter)

# Fold one signal's own two navigation measurements into its slot. The
# `DiscriminatorAccumulator` no-op is in conventional_pll_and_dll.jl: that is
# the scalar fallback, where the navigation filter reads nothing.
@inline _accumulate_signal_measurements(
    acc::PerSignalAccumulator{N},
    slot::Integer,
    dll_discriminator,
    fll_discriminator,
) where {N} = PerSignalAccumulator{N}(
    _accumulate_one_discr(acc.code_discr_acc, slot, dll_discriminator),
    _accumulate_one_discr(acc.carrier_discr_acc, slot, fll_discriminator),
)

@inline _accumulate_signal_measurements(
    acc::VectorLoopAccumulator,
    slot::Integer,
    dll_discriminator,
    fll_discriminator,
) = VectorLoopAccumulator(
    acc.combined,
    _accumulate_signal_measurements(
        acc.measurements,
        slot,
        dll_discriminator,
        fll_discriminator,
    ),
)

# One passenger record's two navigation measurements into its own slot, shared
# by both vector accumulators — collecting them follows `vt_on`, which both have
# set, and not what else the accumulator does with the record. Nothing reads the
# differential group delay: referring a code measurement to one ranging datum is
# the navigation filter's job here, and it is the party that knows which datum
# it ranges on.
@inline function _measure_passenger_record(
    acc::AbstractDiscriminatorAccumulator,
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    filtered_correlator,
    previous_prompt,
    ctx::PassengerFoldContext,
)
    signal = tracked_signal.signal
    integration_time = output.integrated_samples / ctx.sampling_frequency
    _accumulate_signal_measurements(
        acc,
        ctx.signal_index,
        dll_disc(signal, filtered_correlator, ctx.code_doppler, ctx.sampling_frequency),
        fll_disc(signal, filtered_correlator, previous_prompt, integration_time),
    )
end

# A vector-closed satellite that does not combine: the record's two navigation
# measurements are all it contributes, and the loops stay the driver's own. The
# pre-sync gate excludes an invalidated record for the same reason it does when
# combining — a replica the sync had not yet corrected carries no usable phase or
# delay information, whether the value is about to be averaged or shipped. It
# rarely fires here (it needs a satellite admitted to the vector loop with a
# secondary-coded signal still unsynced), but it belongs to the record rather
# than to the mode.
@inline function _accumulate_passenger_discriminators(
    acc::PerSignalAccumulator,
    tracked_signal::TrackedSignal,
    output::CorrelatorOutput,
    filtered_correlator,
    previous_prompt,
    ctx::PassengerFoldContext,
    correlated_pre_sync::Bool,
)
    _passenger_record_invalidated_by_sync(tracked_signal.signal, correlated_pre_sync) &&
        return acc
    _measure_passenger_record(
        acc,
        tracked_signal,
        output,
        filtered_correlator,
        previous_prompt,
        ctx,
    )
end

# One passenger record under vector closure with combining on: its carrier phase
# measurement into the combination, its other two into its own accumulation.
@inline function _accumulate_passenger_discriminators(
    acc::VectorLoopAccumulator,
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
        _weighted_pll_contribution(
            signal,
            filtered_correlator,
            output.integrated_samples,
            ctx.derotation,
        )...,
    )
    _measure_passenger_record(
        VectorLoopAccumulator(combined, acc.measurements),
        tracked_signal,
        output,
        filtered_correlator,
        previous_prompt,
        ctx,
    )
end

# Carrier phase combined, the other two this signal's own — and formed by the
# plain discriminators rather than by `_weighted_record_contribution`, so no
# weight, no noise gain and no group delay is computed for a value nothing will
# weigh. The conventional estimator's method, which combines all three, is in
# conventional_pll_and_dll.jl.
#
# Both methods here form the FLL discriminator whatever the carrier loop filter
# is: this accumulator is collecting it for the navigation filter, so the filter
# is not its only reader. `carrier_loop_filter` is therefore unread — only the
# `nothing` accumulator gates on it (see `_driver_fll_discriminator`).
@inline function _combined_loop_discriminators(
    acc::VectorLoopAccumulator,
    signal::AbstractGNSSSignal,
    ::AbstractLoopFilter,
    folded,
    integrated_samples::Integer,
    code_doppler,
    sampling_frequency,
)
    combined = acc.combined
    own_pll, own_pll_weight = _weighted_pll_contribution(
        signal,
        folded.filtered_correlator,
        integrated_samples,
        one(ComplexF64),
    )
    (;
        pll = _combine_with_own(
            combined.pll_sum,
            combined.pll_weight,
            own_pll,
            own_pll_weight,
        ),
        fll = fll_disc(
            signal,
            folded.filtered_correlator,
            folded.previous_prompt,
            folded.integration_time,
        ),
        dll = dll_disc(
            signal,
            folded.filtered_correlator,
            code_doppler,
            sampling_frequency,
        ),
    )
end

# Nothing combines: every loop closes on the driver's own discriminators, so a
# satellite that collects measurements without combining tracks bit-identically
# to one that collects nothing at all.
@inline _combined_loop_discriminators(
    ::PerSignalAccumulator,
    signal::AbstractGNSSSignal,
    ::AbstractLoopFilter,
    folded,
    ::Integer,
    code_doppler,
    sampling_frequency,
) = (;
    pll = pll_disc(signal, folded.filtered_correlator),
    fll = fll_disc(
        signal,
        folded.filtered_correlator,
        folded.previous_prompt,
        folded.integration_time,
    ),
    dll = dll_disc(signal, folded.filtered_correlator, code_doppler, sampling_frequency),
)

# Under vector tracking, `vt_on` and `discriminator_combining` are independent (see
# `VectorPLLAndDLL`), and their four combinations are four seeds of one
# traversal: measurements only (`PerSignalAccumulator`), combining only
# (`DiscriminatorAccumulator`), both (`VectorLoopAccumulator`), or neither —
# `nothing`, which folds the whole combining and measurement apparatus away and
# leaves the plain scalar fallback.
#
# Both flags are read off the per-satellite state, not the shared estimator, and
# the combining branch is a runtime one for the reason given at the
# `SatConventionalPLLAndDLL` method. The traversal specialises on each seed
# behind it.

# No passengers: there is nothing to combine, so only `vt_on` forks — the same
# two seeds this estimator took before passengers existed.
@inline _seed_and_fold(
    driver::TrackedSignal,
    passengers::Tuple{},
    sat::TrackedSat,
    state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
) = _fold_satellite_chunk(
    driver,
    passengers,
    sat,
    state,
    sampling_frequency,
    noise,
    driver_carrier_phase_offset,
    state.vt_on ? PerSignalAccumulator(state) : nothing,
)

@inline function _seed_and_fold(
    driver::TrackedSignal,
    passengers::Tuple{TrackedSignal,Vararg{TrackedSignal}},
    sat::TrackedSat,
    state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase_offset::Real,
)
    fold(acc) = _fold_satellite_chunk(
        driver,
        passengers,
        sat,
        state,
        sampling_frequency,
        noise,
        driver_carrier_phase_offset,
        acc,
    )
    if state.discriminator_combining
        state.vt_on ? fold(VectorLoopAccumulator(state)) :
        fold(state.pending_discriminators)
    else
        state.vt_on ? fold(PerSignalAccumulator(state)) : fold(nothing)
    end
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
    SatVectorPLLAndDLL(state; vt_on, pending_discriminators = DiscriminatorAccumulator())
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
    differential payload group delay on top of the error being measured, and
    nothing in the returned number says so. Subtract it before fusing:

    ```julia
    discr = mean_code_discr(sat, i)
    delay = get_driver_relative_group_delay(sat, i)
    isnothing(delay) && continue          # unreferable: withhold, do not assume 0
    referred = discr - delay * code_frequency   # chips
    ```

    [`get_driver_relative_group_delay`](@ref) is the difference already formed, so
    there is no special case for slot 1 — it answers `0.0s` and the subtraction is
    a no-op there. `nothing` means this signal's value or the satellite's datum is
    missing; withhold the measurement rather than substituting `0.0`, which is the
    one unsafe direction (see [`set_differential_group_delay!`](@ref)).

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
