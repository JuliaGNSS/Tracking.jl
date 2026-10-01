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
    **One pair per signal**, in `sat.signals` order, each holding that
    signal's own raw discriminator, so the navigation filter can fuse the
    signals itself. Read them with [`mean_code_discr`](@ref) /
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
  - `discriminator_combining`: whether the passenger signals aid the loops
    this package closes — see [`VectorPLLAndDLL`](@ref).

The accumulators have no default: their length is the satellite's signal
count. Build the state from the satellite (`SatVectorPLLAndDLL(sat, …)`) or
through [`init_estimator_state`](@ref).
"""
@kwdef struct SatVectorPLLAndDLL{CA<:AbstractLoopFilter,CO<:AbstractLoopFilter,N}
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_loop_filter::CA = ThirdOrderAssistedBilinearLF()
    code_loop_filter::CO = SecondOrderBilinearLF()
    carrier_loop_filter_bandwidth::typeof(1.0Hz) = 18.0Hz
    code_loop_filter_bandwidth::typeof(1.0Hz) = 1.0Hz
    code_discr_acc::NTuple{N,Tuple{Int,Float64}}
    code_freq_update::typeof(0.0Hz) = 0.0Hz
    carrier_discr_acc::NTuple{N,Tuple{Int,typeof(0.0Hz)}}
    carrier_freq_update::typeof(0.0Hz) = 0.0Hz
    vt_on::Bool = false
    discriminator_combining::Bool = false
end

# One zeroed accumulator per signal. Mapping over the signal tuple keeps the
# count a compile-time constant, so every satellite of a group has the same
# state type.
@inline _zero_code_discr_acc(sat::TrackedSat) = map(_ -> (0, 0.0), sat.signals)
@inline _zero_carrier_discr_acc(sat::TrackedSat) = map(_ -> (0, 0.0Hz), sat.signals)

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
        code_discr_acc = _zero_code_discr_acc(sat),
        carrier_discr_acc = _zero_carrier_discr_acc(sat),
        discriminator_combining,
    )
end

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
    discriminator_combining::Maybe{Bool} = nothing,
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
        isnothing(discriminator_combining) ?
        sat_vector_pll_and_dll.discriminator_combining : discriminator_combining,
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

  - This estimator accumulates each signal's DLL / FLL discriminator
    outputs for the navigation filter to consume (and reset via
    [`reset_code_discr_acc!`](@ref) / [`reset_carrier_discr_acc!`](@ref)) —
    one raw measurement per signal, never a combination: fusing them is the
    navigation filter's job.
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

`discriminator_combining = true` combines a multi-signal satellite's passengers
into the loops this package still closes, by the rules of
[`ConventionalPLLAndDLL`](@ref) (see [Discriminator combining](@ref)). Which
loops those are follows `vt_on`:

  - `vt_on = false`: all three, exactly as [`ConventionalPLLAndDLL`](@ref)
    combines them, so the scalar fallback is the conventional estimator.
  - `vt_on = true`: the carrier phase loop only. The code and carrier frequency
    loops are the navigation filter's, which receives every signal's own raw
    measurement instead — with no group-delay difference applied, since the
    filter applies its own inter-signal biases.

The per-signal accumulators are filled whatever the flag says.
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

# Kwarg-update constructor for tweaking the configuration in place.
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
    SatVectorPLLAndDLL(
        sat,
        carrier_loop_filter,
        code_loop_filter;
        carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth,
        estimator.discriminator_combining,
    )
end

# Re-seed hook used by `reset_loop_filters!`: zero the loop-filter
# integrators, the discriminator accumulators, and the NCO corrections, and
# re-seed the init Dopplers from the sat's current (converged) Dopplers.
# The NCO corrections must be zeroed together with the init Dopplers: the
# current Dopplers already contain the last correction, so keeping it would
# apply it twice after the re-seed. Per-sat bandwidth overrides, the `vt_on`
# flag and the combining flag survive the reset.
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
        discriminator_combining = state.discriminator_combining,
    )
end

# The vector estimator's own fold over a satellite's chunk, replacing the
# conventional one by dispatch on the per-sat state. It runs the same
# combining as the conventional fold — the same coincidence rule, contexts,
# weights and per-record passenger step, all shared — and differs in two
# places:
#
#   - Loop closure. While `vt_on`, the navigation filter's NCO corrections drive
#     the loops: `code_freq_update` directly (code loop filter bypassed) and
#     `carrier_freq_update` through the FLL branch of the carrier loop filter,
#     whose PLL branch still runs on the (combined) phase discriminator. With
#     `vt_on` unset it is the conventional scalar PLL/DLL, so this fold and the
#     conventional one must close the loops identically; a test pins that.
#   - Measurements. While `vt_on`, every record of every signal adds its own raw
#     DLL / FLL discriminator to that signal's accumulator, combined or not. So
#     the passengers that are not combined are folded here too, rather than
#     left to the shared passenger fold, which forms no discriminators.
@inline function _fold_driver(
    tracked_signal::TrackedSignal,
    passengers::Tuple,
    sat::TrackedSat,
    pll_and_dll_state::SatVectorPLLAndDLL,
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase::Real,
)
    noise_density, noise_density_ready = first(noise)
    signal = tracked_signal.signal
    vt_on = pll_and_dll_state.vt_on
    # Captures `tracked_signal`, never `ts`: a captured variable that is
    # reassigned later is boxed, and the fold would allocate per call.
    contexts = map(passengers, Base.tail(noise)) do passenger, passenger_noise
        _passenger_context(
            tracked_signal,
            passenger,
            passenger_noise,
            pll_and_dll_state.discriminator_combining,
            sat.code_doppler,
            driver_carrier_phase,
        )
    end
    ts = tracked_signal
    outputs = ts.correlator_outputs
    carrier_loop_filter = pll_and_dll_state.carrier_loop_filter
    code_loop_filter = pll_and_dll_state.code_loop_filter
    code_discr_acc = pll_and_dll_state.code_discr_acc
    carrier_discr_acc = pll_and_dll_state.carrier_discr_acc
    carrier_doppler = sat.carrier_doppler
    code_doppler = sat.code_doppler
    found_before_fold = has_bit_or_secondary_code_been_found(ts.bit_buffer)
    @inbounds for k in eachindex(outputs)
        output = outputs[k]
        # Every passenger's FLL reading is formed, whatever the carrier loop
        # filter: under `vt_on` it is a measurement for the navigation filter.
        passengers, contributions = _apply_coincident_passenger_records(
            passengers,
            contexts,
            k,
            sat.prn,
            sampling_frequency,
            sat.code_doppler,
            driver_carrier_phase,
            true,
        )
        passenger_sums = _passenger_sums(contributions, contexts)
        # FLL needs the previous record's filtered prompt; the first record of
        # the chunk chains from the sat's carried-over
        # `last_fully_integrated_filtered_prompt` (the previous chunk's last).
        # Read it off `ts` BEFORE the advance overwrites it.
        previous_prompt = get_last_fully_integrated_filtered_prompt(ts)
        # Per-record integration time — the block time, NOT the chunk time.
        integration_time = output.integrated_samples / sampling_frequency
        # A record that follows a sync detected earlier in THIS fold was
        # correlated with pre-sync replicas — its blocks still count towards the
        # bit, only its prompt may have to be dropped (see
        # `_apply_correlator_output`).
        synced_earlier_in_fold =
            !found_before_fold && has_bit_or_secondary_code_been_found(ts.bit_buffer)
        # The driver de-rotates against itself (offset 0), so the derotation is a
        # no-op for it; passed for symmetry with the passenger path.
        ts, filtered_correlator, integrated_code_blocks = _apply_correlator_output(
            ts,
            output,
            sat.prn,
            sampling_frequency,
            noise_density,
            noise_density_ready,
            driver_carrier_phase;
            correlated_pre_sync = synced_earlier_in_fold,
        )

        # The configured CARRIER bandwidth is referenced to a
        # one-primary-code-period integration. Coherently integrating N periods
        # grows the loop update interval by that factor, so scale the effective
        # bandwidth by 1/N to hold the loop's BL·Δt stability product at its
        # single-period value. N (`integrated_code_blocks`) is the number of
        # blocks this record ACTUALLY covered — recovered from its sample count —
        # not the intended integration length: the bandwidth pairs with the
        # record's true `integration_time`, so it only switches once the
        # integrations really lengthen. Scaling by the intended length instead
        # would under-gain the loop for single-block records folded after a
        # mid-fold sync detection (correlated pre-sync, but the live bit buffer
        # already reports the post-sync length) and for the first post-sync
        # integration, which is truncated to land on the data-bit boundary. For
        # the N=1 path this divides by 1 and is bit-identical to before.
        carrier_bandwidth =
            pll_and_dll_state.carrier_loop_filter_bandwidth / integrated_code_blocks
        # The DLL's is an absolute bandwidth, so it is capped by its own
        # stability product against this record's integration time instead of
        # scaled by N — see `effective_code_loop_filter_bandwidth`.
        code_bandwidth = effective_code_loop_filter_bandwidth(
            pll_and_dll_state.code_loop_filter_bandwidth,
            integration_time,
        )

        # The driver's own weight is its power share, like every passenger's;
        # its FLL weight is zero where there is no previous prompt yet.
        weight = get_relative_power(signal)
        pll_discriminator = pll_disc(signal, filtered_correlator)
        # Formed whatever the carrier loop filter: under `vt_on` it is a
        # measurement for the navigation filter.
        fll_discriminator =
            fll_disc(signal, filtered_correlator, previous_prompt, integration_time)
        # `dll_disc` is fed the chunk-fixed `sat.code_doppler` — the code Doppler
        # that actually generated this chunk's replicas — for every record;
        # only the loop-filter *state* threads across records.
        dll_discriminator =
            dll_disc(signal, filtered_correlator, sat.code_doppler, sampling_frequency)
        pll = _combined_discriminator(
            pll_discriminator,
            weight,
            passenger_sums.pll,
            passenger_sums.pll_weight,
        )

        if vt_on
            code_discr_acc = (
                _add_discr(first(code_discr_acc), dll_discriminator),
                _add_measured(Base.tail(code_discr_acc), contributions, :dll)...,
            )
            carrier_discr_acc = (
                _add_discr(first(carrier_discr_acc), fll_discriminator),
                _add_measured(Base.tail(carrier_discr_acc), contributions, :fll)...,
            )
            # The code loop and the carrier loop's FLL branch follow the
            # navigation filter's corrections.
            code_freq_update = pll_and_dll_state.code_freq_update
            fll = pll_and_dll_state.carrier_freq_update
        else
            fll = _combined_discriminator(
                fll_discriminator,
                iszero(previous_prompt) ? 0.0 : weight,
                passenger_sums.fll,
                passenger_sums.fll_weight,
            )
            dll = _combined_discriminator(
                dll_discriminator,
                weight,
                passenger_sums.dll,
                passenger_sums.dll_weight,
            )
            code_freq_update, code_loop_filter =
                filter_loop(code_loop_filter, dll, integration_time, code_bandwidth)
        end
        carrier_freq_update, carrier_loop_filter = filter_loop(
            carrier_loop_filter,
            _carrier_loop_input(carrier_loop_filter, pll, fll),
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
    end
    empty!(outputs)
    folded = map(
        passengers,
        contexts,
        Base.tail(code_discr_acc),
        Base.tail(carrier_discr_acc),
    ) do passenger, context, code_acc, carrier_acc
        _fold_uncombined_passenger(
            passenger,
            context,
            code_acc,
            carrier_acc,
            vt_on,
            sat,
            sampling_frequency,
            driver_carrier_phase,
        )
    end
    passengers = map(first, folded)
    code_discr_acc = (first(code_discr_acc), map(f -> f[2], folded)...)
    carrier_discr_acc = (first(carrier_discr_acc), map(f -> f[3], folded)...)
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
    return ts, new_doppler_estimator_state, carrier_doppler, code_doppler, passengers
end

@inline _add_discr(acc::Tuple{Int,T}, value) where {T} = (acc[1] + 1, acc[2] + value)

# Add each combined passenger's raw reading of one loop to its own accumulator;
# a passenger whose record was not combined or was excluded adds nothing here.
@inline _add_measured(accs::Tuple, contributions::Tuple, loop::Symbol) =
    map(accs, contributions) do acc, contribution
        contribution.measured ? _add_discr(acc, getfield(contribution, loop)) : acc
    end

# A passenger whose records were not combined this chunk: apply them in order,
# and while `vt_on` add each one's raw DLL / FLL reading to its accumulators.
# A combined passenger's records were applied in the driver loop already.
@inline function _fold_uncombined_passenger(
    tracked_signal::TrackedSignal,
    context::NamedTuple,
    code_acc,
    carrier_acc,
    vt_on::Bool,
    sat::TrackedSat,
    sampling_frequency,
    driver_carrier_phase::Real,
)
    outputs = tracked_signal.correlator_outputs
    ts = tracked_signal
    if !context.combined
        @inbounds for k in eachindex(outputs)
            output = outputs[k]
            ts, filtered_correlator, previous_prompt, excluded = _apply_passenger_record(
                ts,
                output,
                context.found_before_fold,
                sat.prn,
                sampling_frequency,
                context.noise_density,
                context.noise_density_ready,
                driver_carrier_phase,
            )
            if vt_on && !excluded
                signal = ts.signal
                code_acc = _add_discr(
                    code_acc,
                    dll_disc(
                        signal,
                        filtered_correlator,
                        sat.code_doppler,
                        sampling_frequency,
                    ),
                )
                carrier_acc = _add_discr(
                    carrier_acc,
                    fll_disc(
                        signal,
                        filtered_correlator,
                        previous_prompt,
                        output.integrated_samples / sampling_frequency,
                    ),
                )
            end
        end
    end
    empty!(outputs)
    ts, code_acc, carrier_acc
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
[`reset_code_discr_acc!`](@ref), in chips, or `nothing` if nothing has been
accumulated yet (the `(count, sum)` accumulator's `count` is 0). This is the
single place the accumulator's averaging convention lives — read it here
rather than dividing `code_discr_acc` by hand.

Every signal of a satellite accumulates its **own** raw discriminator, on its
own code phase: no [`get_group_delay`](@ref) difference is applied, and nothing
is combined across signals. Fusing the signals, and applying their
inter-signal biases, is the navigation filter's job. Every record counts except
a passenger's correlated with a replica that a same-chunk secondary-code sync
had not yet corrected, which skips both accumulators.

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

# A single-signal satellite needs no selector; a multi-signal one is refused as
# an unqualified per-signal read is, rather than silently answering with the
# driver's.
@inline _sole_signal_index(::NTuple{1,Any}) = 1
@inline _sole_signal_index(::Tuple) = _throw_needs_signal_selector()

mean_code_discr(state::SatVectorPLLAndDLL) =
    mean_code_discr(state, _sole_signal_index(state.code_discr_acc))
mean_carrier_discr(state::SatVectorPLLAndDLL) =
    mean_carrier_discr(state, _sole_signal_index(state.carrier_discr_acc))

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
            sig::Union{Integer,Type{<:AbstractGNSSSignal}},
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
