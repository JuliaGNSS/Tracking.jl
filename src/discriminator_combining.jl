# Discriminator combining: a multi-signal satellite's passengers (`signals[2:end]`)
# aid the loops its driver (`signals[1]`) closes. Used by `ConventionalPLLAndDLL`
# with `combine_discriminators`; see "Discriminator combining" in the docs for
# the rules. Also home of the group-delay API.

# The stored form of an unknown group delay.
const _UNKNOWN_GROUP_DELAY = NaN * 1.0s

# A per-sat state that is not Tracking's holds no group delays: the group-delay
# API refuses it by name.
_group_delays(state) = throw(
    ArgumentError(
        "$(nameof(typeof(state))) holds no group delays: only Tracking's own " *
        "estimators store them",
    ),
)

# The shared driver fold of a combining conventional state, with the coincident
# passengers taking part through its passenger hooks (below); the others are
# left to the shared passenger fold.
@inline function _fold_driver(
    tracked_signal::TrackedSignal,
    passengers::Tuple,
    sat::TrackedSat,
    state::SatConventionalPLLAndDLL{<:Any,<:Any,<:Any,true},
    sampling_frequency,
    noise::Tuple,
    driver_carrier_phase::Real,
)
    # Nothing to fold: the shared fold's own early return, without forming the
    # contexts. No passenger coincides with a driver that has no records.
    if isempty(tracked_signal.correlator_outputs)
        return tracked_signal, state, sat.carrier_doppler, sat.code_doppler, passengers
    end
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
end

# The passengers' contexts for this chunk (`_passenger_context`), for any per-sat
# state that holds group delays (`_group_delays`).
@inline function _passenger_contexts(
    tracked_signal::TrackedSignal,
    passengers::Tuple,
    sat::TrackedSat,
    state,
    noise::Tuple,
    driver_carrier_phase::Real,
)
    group_delays = _group_delays(state)
    driver_delay = first(group_delays)
    # Captures `tracked_signal`, never a variable reassigned later: such a
    # capture is boxed, and the fold would allocate per call.
    contexts = map(
        passengers,
        Base.tail(noise),
        Base.tail(group_delays),
    ) do passenger, passenger_noise, delay
        _passenger_context(
            tracked_signal,
            passenger,
            passenger_noise,
            driver_delay,
            delay,
            sat.code_doppler,
            driver_carrier_phase,
        )
    end
    PassengerContexts(contexts)
end

# Whether a passenger's records measure the same intervals as the driver's: as
# many, each ending on the same sample after the same number of samples.
@inline _records_coincide(driver_outputs::Vector, passenger_outputs::Vector) =
    length(driver_outputs) == length(passenger_outputs) &&
    all(zip(driver_outputs, passenger_outputs)) do (driver_output, passenger_output)
        driver_output.sample_index == passenger_output.sample_index &&
        driver_output.integrated_samples == passenger_output.integrated_samples
    end

# Everything one passenger's records need for this chunk, formed once:
#
#   - `coincident`: whether its records coincide with the driver's this chunk;
#   - `found_before_fold`: its sync state at the chunk boundary, for the
#     pre-sync rule of `_apply_correlator_output`;
#   - `derotation`: the rotation onto the driver's carrier phase frame;
#   - `dll_offset`: how far, in chips, this signal's received code sits ahead of
#     the driver's against the shared replica, from the two group delays; its DLL
#     reads the driver's error plus this offset. `NaN` where either delay is
#     unknown, which keeps the passenger out of the code loop only;
#   - its own noise density, never the driver's.
@inline function _passenger_context(
    driver::TrackedSignal,
    passenger::TrackedSignal,
    (noise_density, noise_density_ready)::Tuple{Any,Bool},
    driver_delay::typeof(1.0s),
    delay::typeof(1.0s),
    code_doppler,
    driver_carrier_phase::Real,
)
    # `NaN` for an unknown delay on either end carries through to the offset.
    dll_offset = uconvert(
        NoUnits,
        (driver_delay - delay) * (get_code_frequency(passenger.signal) + code_doppler),
    )
    (;
        coincident = _records_coincide(
            driver.correlator_outputs,
            passenger.correlator_outputs,
        ),
        found_before_fold = has_bit_or_secondary_code_been_found(passenger.bit_buffer),
        derotation = _carrier_phase_derotation(driver_carrier_phase, passenger.signal),
        dll_offset,
        noise_density,
        noise_density_ready,
    )
end

# The passengers' contexts for this chunk, in passenger order: what the shared
# driver fold's passenger hooks dispatch on.
struct PassengerContexts{C<:Tuple}
    contexts::C
end

# What one passenger record contributes to its driver record's loop update: its
# raw readings and their weights. `dll` is the raw reading against the shared
# replica; the referral by `dll_offset` happens where the sums are formed.
# `coincident` is false for a passenger that did not take part (all weights
# zero); an estimator that also accumulates measurements reads the raw values.
struct PassengerContribution
    coincident::Bool
    pll_weight::Float64
    pll::Float64
    fll_weight::Float64
    fll::typeof(0.0Hz)
    dll_weight::Float64
    dll::Float64
    dll_offset::Float64
end

const _NO_CONTRIBUTION = PassengerContribution(false, 0.0, 0.0, 0.0, 0.0Hz, 0.0, 0.0, NaN)

# Per-loop sums over the passengers of one driver record: `Σ wᵢ` and `Σ wᵢ dᵢ`.
struct PassengerSums
    pll_weight::Float64
    pll::Float64
    fll_weight::Float64
    fll::typeof(0.0Hz)
    dll_weight::Float64
    dll::Float64
end

Base.:+(a::PassengerSums, b::PassengerSums) = PassengerSums(
    a.pll_weight + b.pll_weight,
    a.pll + b.pll,
    a.fll_weight + b.fll_weight,
    a.fll + b.fll,
    a.dll_weight + b.dll_weight,
    a.dll + b.dll,
)
Base.zero(::Type{PassengerSums}) = PassengerSums(0.0, 0.0, 0.0, 0.0Hz, 0.0, 0.0)

# One contribution's share of the sums, its DLL reading referred to the driver's
# code phase. A zero code weight (unknown group delay, `NaN` offset) adds nothing.
@inline _passenger_sum(c::PassengerContribution) = PassengerSums(
    c.pll_weight,
    c.pll_weight * c.pll,
    c.fll_weight,
    c.fll_weight * c.fll,
    c.dll_weight,
    iszero(c.dll_weight) ? 0.0 : c.dll_weight * (c.dll - c.dll_offset),
)

# Apply record `k` of every coincident passenger; returns the passengers and each
# one's contribution, in passenger order.
@inline function _apply_passenger_records(
    passenger_contexts::PassengerContexts,
    passengers::Tuple,
    args...,
)
    applied = map(passengers, passenger_contexts.contexts) do passenger, context
        _apply_coincident_passenger_record(passenger, context, args...)
    end
    map(first, applied), map(last, applied)
end

# The driver's discriminators with the passengers' contributions mixed in, loop
# by loop (`_combined_discriminator`); the driver's FLL weight is zero where
# there is no previous prompt yet.
@inline function _mixed_discriminators(
    contributions::Tuple,
    signal::AbstractGNSSSignal,
    own,
    previous_prompt,
)
    sums = foldl(+, map(_passenger_sum, contributions); init = zero(PassengerSums))
    weight = get_relative_power(signal)
    (;
        pll = _combined_discriminator(own.pll, weight, sums.pll, sums.pll_weight),
        fll = _combined_discriminator(
            own.fll,
            iszero(previous_prompt) ? 0.0 : weight,
            sums.fll,
            sums.fll_weight,
        ),
        dll = _combined_discriminator(own.dll, weight, sums.dll, sums.dll_weight),
    )
end

# Hand the passengers back with the coincident ones' records consumed: they were
# applied in the driver fold.
@inline function _finish_passengers(
    passenger_contexts::PassengerContexts,
    passengers::Tuple,
)
    foreach(passengers, passenger_contexts.contexts) do passenger, context
        context.coincident && empty!(passenger.correlator_outputs)
    end
    passengers
end

# One loop's discriminator: the weighted mean of the driver's own value and the
# passengers' sums, each signal weighted by its power share. With no passenger
# weight it is the driver's value as is (`(w·d)/w` is not `d` bit for bit), so
# a satellite with nothing to combine closes its loops as the conventional
# estimator does.
@inline _combined_discriminator(own, own_weight, sum, weight) =
    iszero(weight) ? own : (own_weight * own + sum) / (own_weight + weight)

@inline function _apply_coincident_passenger_record(
    tracked_signal::TrackedSignal,
    context::NamedTuple,
    k::Int,
    prn::Integer,
    sampling_frequency,
    code_doppler,
    driver_carrier_phase::Real,
)
    context.coincident || return tracked_signal, _NO_CONTRIBUTION
    output = @inbounds tracked_signal.correlator_outputs[k]
    ts, filtered_correlator, previous_prompt = _apply_passenger_record(
        tracked_signal,
        output,
        context.found_before_fold,
        prn,
        sampling_frequency,
        context.noise_density,
        context.noise_density_ready,
        driver_carrier_phase,
    )
    signal = ts.signal
    weight = get_relative_power(signal)
    # `fll_disc` forms `conj(p_prev) · p` and `dll_disc` reads tap magnitudes,
    # so a constant rotation cancels from both; only the PLL needs it.
    derotated = @inline apply(tap -> tap * context.derotation, filtered_correlator)
    discriminators = _record_discriminators(
        signal,
        filtered_correlator,
        previous_prompt,
        output.integrated_samples / sampling_frequency,
        code_doppler,
        sampling_frequency,
    )
    ts,
    PassengerContribution(
        true,
        weight,
        pll_disc(signal, derotated),
        iszero(previous_prompt) ? 0.0 : weight,
        discriminators.fll,
        # A `NaN` offset is an unknown group delay: no code contribution.
        isnan(context.dll_offset) ? 0.0 : weight,
        discriminators.dll,
        context.dll_offset,
    )
end

# ---------------------------------------------------------------------------
# Group delays
# ---------------------------------------------------------------------------

# The one conversion from the public form of a group delay to the stored one
# (see `SatConventionalPLLAndDLL`): `nothing` becomes the `NaN` sentinel, any time
# unit becomes `typeof(1.0s)`, and a non-time is refused by `convert`. NaN itself
# is refused, being the stored form of "unknown".
_stored_group_delay(::Nothing) = _UNKNOWN_GROUP_DELAY
function _stored_group_delay(delay::Number)
    converted = convert(typeof(1.0s), delay)
    isnan(converted) && throw(
        ArgumentError("a group delay cannot be NaN; pass `nothing` for an unknown one"),
    )
    converted
end

"""
$(SIGNATURES)

The group delay of one signal on one satellite, as a time, or `nothing` while
unknown (see [`set_group_delay!`](@ref)). Addressed like the setter, by group,
satellite and signal; the per-satellite forms take the signal (by type or
index) or, for the per-satellite state, its index.
"""
get_group_delay(state::SatConventionalPLLAndDLL, signal_index::Integer) =
    _group_delay(state, signal_index)
get_group_delay(sat::TrackedSat, sig::_SignalSelector) =
    _group_delay(sat.doppler_estimator_state, _signal_index(sat.signals, sig))
get_group_delay(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    sig::_SignalSelector,
) = get_group_delay(get_sat_state(track_state, group, sat_id), sig)

"""
$(SIGNATURES)

Set the payload **group delay** of one signal on one satellite, as a time
(`1.2e-9s`, `-0.3u"ns"`), or `nothing` to mark it unknown again. Every signal
starts at `nothing`, and the value survives [`reset_loop_filters!`](@ref). Only
an estimator with `combine_discriminators` reads it (see
[`ConventionalPLLAndDLL`](@ref)); the others store it all the same.

It is read only as a **difference** against the satellite's estimator-driver
signal (`signals[1]`), so the datum is yours to choose: put `0.0s` on the driver
and each passenger's value is its bias relative to it. Positive means delayed —
the signal arrives later and sits at a smaller code phase. A passenger's code
discriminator has `(delay(signals[1]) − delay(passenger)) · f_code` chips
subtracted before it is combined, so the shared `code_phase` keeps meaning the
driver's code phase.

`nothing` is not `0.0s`: a passenger whose delay is unknown stays out of the
combined **code** loop and still aids the carrier loops; an unknown delay on
`signals[1]` keeps every passenger out of the code loop.

```julia
set_group_delay!(ts, :galileo_e1, 11, GalileoE1B, 0.0s)  # (group, prn, signal)
set_group_delay!(ts, :gps_l5, 7, 2, -0.5u"ns")           # (group, prn, index)
```

As for the other per-signal setters, the group is always named, even on a
single-group `TrackState`.

Mutates `track_state` in place and returns it.
"""
function set_group_delay!(
    track_state::TrackState{<:SignalGroups},
    group::Union{Symbol,Integer,Val},
    sat_id::Integer,
    sig::_SignalSelector,
    delay::Maybe{Number},
)
    sats = get_sat_states(track_state, group)
    sats[sat_id] = _set_sat_group_delay(sats[sat_id], _stored_group_delay(delay), sig)
    track_state
end

function _set_sat_group_delay(sat::TrackedSat, delay::typeof(1.0s), sel...)
    i = _signal_index(sat.signals, sel...)
    state = sat.doppler_estimator_state
    TrackedSat(
        sat;
        doppler_estimator_state = _with_group_delays(
            state,
            Base.setindex(_group_delays(state), delay, i),
        ),
    )
end

# Signal `i`'s group delay in its public form: `nothing` for the stored `NaN`.
@inline function _group_delay(state, i::Integer)
    delay = _group_delays(state)[i]
    isnan(delay) ? nothing : delay
end
