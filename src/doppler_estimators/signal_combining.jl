# Signal combining: the passengers (`signals[2:end]`) of a satellite whose
# `SignalGroup` enables it mix their discriminators into the driver's
# (`signals[1]`) before its loop filters read them. See the "Signal combining"
# section of docs/src/loop_filter.md.

"""
$(TYPEDEF)

One loop's weighted sum of discriminators and its summed weights (see
`_discriminator_weight`).

$(TYPEDFIELDS)
"""
struct WeightedSum{S,W}
    """
    Weighted sum of the discriminators
    """
    sum::S
    """
    Summed weights
    """
    weight::W
end

Base.convert(::Type{WeightedSum{S,W}}, ws::WeightedSum) where {S,W} =
    WeightedSum{S,W}(ws.sum, ws.weight)

# One more discriminator `reading` with weight `weight` added.
@inline _accumulated(ws::WeightedSum, reading, weight) =
    WeightedSum(ws.sum + weight * reading, ws.weight + weight)

"""
$(TYPEDEF)

The passengers' weighted discriminators pending for the driver's next record,
one [`Tracking.WeightedSum`](@ref) per loop; carried across chunks.

$(TYPEDFIELDS)
"""
struct SignalCombiningSums
    """
    PLL discriminators, in cycles
    """
    pll::WeightedSum{typeof(1.0s),typeof(1.0s)}
    """
    FLL discriminators
    """
    fll::WeightedSum{typeof(1.0Hz * 1.0s^3),typeof(1.0s^3)}
    """
    DLL discriminators referred to the driver, in chips
    """
    dll::WeightedSum{typeof(1.0s),typeof(1.0s)}
end

SignalCombiningSums() = SignalCombiningSums(
    WeightedSum(0.0s, 0.0s),
    WeightedSum(0.0Hz * 0.0s^3, 0.0s^3),
    WeightedSum(0.0s, 0.0s),
)

# At the sync phase snap: passenger records left pending ended before the driver's
# re-integration from the snapped phase starts, so their sums are dropped with
# the driver's in-flight integration.
@inline _drop_pending_combining(state, combine_signals::Bool) =
    combine_signals ?
    _with_loop_state(state; signal_combining_sums = SignalCombiningSums()) : state

# Which of a satellite's loops passengers may be combined into, by its per-sat
# state: all three unless the state's estimator closes some of them otherwise.
@inline _locally_closed_loops(state) = (pll = true, fll = true, dll = true)

# The loops the passengers are combined into for one fold, the one place this is
# decided: none without `combine_signals`, else those `_locally_closed_loops`
# leaves open, the FLL only while it is formed. The passengers read the
# two-quadrant discriminators, so the loop closures combine their sums into a
# four-quadrant one only within that range (see `_gated_mean`).
@inline function _loops_to_combine(state, combine_signals::Bool)
    loops = _locally_closed_loops(state)
    (
        pll = combine_signals && loops.pll,
        fll = combine_signals &&
              loops.fll &&
              _uses_fll(state.carrier_loop_filter) &&
              !state.frequency_lock.locked,
        dll = combine_signals && loops.dll,
    )
end

# Whether the estimator also measures on the passengers' records, which then come
# through the driver fold even without combining, and its per-passenger
# `(code, carrier)` `(count, sum)` accumulator pairs. By default it does not and
# has none.
@inline _measures_passengers(state) = false
@inline _passenger_measurement_accs(state, passengers::Tuple) =
    map(_ -> nothing, passengers)
@inline _with_passenger_measurements(state, _) = state

# One passenger record's raw readings added to its accumulator pair: the DLL
# against the shared replica, and the FLL four-quadrant where its prompt is wiped
# off. Nothing without accumulators or while not measuring.
@inline _measured(::Nothing, _...) = nothing
@inline function _measured(
    (code_acc, carrier_acc)::Tuple,
    measure::Bool,
    signal::AbstractGNSSSignal,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    context,
    code_doppler,
    sampling_frequency,
)
    measure || return (code_acc, carrier_acc)
    (
        _accumulated(
            code_acc,
            dll_disc(signal, correlator, code_doppler, sampling_frequency),
        ),
        _accumulated_fll(
            carrier_acc,
            fll_disc(
                signal,
                correlator,
                previous_prompt,
                integration_time;
                four_quadrant = context.wiped_off,
            ),
            previous_prompt,
        ),
    )
end

# A weighted mean of the driver's own discriminator and the passengers' sum. With
# no passenger weight it is the driver's own as is: `(w · d) / w` is not `d` bit
# for bit.
@inline _weighted_mean(own, own_weight, pending::WeightedSum) =
    iszero(pending.weight) ? own :
    (own_weight * own + pending.sum) / (own_weight + pending.weight)

# The weighted mean, but only while the driver's own reading lies within `range`,
# the passengers' two-quadrant range: there a four-quadrant driver reads the error
# as they do, beyond it they would fold it by half a cycle. A two-quadrant
# driver's reading never leaves the range.
@inline _gated_mean(own, own_weight, pending::WeightedSum, range) =
    abs(own) < range ? _weighted_mean(own, own_weight, pending) : own

# The two-quadrant discriminators' ranges, which the passengers read and which
# `_gated_mean` keeps a four-quadrant driver within: ±1/4 cycle for the PLL's
# `atan(Q / I)`, ±1/(4T) for the FLL's `atan(cross / dot)`.
const _TWO_QUADRANT_PLL_RANGE = 0.25
@inline _two_quadrant_fll_range(integration_time) = uconvert(Hz, 1 / (4 * integration_time))

# The weight of one record's discriminator: its signal's ICD power share times its
# integration time. The FLL's is the integration time cubed, as its noise
# variance falls with its cube.
@inline _discriminator_weight(signal::AbstractGNSSSignal, integration_time) =
    get_relative_power(signal) * uconvert(s, integration_time)
@inline _fll_discriminator_weight(signal::AbstractGNSSSignal, integration_time) =
    get_relative_power(signal) * uconvert(s, integration_time)^3

# Per-passenger constants of one fold. `differential_group_delay_chips`, the
# passenger's group delay minus the driver's in chips, refers the passenger's DLL
# discriminator to the driver's code phase: a passenger with the larger group
# delay arrives later, so it reads the driver's code error minus that difference.
# It is `NaN` where either group delay is unknown.
@inline function _passenger_context(
    passenger::TrackedSignal,
    (noise_density, noise_density_ready),
    driver::TrackedSignal,
    driver_carrier_phase_offset::Real,
    code_frequency,
)
    known = !isnothing(passenger.group_delay) && !isnothing(driver.group_delay)
    differential_group_delay_chips =
        known ?
        ustrip(NoUnits, (passenger.group_delay - driver.group_delay) * code_frequency) : NaN
    found_before_fold = has_bit_or_secondary_code_been_found(passenger.bit_buffer)
    (;
        found_before_fold,
        wiped_off = _is_wiped_off(passenger.signal, found_before_fold),
        noise_density,
        noise_density_ready,
        derotation = _carrier_phase_derotation(
            driver_carrier_phase_offset,
            passenger.signal,
        ),
        differential_group_delay_chips,
    )
end

# Apply each passenger's records that end by sample `until`, from its cursor on,
# adding their discriminators to `sums` and their raw readings to its
# measurement accumulators. Recursive over the tuples for type stability.
@inline _advance_passengers(
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Tuple{},
    ::Integer,
    sums::SignalCombiningSums,
    ::Vararg{Any,N},
) where {N} = ((), (), (), sums)
@inline function _advance_passengers(
    passengers::Tuple,
    cursors::Tuple,
    measurements::Tuple,
    contexts::Tuple,
    until::Integer,
    sums::SignalCombiningSums,
    args::Vararg{Any,N},
) where {N}
    passenger, cursor, measurement, sums = _advance_passenger(
        first(passengers),
        first(cursors),
        first(measurements),
        first(contexts),
        until,
        sums,
        args...,
    )
    rest, rest_cursors, rest_measurements, sums = _advance_passengers(
        Base.tail(passengers),
        Base.tail(cursors),
        Base.tail(measurements),
        Base.tail(contexts),
        until,
        sums,
        args...,
    )
    (passenger, rest...),
    (cursor, rest_cursors...),
    (measurement, rest_measurements...),
    sums
end

@inline function _advance_passenger(
    tracked_signal::TrackedSignal,
    cursor::Int,
    measurement,
    context,
    until::Integer,
    sums::SignalCombiningSums,
    loops_to_combine,
    measure::Bool,
    prn::Integer,
    sampling_frequency,
    driver_carrier_phase_offset::Real,
    code_doppler,
)
    ts = tracked_signal
    outputs = ts.correlator_outputs
    @inbounds while cursor <= length(outputs) && outputs[cursor].sample_index <= until
        output = outputs[cursor]
        previous_prompt = _fll_previous_prompt(ts, output, sampling_frequency)
        integration_time = output.integrated_samples / sampling_frequency
        ts, filtered_correlator = _apply_passenger_record(
            ts,
            output,
            context.found_before_fold,
            prn,
            sampling_frequency,
            context.noise_density,
            context.noise_density_ready,
            driver_carrier_phase_offset,
        )
        sums = _add_passenger_discriminators(
            sums,
            ts.signal,
            filtered_correlator,
            previous_prompt,
            integration_time,
            context,
            loops_to_combine,
            code_doppler,
            sampling_frequency,
        )
        measurement = _measured(
            measurement,
            measure,
            ts.signal,
            filtered_correlator,
            previous_prompt,
            integration_time,
            context,
            code_doppler,
            sampling_frequency,
        )
        cursor += 1
    end
    ts, cursor, measurement, sums
end

# One passenger record's weighted discriminators, formed only for the loops it
# is combined into; its PLL is read on the driver's carrier phase frame. Its
# two-quadrant carrier discriminators are blind to a sign flip of a whole record,
# so a record correlated before the passenger's own sync counts like any other.
@inline function _add_passenger_discriminators(
    sums::SignalCombiningSums,
    signal::AbstractGNSSSignal,
    correlator::AbstractCorrelator,
    previous_prompt::Complex,
    integration_time,
    context,
    loops_to_combine,
    code_doppler,
    sampling_frequency,
)
    weight = _discriminator_weight(signal, integration_time)
    pll_weight = loops_to_combine.pll ? weight : zero(weight)
    pll =
        loops_to_combine.pll ?
        pll_disc(
            signal,
            update_accumulator(
                correlator,
                get_accumulators(correlator) .* context.derotation,
            ),
        ) : 0.0
    fll_weight =
        loops_to_combine.fll && !iszero(previous_prompt) ?
        _fll_discriminator_weight(signal, integration_time) : 0.0s^3
    fll =
        iszero(fll_weight) ? 0.0Hz :
        fll_disc(signal, correlator, previous_prompt, integration_time)
    dll_weight =
        loops_to_combine.dll && !isnan(context.differential_group_delay_chips) ? weight :
        zero(weight)
    dll =
        iszero(dll_weight) ? 0.0 :
        dll_disc(signal, correlator, code_doppler, sampling_frequency) +
        context.differential_group_delay_chips
    SignalCombiningSums(
        _accumulated(sums.pll, pll, pll_weight),
        _accumulated(sums.fll, fll, fll_weight),
        _accumulated(sums.dll, dll, dll_weight),
    )
end
