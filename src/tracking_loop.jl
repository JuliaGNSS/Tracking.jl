
"""
$(SIGNATURES)

Track the signal `signal` based on the current tracking `state`, the sample frequency
`sample_frequency` and PRN `prn`. Optional configurations are:
- Post correlation filter `post_corr_filter` defaults to `get_default_post_corr_filter(...)`
- Intermediate frequency `intermediate_frequency` defaults to 0Hz
- Maximal integration time `max_integration_time` defaults to 1ms. The actual integration
  time can be less in order to find the secondary code or the bit shift. Once this is found
  the integration time is identical to the maximal integration time.
- Minimal integration time `min_integration_time` defaults to 0.75ms. It's the minimal
  integration time, that leads to a valid correlation result. It is only used for the first
  integration period.
- Sample shift between early and late `early_late_sample_shift` defaults to
  `get_early_late_sample_shift(...)`
- Bandwidth of the carrier loop `carrier_loop_filter_bandwidth` defaults to 18Hz
- Bandwidth of the code loop `code_loop_filter_bandwidth` defaults to 1Hz
- Velocity aiding `velocity_aiding` defaults to 0Hz
"""
function track(
        signal,
        state::TrackingState{S, C, CALF, COLF, CN},
        prn::Integer,
        sample_frequency;
        post_corr_filter = get_default_post_corr_filter(get_correlator(state)),
        intermediate_frequency = 0.0Hz,
        max_integration_time::typeof(1ms) = 1ms,
        min_integration_time::typeof(1.0ms) = 0.75ms,
        early_late_sample_shift = get_early_late_sample_shift(S,
            get_correlator(state), sample_frequency, 0.5),
        carrier_loop_filter_bandwidth = 18Hz,
        code_loop_filter_bandwidth = 1Hz,
        velocity_aiding = 0Hz
) where {
    S <: AbstractGNSSSystem,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator
}
    if get_data_frequency(S) != 0Hz
        @assert rem(1 / get_data_frequency(S), max_integration_time) == 0ms
    end
    init_carrier_doppler = get_init_carrier_doppler(state)
    init_code_doppler = get_init_code_doppler(state)
    carrier_doppler = get_carrier_doppler(state)
    code_doppler = get_code_doppler(state)
    carrier_phase = get_carrier_phase(state)
    code_phase = get_code_phase(state)
    sc_bit_detector = get_sc_bit_detector(state)
    correlator = get_correlator(state)
    carrier_loop_filter = get_carrier_loop_filter(state)
    code_loop_filter = get_code_loop_filter(state)
    prompt_accumulator = get_prompt_accumulator(state)
    integrated_samples = get_integrated_samples(state)
    cn0_estimator = get_cn0_estimator(state)
    signal_start_sample = 1
    bit_buffer = BitBuffer()
    valid_correlator = zero(correlator)
    got_correlator = false
    while true
        num_samples_left_to_integrate = get_num_samples_left_to_integrate(
            S,
            max_integration_time,
            sample_frequency,
            code_doppler,
            code_phase,
            found(sc_bit_detector)
        )
        signal_samples_left = get_num_samples(signal) - signal_start_sample + 1
        num_samples_left = min(num_samples_left_to_integrate, signal_samples_left)
        carrier_frequency = get_current_carrier_frequency(
            intermediate_frequency,
            carrier_doppler
        )
        code_frequency = get_current_code_frequency(S, code_doppler)
        correlator = correlate(
            S,
            correlator,
            signal,
            prn,
            early_late_sample_shift,
            signal_start_sample,
            num_samples_left,
            carrier_frequency,
            code_frequency,
            sample_frequency,
            carrier_phase,
            code_phase
        )
        integrated_samples += num_samples_left
        carrier_phase = update_carrier_phase(
            num_samples_left,
            carrier_frequency,
            sample_frequency,
            carrier_phase
        )
        prev_code_phase = code_phase
        code_phase = update_code_phase(
            S,
            num_samples_left,
            code_frequency,
            sample_frequency,
            code_phase,
            found(sc_bit_detector)
        )
        integration_time = integrated_samples / sample_frequency
        if num_samples_left == num_samples_left_to_integrate &&
                integration_time >= min_integration_time
            got_correlator = true

            correlator = normalize(correlator, integrated_samples)
            valid_correlator = correlator
            filtered_correlator = filter(post_corr_filter, correlator)
            pll_discriminator = pll_disc(S, filtered_correlator)
            dll_discriminator = dll_disc(
                S,
                filtered_correlator,
                early_late_sample_shift,
                code_frequency / sample_frequency
            )
            carrier_freq_update, carrier_loop_filter = filter_loop(
                carrier_loop_filter,
                pll_discriminator,
                integration_time,
                carrier_loop_filter_bandwidth
            )
            code_freq_update, code_loop_filter = filter_loop(
                code_loop_filter,
                dll_discriminator,
                integration_time,
                code_loop_filter_bandwidth
            )
            carrier_doppler, code_doppler = aid_dopplers(
                S,
                init_carrier_doppler,
                init_code_doppler,
                carrier_freq_update,
                code_freq_update,
                velocity_aiding
            )
            cn0_estimator = update(cn0_estimator, get_prompt(filtered_correlator))
            bit_buffer, prompt_accumulator = buffer(
                S,
                bit_buffer,
                prompt_accumulator,
                found(sc_bit_detector),
                prev_code_phase,
                code_phase,
                max_integration_time,
                get_prompt(filtered_correlator)
            )
            sc_bit_detector = find(S, sc_bit_detector, get_prompt(filtered_correlator))
            correlator = zero(correlator)
            integrated_samples = 0
        end

        num_samples_left == signal_samples_left && break
        signal_start_sample += num_samples_left
    end
    next_state = TrackingState{S, C, CALF, COLF, CN}(
        init_carrier_doppler,
        init_code_doppler,
        carrier_doppler,
        code_doppler,
        carrier_phase,
        code_phase,
        correlator,
        carrier_loop_filter,
        code_loop_filter,
        sc_bit_detector,
        integrated_samples,
        prompt_accumulator,
        cn0_estimator
    )
    estimated_cn0 = estimate_cn0(cn0_estimator, max_integration_time)
    TrackingResults(next_state, valid_correlator, got_correlator, bit_buffer, estimated_cn0)
end

"""
$(SIGNATURES)

Default post corr filter
"""
function get_default_post_corr_filter(correlator::AbstractCorrelator{T}) where T <: SVector
    x -> x[1]
end

function get_default_post_corr_filter(correlator)
    x -> x
end

"""
$(SIGNATURES)

Returns the appropiate integration time. It will be the maximum integration time once the
secondary code or the bit shift has been found.
"""
function get_integration_time(
    ::Type{S},
    max_integration_time,
    secondary_code_or_bit_found::Bool
) where S <: AbstractGNSSSystem
    ifelse(
        secondary_code_or_bit_found,
        max_integration_time,
        min(
            convert(typeof(1ms), get_code_length(S) / get_code_frequency(S)),
            max_integration_time
        )
    )
end

"""
$(SIGNATURES)

Calculates the number of chips to integrate.
"""
function get_num_chips_to_integrate(
    ::Type{S},
    max_integration_time,
    current_code_phase,
    secondary_code_or_bit_found
) where S <: AbstractGNSSSystem
    max_phase = Int(upreferred(get_code_frequency(S) *
        get_integration_time(S, max_integration_time, secondary_code_or_bit_found)))
    current_phase_mod_max_phase = mod(current_code_phase, max_phase)
    max_phase - current_phase_mod_max_phase
end

"""
$(SIGNATURES)

Calculates the number of samples to integrate.
"""
function get_num_samples_left_to_integrate(
    ::Type{S},
    max_integration_time,
    sample_frequency,
    current_code_doppler,
    current_code_phase,
    secondary_code_or_bit_found
) where S <: AbstractGNSSSystem
    phase_to_integrate = get_num_chips_to_integrate(
        S,
        max_integration_time,
        current_code_phase,
        secondary_code_or_bit_found
    )
    code_frequency = get_code_frequency(S) + current_code_doppler
    ceil(Int, phase_to_integrate * sample_frequency / code_frequency)
end

"""
$(SIGNATURES)

Calculates the current carrier frequency.
"""
function get_current_carrier_frequency(intermediate_frequency, carrier_doppler)
    carrier_doppler + intermediate_frequency
end

"""
$(SIGNATURES)

Calculates the current code frequency.
"""
function get_current_code_frequency(::Type{S}, code_doppler) where S <: AbstractGNSSSystem
    code_doppler + get_code_frequency(S)
end

"""
$(SIGNATURES)

Updates the carrier phase.
"""
function update_carrier_phase(
    num_samples,
    carrier_frequency,
    sample_frequency,
    start_carrier_phase
)
    mod(
        carrier_frequency * num_samples / sample_frequency + start_carrier_phase + 0.5,
        1
    ) - 0.5
end

"""
$(SIGNATURES)

Updates the code phase.
"""
function update_code_phase(
    ::Type{S},
    num_samples,
    code_frequency,
    sample_frequency,
    start_code_phase,
    secondary_code_or_bit_found
) where S <: AbstractGNSSSystem
    secondary_code_or_bit_length =
        get_data_frequency(S) == 0Hz ?
        get_secondary_code_length(S) :
        Int(get_code_frequency(S) / (get_data_frequency(S) * get_code_length(S)))
    code_length = get_code_length(S) *
        (secondary_code_or_bit_found ? secondary_code_or_bit_length : 1)
    mod(code_frequency * num_samples / sample_frequency + start_code_phase, code_length)
end

"""
$(SIGNATURES)

Correlate the signal the the replicas.
"""
function correlate(
    ::Type{S},
    correlator,
    signal,
    prn::Integer,
    early_late_sample_shift::Integer,
    start_sample::Integer,
    num_samples::Integer,
    carrier_frequency,
    code_frequency,
    sample_frequency,
    start_carrier_phase,
    start_code_phase
) where S <: AbstractGNSSSystem
    total_code_length = get_code_length(S) * get_secondary_code_length(S)
    start_code_phase = mod(start_code_phase, total_code_length)
    carrier_phase_wrap = 0
    code_phase_wrap = 0
    code_register = 0
    max_early_late_sample_shift = maximum(early_late_sample_shift)
    for i = -max_early_late_sample_shift:max_early_late_sample_shift-1
        code = get_code(S, code_frequency * i / sample_frequency + start_code_phase, prn)
        code_register = code_register << 1 + (code + 1) >> 1
    end
    @inbounds @fastmath for i = 0:num_samples - 1
        carrier_phase = carrier_frequency * i / sample_frequency + start_carrier_phase
        carrier_phase_wrap += (carrier_phase - carrier_phase_wrap >= 0.5) * 1
        carrier = cis(-2π * (carrier_phase - carrier_phase_wrap))
        current_signal = get_signal(correlator, signal, start_sample + i)
        earliest_code_phase = code_frequency * (i + max_early_late_sample_shift) /
            sample_frequency + start_code_phase
        code_phase_wrap += (earliest_code_phase - code_phase_wrap >= total_code_length) *
            total_code_length
        earliest_code = get_code_unsafe(S, earliest_code_phase - code_phase_wrap, prn)
        code_register = code_register << 1 + (earliest_code + 1) >> 1
        correlator = correlate_iteration(
            S,
            correlator,
            current_signal,
            carrier,
            code_register,
            early_late_sample_shift,
            earliest_code
        )
    end
    correlator
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(
    ::Type{S},
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
    velocity_aiding
) where S <: AbstractGNSSSystem
    carrier_doppler = carrier_freq_update + velocity_aiding
    code_doppler = code_freq_update + carrier_doppler * get_code_center_frequency_ratio(S)
    init_carrier_doppler + carrier_doppler, init_code_doppler + code_doppler
end

"""
$(SIGNATURES)

Get the signal at sample `sample`.
"""
Base.@propagate_inbounds @inline function get_signal(
    correlator::AbstractCorrelator{T},
    signal,
    sample::Integer
) where T <: SVector
    T(view(signal, :, sample))
end

Base.@propagate_inbounds @inline function get_signal(correlator, signal, sample::Integer)
    signal[sample]
end

"""
$(SIGNATURES)

Get the number of samples in the signal.
"""
@inline function get_num_samples(signal)
    length(signal)
end

@inline function get_num_samples(signal::AbstractMatrix)
    size(signal, 2)
end
