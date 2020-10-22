"""
$(SIGNATURES)

Track the signal `signal` based on the current tracking `state`, the sampling frequency
`sampling_frequency` and PRN `prn`. Optional configurations are:
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
        gain_controlled_signal::GainControlledSignal,
        state::TrackingState{S, C, CALF, COLF, CN, DS},
        prn::Integer,
        sampling_frequency;
        post_corr_filter = get_default_post_corr_filter(get_correlator(state)),
        intermediate_frequency = 0.0Hz,
        max_integration_time::typeof(1ms) = 1ms,
        min_integration_time::typeof(1.0ms) = 0.75ms,
        early_late_sample_shift = get_early_late_sample_shift(S,
            get_correlator(state), sampling_frequency, 0.5),
        carrier_loop_filter_bandwidth = 18Hz,
        code_loop_filter_bandwidth = 1Hz,
        velocity_aiding = 0Hz,
        carrier_amplitude_power::Val{N} = Val(5)
) where {
    S <: AbstractGNSSSystem,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator,
    DS <: StructArray,
    N
}
    if get_data_frequency(S) != 0Hz
        @assert rem(1 / get_data_frequency(S), max_integration_time) == 0ms
    end
    N > 7 && throw(ArgumentError("The carrier amplitude power should be less than 8 to stay within 16 bits."))
    (get_amplitude_power(gain_controlled_signal) + N) > 16 && throw(ArgumentError("The AGC amplitude + carrier replica amplitude should not exceed 16 bits"))
    signal = get_signal(gain_controlled_signal)
    correlator = get_correlator(state)
    num_ants = get_num_ants(correlator)
    size(signal, 2) == num_ants || throw(ArgumentError("The second dimension of the signal should be equal to the number of antennas specified by num_ants = NumAnts(N) in the TrackingState."))
    agc_amplitude_power = get_amplitude_power(gain_controlled_signal)
    agc_attenuation = get_attenuation(gain_controlled_signal)
    downconverted_signal = resize!(get_downconverted_signal(state), size(signal, 1))
    carrier_replica = resize!(get_carrier(state), size(signal, 1))
    code_replica = resize!(get_code(state), size(signal, 1) + 2 * maximum(early_late_sample_shift))
    init_carrier_doppler = get_init_carrier_doppler(state)
    init_code_doppler = get_init_code_doppler(state)
    carrier_doppler = get_carrier_doppler(state)
    code_doppler = get_code_doppler(state)
    carrier_phase = get_carrier_phase(state)
    code_phase = get_code_phase(state)
    sc_bit_detector = get_sc_bit_detector(state)
    carrier_loop_filter = get_carrier_loop_filter(state)
    code_loop_filter = get_code_loop_filter(state)
    prompt_accumulator = get_prompt_accumulator(state)
    integrated_samples = get_integrated_samples(state)
    cn0_estimator = get_cn0_estimator(state)
    signal_start_sample = 1
    bit_buffer = BitBuffer()
    valid_correlator = zero(correlator)
    valid_correlator_carrier_phase = 0.0
    valid_correlator_carrier_frequency = 0.0Hz
    got_correlator = false
    while true
        num_samples_left_to_integrate = get_num_samples_left_to_integrate(
            S,
            max_integration_time,
            sampling_frequency,
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
        code_replica = gen_code_replica!(
            code_replica,
            S,
            code_frequency,
            sampling_frequency,
            code_phase,
            signal_start_sample,
            num_samples_left,
            early_late_sample_shift,
            prn
        )
        carrier_replica = gen_carrier_replica!(
            carrier_replica,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
            carrier_amplitude_power,
            signal_start_sample,
            num_samples_left
        )
        downconverted_signal = downconvert!(
            downconverted_signal,
            signal,
            carrier_replica,
            signal_start_sample,
            num_samples_left
        )
        correlator = correlate(
            correlator,
            downconverted_signal,
            code_replica,
            early_late_sample_shift,
            signal_start_sample,
            num_samples_left,
            agc_attenuation,
            agc_amplitude_power,
            carrier_amplitude_power
        )
        integrated_samples += num_samples_left
        carrier_phase = update_carrier_phase(
            num_samples_left,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
            carrier_amplitude_power
        )
        prev_code_phase = code_phase
        code_phase = update_code_phase(
            S,
            num_samples_left,
            code_frequency,
            sampling_frequency,
            code_phase,
            found(sc_bit_detector)
        )
        integration_time = integrated_samples / sampling_frequency
        if num_samples_left == num_samples_left_to_integrate &&
                integration_time >= min_integration_time
            got_correlator = true

            correlator = normalize(correlator, integrated_samples)
            valid_correlator = correlator
            valid_correlator_carrier_phase = carrier_phase
            valid_correlator_carrier_frequency = carrier_frequency
            filtered_correlator = filter(post_corr_filter, correlator)
            pll_discriminator = pll_disc(S, filtered_correlator)
            dll_discriminator = dll_disc(
                S,
                filtered_correlator,
                early_late_sample_shift,
                code_frequency / sampling_frequency
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
    next_state = TrackingState{S, C, CALF, COLF, CN, DS}(
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
        cn0_estimator,
        downconverted_signal,
        carrier_replica,
        code_replica
    )
    estimated_cn0 = estimate_cn0(cn0_estimator, max_integration_time)
    TrackingResults(
        next_state,
        valid_correlator,
        valid_correlator_carrier_frequency,
        valid_correlator_carrier_phase,
        got_correlator,
        bit_buffer,
        estimated_cn0
    )
end

@inline function track(
        signal::AbstractArray,
        state::TrackingState{S, C, CALF, COLF, CN, DS},
        prn::Integer,
        sampling_frequency;
        post_corr_filter = get_default_post_corr_filter(get_correlator(state)),
        intermediate_frequency = 0.0Hz,
        max_integration_time::typeof(1ms) = 1ms,
        min_integration_time::typeof(1.0ms) = 0.75ms,
        early_late_sample_shift = get_early_late_sample_shift(S,
            get_correlator(state), sampling_frequency, 0.5),
        carrier_loop_filter_bandwidth = 18Hz,
        code_loop_filter_bandwidth = 1Hz,
        velocity_aiding = 0Hz,
        carrier_amplitude_power::Val{N} = Val(5)
) where {
    S <: AbstractGNSSSystem,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator,
    DS <: StructArray,
    N
}
    correlator = get_correlator(state)
    num_ants = get_num_ants(correlator)
    size(signal, 2) == num_ants || throw(ArgumentError("The second dimension of the signal should be equal to the number of antennas specified by num_ants = NumAnts(N) in the TrackingState."))
    track(
        GainControlledSignal(signal),
        state,
        prn,
        sampling_frequency,
        post_corr_filter = post_corr_filter,
        intermediate_frequency = intermediate_frequency,
        max_integration_time = max_integration_time,
        min_integration_time = min_integration_time,
        early_late_sample_shift = early_late_sample_shift,
        carrier_loop_filter_bandwidth = carrier_loop_filter_bandwidth,
        code_loop_filter_bandwidth = code_loop_filter_bandwidth,
        velocity_aiding = velocity_aiding,
        carrier_amplitude_power = carrier_amplitude_power
    )
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
    sampling_frequency,
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
    ceil(Int, phase_to_integrate * sampling_frequency / code_frequency)
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

Get the number of samples in the signal.
"""
@inline function get_num_samples(signal)
    length(signal)
end

@inline function get_num_samples(signal::AbstractMatrix)
    size(signal, 1)
end

function resize!(A::StructArray{Complex{T}, 2}, b::Integer) where T
    if size(A, 1) == b
        return A
    end
    num_ants = size(A, 2)
    StructArray{Complex{T}}((Matrix{T}(undef, b, num_ants), Matrix{T}(undef, b, num_ants)))
end
