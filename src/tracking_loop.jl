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
- Sample shift of the correlator replica `correlator_sample_shifts` defaults to
  `get_correlator_sample_shifts(...)`
- Bandwidth of the carrier loop `carrier_loop_filter_bandwidth` defaults to 18Hz
- Bandwidth of the code loop `code_loop_filter_bandwidth` defaults to 1Hz
- Velocity aiding `velocity_aiding` defaults to 0Hz
"""
function track(
        signal,
        state::TrackingState{S, C, CALF, COLF, CN, DS, CAR, COR},
        sampling_frequency;
        post_corr_filter = get_default_post_corr_filter(get_correlator(state)),
        intermediate_frequency = 0.0Hz,
        max_integration_time::typeof(1ms) = 1ms,
        min_integration_time::typeof(1.0ms) = 0.75ms,
        correlator_sample_shifts = get_correlator_sample_shifts(get_system(state),
            get_correlator(state), sampling_frequency, 0.5),
        early_late_index_shift = get_early_late_index_shift(get_system(state),
            correlator_sample_shifts, get_correlator(state), sampling_frequency, 0.5),
        carrier_loop_filter_bandwidth = 18Hz,
        code_loop_filter_bandwidth = 1Hz,
        velocity_aiding = 0Hz,
) where {
    S <: AbstractGNSS,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator,
    DS,
    CAR,
    COR
}
    prn = get_prn(state)
    system = get_system(state)
    if get_data_frequency(system) != 0Hz
        @assert rem(1 / get_data_frequency(system), max_integration_time) == 0ms
    end
    correlator = get_correlator(state)
    num_ants = get_num_ants(correlator)
    size(signal, 2) == num_ants || throw(ArgumentError("The second dimension of the signal should be equal to the number of antennas specified by num_ants = NumAnts(N) in the TrackingState."))
    if typeof(get_codes(system)) <: CuMatrix
        typeof(signal) <: StructArray || throw(ArgumentError("Signal is not a StructArray, initialize the signal properly and try again."))
        typeof(signal.re) <: CuArray && typeof(get_codes(system)) <: CuArray || throw(ArgumentError("Signal and GNSS codes are not of the same type. Please check if CPU or GPU is used."))
    end
    downconverted_signal_temp = get_downconverted_signal(state)
    downconverted_signal = resize!(downconverted_signal_temp, size(signal, 1), signal)
    carrier_replica = get_carrier(state)
    resize!(choose(carrier_replica, signal), size(signal, 1))
    code_replica = get_code(state)
    resize!(code_replica, size(signal, 1) + correlator_sample_shifts[end]-correlator_sample_shifts[1])
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
            system,
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
        code_frequency = get_current_code_frequency(system, code_doppler)
        correlator = downconvert_and_correlate!(
            system,
            signal,
            correlator,
            code_replica,
            code_phase,
            carrier_replica,
            carrier_phase,
            downconverted_signal,
            code_frequency,
            correlator_sample_shifts,
            carrier_frequency,
            sampling_frequency,
            signal_start_sample,
            num_samples_left,
            prn
        )
        integrated_samples += num_samples_left
        carrier_phase = update_carrier_phase(
            num_samples_left,
            carrier_frequency,
            sampling_frequency,
            carrier_phase,
        )
        prev_code_phase = code_phase
        code_phase = update_code_phase(
            system,
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
            pll_discriminator = pll_disc(
                system,
                filtered_correlator,
                correlator_sample_shifts
            )
            dll_discriminator = dll_disc(
                system,
                filtered_correlator,
                correlator_sample_shifts,
                early_late_index_shift,
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
                system,
                init_carrier_doppler,
                init_code_doppler,
                carrier_freq_update,
                code_freq_update,
                velocity_aiding
            )
            cn0_estimator = update(
                cn0_estimator,
                get_prompt(filtered_correlator, correlator_sample_shifts)
            )
            bit_buffer, prompt_accumulator = buffer(
                system,
                bit_buffer,
                prompt_accumulator,
                found(sc_bit_detector),
                prev_code_phase,
                code_phase,
                max_integration_time,
                get_prompt(filtered_correlator, correlator_sample_shifts)
            )
            sc_bit_detector = find(
                system,
                sc_bit_detector,
                get_prompt(filtered_correlator, correlator_sample_shifts)
            )
            correlator = zero(correlator)
            integrated_samples = 0
        end

        num_samples_left == signal_samples_left && break
        signal_start_sample += num_samples_left
    end
    next_state = TrackingState{S, C, CALF, COLF, CN, DS, CAR, COR}(
        prn,
        system,
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
        correlator_sample_shifts,
        early_late_index_shift,
        valid_correlator_carrier_frequency,
        valid_correlator_carrier_phase,
        got_correlator,
        bit_buffer,
        estimated_cn0
    )
end

function downconvert_and_correlate!(
    system,
    signal,
    correlator,
    code_replica,
    code_phase,
    carrier_replica,
    carrier_phase,
    downconverted_signal,
    code_frequency,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn
)
    gen_code_replica!(
        code_replica,
        system,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        num_samples_left,
        correlator_sample_shifts,
        prn
    )
    gen_carrier_replica!(
        choose(carrier_replica, signal),
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        num_samples_left
    )
    downconvert!(
        choose(downconverted_signal, signal),
        signal,
        choose(carrier_replica, signal),
        signal_start_sample,
        num_samples_left
    )
    correlate(
        correlator,
        choose(downconverted_signal, signal),
        code_replica,
        correlator_sample_shifts,
        signal_start_sample,
        num_samples_left
    )
end

# CUDA downconvert_and_correlate for num_ants > 1
function downconvert_and_correlate!(
    system::AbstractGNSS{C},
    signal::AbstractMatrix,
    correlator::T,
    code_replica,
    code_phase,
    carrier_replica,
    carrier_phase,
    downconverted_signal,
    code_frequency,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn
) where {C <: CuMatrix, T <: AbstractCorrelator}
    accumulator_result = downconvert_and_correlate_kernel_wrapper(
        system,
        view(signal, signal_start_sample:signal_start_sample - 1 + num_samples_left,:),
        correlator,
        code_phase,
        carrier_phase,
        code_frequency,
        correlator_sample_shifts,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        num_samples_left,
        prn
    )
    return T(map(+, get_accumulators(correlator), eachcol(Array(accumulator_result[1,:,:]))))
end

# CUDA downconvert_and_correlate for num_ants = 1
function downconvert_and_correlate!(
    system::AbstractGNSS{C},
    signal::AbstractVector,
    correlator::T,
    code_replica,
    code_phase,
    carrier_replica,
    carrier_phase,
    downconverted_signal,
    code_frequency,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn
) where {C <: CuMatrix, T <: AbstractCorrelator}
    accumulator_result = downconvert_and_correlate_kernel_wrapper(
        system,
        view(signal, signal_start_sample:signal_start_sample - 1 + num_samples_left),
        correlator,
        code_phase,
        carrier_phase,
        code_frequency,
        correlator_sample_shifts,
        carrier_frequency,
        sampling_frequency,
        signal_start_sample,
        num_samples_left,
        prn
    )
    addition(a,b) = a + first(b)
    return T(map(addition, get_accumulators(correlator), eachcol(Array(accumulator_result[1,:,:]))))
end

function choose(replica::CarrierReplicaCPU, signal::AbstractArray{Complex{Float64}})
    replica.carrier_f64
end
function choose(replica::CarrierReplicaCPU, signal::AbstractArray{Complex{T}}) where T <: Number
    replica.carrier_f32
end
function choose(replica::DownconvertedSignalCPU, signal::AbstractArray{Complex{Float64}})
    replica.downconverted_signal_f64
end
function choose(replica::DownconvertedSignalCPU, signal::AbstractArray{Complex{T}}) where T <: Number
    replica.downconverted_signal_f32
end
function choose(replica::Nothing, signal::AbstractArray)
    nothing
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
    system::AbstractGNSS,
    max_integration_time,
    secondary_code_or_bit_found::Bool
)
    ifelse(
        secondary_code_or_bit_found,
        max_integration_time,
        min(
            ceil(typeof(1ms), get_code_length(system) / get_code_frequency(system)),
            max_integration_time
        )
    )
end

"""
$(SIGNATURES)

Calculates the number of chips to integrate.
"""
function get_num_chips_to_integrate(
    system::AbstractGNSS,
    max_integration_time,
    current_code_phase,
    secondary_code_or_bit_found
)
    max_phase = Int(upreferred(get_code_frequency(system) *
        get_integration_time(system, max_integration_time, secondary_code_or_bit_found)))
    current_phase_mod_max_phase = mod(current_code_phase, max_phase)
    return max_phase - current_phase_mod_max_phase
end

"""
$(SIGNATURES)

Calculates the number of samples to integrate.
"""
function get_num_samples_left_to_integrate(
    system::AbstractGNSS,
    max_integration_time,
    sampling_frequency,
    current_code_doppler,
    current_code_phase,
    secondary_code_or_bit_found
)
    phase_to_integrate = get_num_chips_to_integrate(
        system,
        max_integration_time,
        current_code_phase,
        secondary_code_or_bit_found
    )
    code_frequency = get_code_frequency(system) + current_code_doppler
    ceil(Int, phase_to_integrate * sampling_frequency / code_frequency)
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(
    system::AbstractGNSS,
    init_carrier_doppler,
    init_code_doppler,
    carrier_freq_update,
    code_freq_update,
    velocity_aiding
)
    carrier_doppler = carrier_freq_update + velocity_aiding
    code_doppler = code_freq_update + carrier_doppler * get_code_center_frequency_ratio(system)
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

function resize!(ds::DownconvertedSignalCPU, b::Integer, signal::AbstractVector)
    resize!(choose(ds, signal), b)
    ds
end

function resize!(ds::DownconvertedSignalCPU, b::Integer, signal::AbstractMatrix{Complex{Float64}})
    num_ants = size(signal, 2)
    DownconvertedSignalCPU(
        ds.downconverted_signal_f32,
        size(ds.downconverted_signal_f64, 1) == b ?
            ds.downconverted_signal_f64 :
            StructArray{Complex{Float64}}((Matrix{Float64}(undef, b, num_ants), Matrix{Float64}(undef, b, num_ants)))
    )
end

function resize!(ds::DownconvertedSignalCPU, b::Integer, signal::AbstractMatrix{Complex{T}}) where T <: Number
    num_ants = size(signal, 2)
    DownconvertedSignalCPU(
        size(ds.downconverted_signal_f32, 1) == b ?
            ds.downconverted_signal_f32 :
            StructArray{Complex{Float32}}((Matrix{Float32}(undef, b, num_ants), Matrix{Float32}(undef, b, num_ants))),
        ds.downconverted_signal_f64
    )
end

# No need for resizing when dealing with GPU signals
function resize!(ds::Nothing, b::Integer, signal::AbstractArray)
    return ds
end
# No need for resizing the GPU GNSS codes
function resize!(codes::Nothing, b::Integer)
    return codes
end