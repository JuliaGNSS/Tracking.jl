struct CPUDownconvertAndCorrelator{MESF,B} <: AbstractDownconvertAndCorrelator
    buffer::B
end

function CPUDownconvertAndCorrelator(
    maximum_expected_sampling_frequency::Val{MESF},
    buffer::B = default_buffer(),
) where {MESF,B}
    CPUDownconvertAndCorrelator{MESF,B}(buffer)
end

function get_downconvert_signal_buffer(
    ::Type{T},
    num_samples::Int,
    correlator::AbstractCorrelator{1},
) where {T}
    StructVector{Complex{T}}(undef, num_samples)
end
function get_downconvert_signal_buffer(
    ::Type{T},
    num_samples::Int,
    correlator::AbstractCorrelator{M},
) where {T,M}
    StructArray{Complex{T}}(undef, num_samples, M)
end

"""
$(SIGNATURES)

Downconvert und correlate all available satellites on the CPU.
"""
function downconvert_and_correlate(
    downconvert_and_correlator::CPUDownconvertAndCorrelator{MESF},
    signal,
    track_state::TrackState,
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
) where {MESF}
    num_samples_signal = get_num_samples(signal)
    new_multiple_system_sats_state =
        map(track_state.multiple_system_sats_state) do system_sats_state
            new_sat_states = map(system_sats_state.states) do sat_state
                signal_samples_to_integrate, is_integration_completed =
                    calc_signal_samples_to_integrate(
                        system_sats_state.system,
                        sat_state.signal_start_sample,
                        sampling_frequency,
                        sat_state.code_doppler,
                        sat_state.code_phase,
                        preferred_num_code_blocks_to_integrate,
                        has_bit_or_secondary_code_been_found(sat_state),
                        num_samples_signal,
                    )
                if signal_samples_to_integrate == 0
                    return sat_state
                end
                carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
                code_frequency =
                    sat_state.code_doppler + get_code_frequency(system_sats_state.system)

                sample_shifts = get_correlator_sample_shifts(
                    sat_state.correlator,
                    sampling_frequency,
                    code_frequency,
                )
                @no_escape downconvert_and_correlator.buffer begin
                    code_replica_buffer = @alloc(
                        get_code_type(system_sats_state.system),
                        num_samples_signal + maximum(sample_shifts) -
                        minimum(sample_shifts)
                    )
                    carrier_replica_buffer_re = @alloc(Float32, num_samples_signal)
                    carrier_replica_buffer_im = @alloc(Float32, num_samples_signal)
                    carrier_replica_buffer = StructArray{ComplexF32}((
                        carrier_replica_buffer_re,
                        carrier_replica_buffer_im,
                    ))
                    downconvert_signal_buffer_re = @alloc(Float32, size(signal)...)
                    downconvert_signal_buffer_im = @alloc(Float32, size(signal)...)
                    downconvert_signal_buffer = StructArray{ComplexF32}((
                        downconvert_signal_buffer_re,
                        downconvert_signal_buffer_im,
                    ))
                    new_correlator = downconvert_and_correlate!(
                        system_sats_state.system,
                        signal,
                        sat_state.correlator,
                        code_replica_buffer,
                        sat_state.code_phase,
                        carrier_replica_buffer,
                        sat_state.carrier_phase,
                        downconvert_signal_buffer,
                        code_frequency,
                        carrier_frequency,
                        sampling_frequency,
                        sat_state.signal_start_sample,
                        signal_samples_to_integrate,
                        sat_state.prn,
                        Val{MESF}(),
                    )::typeof(sat_state.correlator)
                end
                return update(
                    system_sats_state.system,
                    sat_state,
                    signal_samples_to_integrate,
                    intermediate_frequency,
                    sampling_frequency,
                    new_correlator,
                    is_integration_completed,
                )
            end
            return SystemSatsState(system_sats_state, new_sat_states)
        end
    return TrackState(
        track_state;
        multiple_system_sats_state = new_multiple_system_sats_state,
    )
end

"""
$(SIGNATURES)

Downconvert and correlate a single satellite on the CPU.
"""
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
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn,
    maximum_expected_sampling_frequency,
)
    sample_shifts =
        get_correlator_sample_shifts(correlator, sampling_frequency, code_frequency)
    gen_code_replica!(
        code_replica,
        system,
        code_frequency,
        sampling_frequency,
        code_phase,
        signal_start_sample,
        num_samples_left,
        sample_shifts,
        prn,
        maximum_expected_sampling_frequency,
    )
    gen_carrier_replica!(
        carrier_replica,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        num_samples_left,
    )
    downconvert!(
        downconverted_signal,
        signal,
        carrier_replica,
        signal_start_sample,
        num_samples_left,
    )
    correlate(
        correlator,
        downconverted_signal,
        sample_shifts,
        code_replica,
        signal_start_sample,
        num_samples_left,
    )
end

#=
# This is currently slower than splitting the loop.
# See https://github.com/JuliaSIMD/LoopVectorization.jl/issues/284
function downconvert_and_correlate(
    signal::StructArray{Complex{T}},
    correlator::C,
    code,
    correlator_sample_shifts,
    carrier_frequency,
    sampling_frequency,
    start_phase,
    start_sample,
    num_samples
) where {T, C <: AbstractCorrelator}
    s_re = signal.re; s_im = signal.im
    accumulators = zero_accumulators(get_accumulators(correlator), signal)
    a_re = real.(accumulators)
    a_im = imag.(accumulators)
    @avx for i = start_sample:start_sample + num_samples - 1
        c_im, c_re = sincos(T(2π) * ((i - start_sample) * T(upreferred(carrier_frequency / Hz)) / T(upreferred(sampling_frequency / Hz)) + T(start_phase)))
        d_re = s_re[i] * c_re + s_im[i] * c_im
        d_im = s_im[i] * c_re - s_re[i] * c_im
        for j = 1:length(a_re)
            sample_shift = correlator_sample_shifts[j] - correlator_sample_shifts[1]
            a_re[j] += d_re * code[i + sample_shift]
            a_im[j] += d_im * code[i + sample_shift]
        end
    end
    accumulators_result = complex.(a_re, a_im)
    C(map(+, get_accumulators(correlator), accumulators_result))
end
=#
