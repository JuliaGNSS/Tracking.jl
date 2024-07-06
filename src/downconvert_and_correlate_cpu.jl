struct CPUSystemDownconvertAndCorrelator <: AbstractSystemDownconvertAndCorrelator end

"""
$(SIGNATURES)

A buffer that holds CPU buffers for necessary replicas and downconverted
signal.
"""
struct CPUSatDownconvertAndCorrelator{T,CT,DS<:StructVecOrMat{Complex{T}}} <:
       AbstractSatDownconvertAndCorrelator
    code_replica_buffer::Vector{CT}
    carrier_replica_buffer::StructVector{Complex{T}}
    downconvert_signal_buffer::DS
end

"""
$(SIGNATURES)

Convenient constructor to initialize buffers for the CPU with the correct lengths for a single
satellite.
"""
function CPUSatDownconvertAndCorrelator(
    ::Type{T},
    system::AbstractGNSS,
    correlator::AbstractCorrelator,
    num_samples,
) where {T}
    code_shifts = get_shifts(correlator)
    CPUSatDownconvertAndCorrelator(
        Vector{get_code_type(system)}(
            undef,
            num_samples + maximum(code_shifts) - minimum(code_shifts),
        ),
        StructVector{Complex{T}}(undef, num_samples),
        get_downconvert_signal_buffer(T, num_samples, correlator),
    )
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

Convenient constructor to initialize buffers for the CPU with the correct lengths for a single
satellite. This constructor uses Float32 as the sample data type.
"""
function CPUSatDownconvertAndCorrelator(
    system::AbstractGNSS,
    correlator::AbstractCorrelator,
    num_samples,
)
    CPUSatDownconvertAndCorrelator(Float32, system, correlator, num_samples)
end

"""
$(SIGNATURES)

Downconvert und correlate all available satellites on the CPU.
"""
function downconvert_and_correlate(
    signal,
    track_state::TrackState{
        <:MultipleSystemSatsState{
            N,
            <:AbstractGNSS,
            <:SatState{
                <:AbstractCorrelator,
                <:AbstractSatDopplerEstimator,
                <:CPUSatDownconvertAndCorrelator,
            },
            <:AbstractSystemDopplerEstimator,
            <:CPUSystemDownconvertAndCorrelator,
        },
    },
    sample_params::TupleLike{<:NTuple{N,Dictionary{I,SampleParams}}},
    sampling_frequency,
    intermediate_frequency,
    num_samples_signal::Int,
) where {I,N}
    new_multiple_system_sats_state = map(
        track_state.multiple_system_sats_state,
        sample_params,
    ) do system_sats_state, system_sat_params
        new_sat_states =
            map(system_sat_params, system_sats_state.states) do sat_params, sat_state
                if sat_params.signal_samples_to_integrate == 0
                    return sat_state
                end
                carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
                code_frequency =
                    sat_state.code_doppler + get_code_frequency(system_sats_state.system)

                new_correlator = downconvert_and_correlate!(
                    system_sats_state.system,
                    signal,
                    sat_state.correlator,
                    sat_state.downconvert_and_correlator.code_replica_buffer,
                    sat_state.code_phase,
                    sat_state.downconvert_and_correlator.carrier_replica_buffer,
                    sat_state.carrier_phase,
                    sat_state.downconvert_and_correlator.downconvert_signal_buffer,
                    code_frequency,
                    carrier_frequency,
                    sampling_frequency,
                    sat_params.signal_start_sample,
                    sat_params.signal_samples_to_integrate,
                    sat_state.prn,
                )::typeof(sat_state.correlator)
                return update(
                    system_sats_state.system,
                    sat_state,
                    sat_params,
                    intermediate_frequency,
                    sampling_frequency,
                    new_correlator,
                    num_samples_signal,
                )
            end
        return SystemSatsState(system_sats_state, new_sat_states)
    end
    return TrackState(track_state, new_multiple_system_sats_state)
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
