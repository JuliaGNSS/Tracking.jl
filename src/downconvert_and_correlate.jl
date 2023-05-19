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
        correlator.shifts,
        prn
    )
    gen_carrier_replica!(
        carrier_replica,
        carrier_frequency,
        sampling_frequency,
        carrier_phase,
        signal_start_sample,
        num_samples_left
    )
    downconvert!(
        downconverted_signal,
        signal,
        carrier_replica,
        signal_start_sample,
        num_samples_left
    )
    correlate(
        correlator,
        downconverted_signal,
        code_replica,
        signal_start_sample,
        num_samples_left
    )
end

abstract type AbstractDownconvertAndCorrelator end

const StructVecOrMat{T} = Union{StructVector{T}, StructArray{T, 2}}

struct CPUBuffers{T, CT, DS <: StructVecOrMat{Complex{T}}}
    code_replica_buffer::Vector{CT}
    carrier_replica_buffer::StructVector{Complex{T}}
    downconvert_signal_buffer::DS
end

function CPUBuffers(
    ::Type{T},
    system::AbstractGNSS,
    state::SatState,
    num_samples,
) where T
    code_shifts = get_correlator(state).shifts
    CPUBuffers(
        Vector{get_code_type(system)}(undef, num_samples + maximum(code_shifts) - minimum(code_shifts)),
        StructVector{Complex{T}}(undef, num_samples),
        get_downconvert_signal_buffer(T, num_samples, state)
    )
end

function CPUBuffers(
    system::AbstractGNSS,
    state::SatState,
    num_samples,
)
    CPUBuffers(Float32, system, state, num_samples)
end

struct CPUDownconvertAndCorrelator{N, B <: TupleLike{<:NTuple{N, Vector{<:CPUBuffers}}}} <: AbstractDownconvertAndCorrelator
    buffers::B
    function CPUDownconvertAndCorrelator(buffers::TupleLike{<:NTuple{N, Vector{<:CPUBuffers}}}) where N
        new{N, typeof(buffers)}(buffers)
    end
end

function CPUDownconvertAndCorrelator(
    ::Type{T},
    system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}},
    num_samples::Integer,
) where {T, N}
    CPUDownconvertAndCorrelator(
        map(system_sats_states) do system_sats_state
            system = system_sats_state.system
            map(system_sats_state.states) do state
                CPUBuffers(T, system, state, num_samples)
            end
        end
    )
end
function CPUDownconvertAndCorrelator(
    system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}},
    num_samples::Integer,
) where N
    CPUDownconvertAndCorrelator(Float32, system_sats_states, num_samples)
end

function get_downconvert_signal_buffer(::Type{T}, num_samples::Int, state::SatState{<: AbstractCorrelator{1}}) where T
    StructVector{Complex{T}}(undef, num_samples)
end
function get_downconvert_signal_buffer(::Type{T}, num_samples::Int, state::SatState{<: AbstractCorrelator{M}}) where {T, M}
    StructArray{Complex{T}}(undef, num_samples, M)
end

function downconvert_and_correlate(
    downconvert_and_correlator::CPUDownconvertAndCorrelator,
    signal,
    sampling_frequency,
    intermediate_frequency,
    system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}},
    params::TupleLike{<:NTuple{N, Vector{SampleParams}}},
) where N
    map(params, system_sats_states, downconvert_and_correlator.buffers) do system_params, system_sats, buffers
        map(system_params, system_sats.states, buffers) do sat_params, sat_state, buffer
            if sat_params.signal_samples_to_integrate == 0
                return sat_state.correlator
            end
            carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
            code_frequency = sat_state.code_doppler + get_code_frequency(system_sats.system)
            downconvert_and_correlate!(
                system_sats.system,
                signal,
                sat_state.correlator,
                buffer.code_replica_buffer,
                sat_state.code_phase,
                buffer.carrier_replica_buffer,
                sat_state.carrier_phase,
                buffer.downconvert_signal_buffer,
                code_frequency,
                carrier_frequency,
                sampling_frequency,
                sat_params.signal_start_sample,
                sat_params.signal_samples_to_integrate,
                sat_state.prn
            )::typeof(sat_state.correlator)
        end
    end
end