"""
$(SIGNATURES)

A buffer that holds CPU buffers for necessary replicas and downconverted
signal.
"""
struct CPUSatDownconvertAndCorrelator{T,CT,DS<:StructVecOrMat{Complex{T}}}
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

const MultipleSystemSatCPUDownconvertAndCorrelator{N,I,T,CT,DS} =
    TupleLike{<:NTuple{N,Dictionary{I,<:CPUSatDownconvertAndCorrelator{<:T,<:CT,<:DS}}}}
struct CPUDownconvertAndCorrelator{
    MESF,
    N,
    I,
    B<:MultipleSystemSatCPUDownconvertAndCorrelator{N,I},
} <: AbstractDownconvertAndCorrelator{N,I}
    buffers::B
end

function CPUDownconvertAndCorrelator(
    maximum_expected_sampling_frequency::Maybe{Val},
    multiple_system_sats_state::MultipleSystemSatsState{N,I},
    num_samples::Int,
) where {N,I}
    if isnothing(maximum_expected_sampling_frequency)
        error(
            "The maximum expected sampling frequency needs to be specified, when downconverting and correlating on the CPU. You can either specify it when initializing the track state or directly when constructing the CPUDownconvertAndCorrelator",
        )
    end
    buffers = map(multiple_system_sats_state) do system_sats_state
        map(system_sats_state.states) do sat_state
            CPUSatDownconvertAndCorrelator(
                Float32,
                system_sats_state.system,
                sat_state.correlator,
                num_samples,
            )
        end
    end
    CPUDownconvertAndCorrelator(maximum_expected_sampling_frequency, buffers)
end

function CPUDownconvertAndCorrelator(
    maximum_expected_sampling_frequency::Val{MESF},
    buffers::B,
) where {MESF,N,I,B<:MultipleSystemSatCPUDownconvertAndCorrelator{N,I}}
    CPUDownconvertAndCorrelator{MESF,N,I,B}(buffers)
end

# See https://github.com/JuliaLang/julia/issues/43155
# Remove when https://github.com/JuliaLang/julia/pull/43338 is merged
function Base.setindex(nt::NamedTuple{names,T}, val::X, idx::Int) where {names,T,X}
    NamedTuple{names}(Base.setindex(Tuple(nt), val, idx))
end

function merge_sats(
    system::AbstractGNSS,
    downconvert_and_correlator::CPUDownconvertAndCorrelator{MESF,N,I,B},
    system_idx,
    sats_state::Dictionary{I,<:SatState},
    num_samples::Int,
) where {MESF,N,I,T,CT,DS,B<:MultipleSystemSatCPUDownconvertAndCorrelator{N,I,T,CT,DS}}
    new_buffers = map(sats_state) do sat_state
        CPUSatDownconvertAndCorrelator(T, system, sat_state.correlator, num_samples)
    end
    new_downconvert_and_correlator_buffers = setindex(
        downconvert_and_correlator.buffers,
        merge(downconvert_and_correlator.buffers[system_idx], new_buffers),
        system_idx,
    )
    CPUDownconvertAndCorrelator{MESF,N,I,B}(new_downconvert_and_correlator_buffers)
end

function filter_out_sats(
    downconvert_and_correlator::CPUDownconvertAndCorrelator{MESF,N,I,B},
    system_idx::Union{Symbol,Integer},
    identifiers,
) where {MESF,N,I,T,CT,DS,B<:MultipleSystemSatCPUDownconvertAndCorrelator{N,I,T,CT,DS}}
    filtered_buffers = map(
        last,
        filter(
            ((id,),) -> !in(id, identifiers),
            pairs(downconvert_and_correlator.buffers[system_idx]),
        ),
    )
    new_downconvert_and_correlator_buffers =
        setindex(downconvert_and_correlator.buffers, filtered_buffers, system_idx)
    CPUDownconvertAndCorrelator{MESF,N,I,B}(new_downconvert_and_correlator_buffers)
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
        <:MultipleSystemSatsState,
        <:AbstractDopplerEstimator,
        <:CPUDownconvertAndCorrelator{MESF},
    },
    preferred_num_code_blocks_to_integrate::Int,
    sampling_frequency,
    intermediate_frequency,
    num_samples_signal::Int,
) where {MESF}
    new_multiple_system_sats_state = map(
        track_state.multiple_system_sats_state,
        track_state.downconvert_and_correlator.buffers,
    ) do system_sats_state, buffers
        new_sat_states = map(system_sats_state.states, buffers) do sat_state, buffer
            signal_samples_to_integrate, is_integration_completed =
                calc_signal_samples_to_integrate(
                    system_sats_state.system,
                    sat_state.signal_start_sample,
                    sampling_frequency,
                    sat_state.code_doppler,
                    sat_state.code_phase,
                    preferred_num_code_blocks_to_integrate,
                    found(sat_state.sc_bit_detector),
                    num_samples_signal,
                )
            if signal_samples_to_integrate == 0
                return sat_state
            end
            carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
            code_frequency =
                sat_state.code_doppler + get_code_frequency(system_sats_state.system)

            new_correlator = downconvert_and_correlate!(
                system_sats_state.system,
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
                sat_state.signal_start_sample,
                signal_samples_to_integrate,
                sat_state.prn,
                Val{MESF}(),
            )::typeof(sat_state.correlator)
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
