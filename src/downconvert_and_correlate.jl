"""
$(SIGNATURES)

Abstract downconverter and correlator type. Structs for
downconversion and correlation must have this abstract type as a
parent.
"""
abstract type AbstractDownconvertAndCorrelator end

function get_downconvert_signal_buffer(
    ::Type{T},
    num_samples::Int,
    state::SatState{<:AbstractCorrelator{1}},
) where {T}
    StructVector{Complex{T}}(undef, num_samples)
end
function get_downconvert_signal_buffer(
    ::Type{T},
    num_samples::Int,
    state::SatState{<:AbstractCorrelator{M}},
) where {T,M}
    StructArray{Complex{T}}(undef, num_samples, M)
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
        prn,
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
        code_replica,
        signal_start_sample,
        num_samples_left,
    )
end
