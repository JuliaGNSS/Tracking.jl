"""
$(SIGNATURES)

Write the code replica of `signal` into `code_replica` over `sample_range`, widened
by the early and late `correlator_sample_offsets`, with the first prompt sample at
`start_code_phase`.
"""
function gen_code_replica!(
    code_replica,
    signal::AbstractGNSSSignal,
    code_frequency,
    sampling_frequency,
    start_code_phase::AbstractFloat,
    sample_range::AbstractUnitRange{<:Integer},
    correlator_sample_offsets::AbstractVector,
    prn::Integer,
)
    start_index = first(sample_range)
    num_samples = length(sample_range)
    earliest_sample_offset = maximum(correlator_sample_offsets)
    latest_sample_offset = minimum(correlator_sample_offsets)
    total_samples = num_samples + earliest_sample_offset - latest_sample_offset
    gen_code!(
        view(code_replica, start_index:(start_index+total_samples-1)),
        signal,
        prn,
        sampling_frequency,
        code_frequency,
        start_code_phase,
        latest_sample_offset,
    )
    code_replica
end
