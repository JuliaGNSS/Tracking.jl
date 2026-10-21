"""
$(SIGNATURES)

Write the code replica of `signal` into `code_replica` from index `start_sample`:
`num_samples` prompt samples, widened by the early and late
`correlator_sample_shifts`, with the first prompt sample at `start_code_phase`.
"""
function gen_code_replica!(
    code_replica,
    signal::AbstractGNSSSignal,
    code_frequency,
    sampling_frequency,
    start_code_phase::AbstractFloat,
    start_sample::Integer,
    num_samples::Integer,
    correlator_sample_shifts::AbstractVector,
    prn::Integer,
)
    earliest_sample_shift = maximum(correlator_sample_shifts)
    latest_sample_shift = minimum(correlator_sample_shifts)
    total_samples = num_samples + earliest_sample_shift - latest_sample_shift
    gen_code!(
        view(code_replica, start_sample:(start_sample+total_samples-1)),
        signal,
        prn,
        sampling_frequency,
        code_frequency,
        start_code_phase,
        latest_sample_shift,
    )
    code_replica
end

"""
$(SIGNATURES)

Updates the code phase.
"""
function update_code_phase(
    signal::AbstractGNSSSignal,
    num_samples,
    code_frequency,
    sampling_frequency,
    start_code_phase,
    secondary_code_or_bit_found,
)
    secondary_code_or_bit_length =
        get_data_frequency(signal) == 0Hz ? get_secondary_code_length(signal) :
        Int(
            get_code_frequency(signal) /
            (get_data_frequency(signal) * get_code_length(signal)),
        )

    code_length =
        get_code_length(signal) *
        (secondary_code_or_bit_found ? secondary_code_or_bit_length : 1)
    mod(code_frequency * num_samples / sampling_frequency + start_code_phase, code_length)
end

"""
$(SIGNATURES)

Calculates the current code frequency.
"""
function get_current_code_frequency(signal::AbstractGNSSSignal, code_doppler)
    code_doppler + get_code_frequency(signal)
end
