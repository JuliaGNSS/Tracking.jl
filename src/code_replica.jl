"""
$(SIGNATURES)

Generate a code replica for a signal from satellite system `S`. The
replica contains `num_samples` prompt samples as well as an additional
number of early and late samples specified by `correlator_sample_shifts`.
The codefrequency is specified by `code_frequency`, while the sampling
rate is given by `sampling_frequency`. The phase of the first prompt
sample is given by `start_code_phase`. The generated signal is returned
in the array `code_replica` with the first generated sample written to
index start_sample.
"""
function gen_code_replica!(
    code_replica,
    system::AbstractGNSS,
    code_frequency,
    sampling_frequency,
    start_code_phase::AbstractFloat,
    start_sample::Integer,
    num_samples::Integer,
    correlator_sample_shifts::AbstractVector,
    prn::Integer
)
    earliest_sample_shift = correlator_sample_shifts[end]
    latest_sample_shift  = correlator_sample_shifts[1]
    total_samples = num_samples + earliest_sample_shift - latest_sample_shift
    gen_code!(
        view(code_replica, start_sample:start_sample + total_samples - 1),
        system,
        prn,
        sampling_frequency,
        code_frequency,
        start_code_phase,
        latest_sample_shift
    )
    code_replica
end

"""
$(SIGNATURES)

Updates the code phase.
"""
function update_code_phase(
    system::AbstractGNSS,
    num_samples,
    code_frequency,
    sampling_frequency,
    start_code_phase,
    secondary_code_or_bit_found
)
    secondary_code_or_bit_length = get_data_frequency(system) == 0Hz ?
        get_secondary_code_length(system) :
        Int(get_code_frequency(system) / (get_data_frequency(system) * get_code_length(system)))

    code_length = get_code_length(system) *
        (secondary_code_or_bit_found ? secondary_code_or_bit_length : 1)
    mod(code_frequency * num_samples / sampling_frequency + start_code_phase, code_length)
#    fixed_point = sizeof(Int) * 8 - 1 - min_bits_for_code_length(S)
#    delta = floor(Int, code_frequency * 1 << fixed_point / sampling_frequency)
#    fixed_point_start_phase = floor(Int, start_code_phase * 1 << fixed_point)
#    phase_fixed_point = delta * num_samples + fixed_point_start_phase
#    mod(phase_fixed_point / 1 << fixed_point, code_length)
end

"""
$(SIGNATURES)

Calculates the current code frequency.
"""
function get_current_code_frequency(system::AbstractGNSS, code_doppler)
    code_doppler + get_code_frequency(system)
end
