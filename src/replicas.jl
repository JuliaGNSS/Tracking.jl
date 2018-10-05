"""
$(SIGNATURES)

Initializes a general replica generator given by the phase `phase`, sample frequency `sampling_freq`
and the two functions `calc_signal` and `calc_phase` to generate the signal and calculate the phase respectively.
A function is returned which depends on the number of signals of the replica `num_samples` and the current
frequency `freq`.
"""
function init_replica(init_freq, phase, sampling_freq, calc_signal, calc_phase)
    (num_samples, freq_update, phase_shift = 0.0) -> _replica(num_samples, phase + phase_shift, init_freq, freq_update, sampling_freq, calc_signal, calc_phase)
end

"""
$(SIGNATURES)

Generates the replica with `num_samples` number of samples depending on the phase `phase`, frequency `freq`,
sampling frequency `sampling_freq` and the two functions `calc_signal` and `calc_phase` to generate
the signal and calculate the phase respectively.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function _replica(num_samples, phase, init_freq, freq_update, sample_freq, calc_signal, calc_phase)
    replica = calc_signal(1:num_samples, init_freq + freq_update, phase, sample_freq)
    next_phase = calc_phase(num_samples, init_freq + freq_update, phase, sample_freq)
    (next_num_samples, next_freq_update, phase_shift = 0.0) -> _replica(next_num_samples, next_phase + phase_shift, init_freq, next_freq_update, sample_freq, calc_signal, calc_phase), replica, next_phase
end

"""
$(SIGNATURES)

Initializes the code replica generator based on the initial code phase `init_code_phase`, sampling frequency
`sampling_freq`, satellite pseudo random noise number `sat_prn`.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function init_code_replica(system, init_freq, init_code_phase, sample_freq, sat_prn)
    early_prompt_late_phase = [-0.5, 0.0, 0.5]
    gen_replica_code(samples, freq, phase, used_sample_freq) = begin
      sample_shifts = round.(Int, early_prompt_late_phase .* used_sample_freq / freq)
      sampled_code = gen_code(system, minimum(samples) + sample_shifts[1]:maximum(samples) + sample_shifts[3], freq, phase, used_sample_freq, sat_prn)
      map(sample_shift -> sampled_code[(minimum(samples) - sample_shifts[1] + sample_shift):(maximum(samples) - sample_shifts[1] + sample_shift)], sample_shifts)
    end
    calc_replica_phase(sample, freq, phase, used_sample_freq) = calc_code_phase(sample, freq, phase, used_sample_freq, system.code_length)
    init_replica(init_freq, init_code_phase, sample_freq, gen_replica_code, calc_replica_phase)
end


"""
$(SIGNATURES)

Initializes the carrier replica generator based on the initial carrier phase `init_carrier_phase`, sampling frequency
`sampling_freq`.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function init_carrier_replica(init_freq, init_carrier_phase, sample_freq)
    init_replica(init_freq, init_carrier_phase, sample_freq, gen_carrier, get_carrier_phase)
end
