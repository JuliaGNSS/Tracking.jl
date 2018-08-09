"""
$(SIGNATURES)

Initializes a general replica generator given by the phase `phase`, sample frequency `sampling_freq`
and the two functions `calc_signal` and `calc_phase` to generate the signal and calculate the phase respectively.
A function is returned which depends on the number of signals of the replica `num_samples` and the current
frequency `freq`.
"""
function init_replica(phase, sampling_freq, calc_signal, calc_phase)
    (num_samples, freq) -> _replica(num_samples, phase, freq, sampling_freq, calc_signal, calc_phase)
end

"""
$(SIGNATURES)

Generates the replica with `num_samples` number of samples depending on the phase `phase`, frequency `freq`,
sampling frequency `sampling_freq` and the two functions `calc_signal` and `calc_phase` to generate
the signal and calculate the phase respectively.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function _replica(num_samples, phase, freq, sampling_freq, calc_signal, calc_phase)
    replica = calc_signal(1:num_samples, freq, phase, sampling_freq)
    next_phase = calc_phase(num_samples, freq, phase, sampling_freq)
    (num_samples, freq) -> _replica(num_samples, next_phase, freq, sampling_freq, calc_signal, calc_phase), replica, next_phase
end

"""
$(SIGNATURES)

Initializes the code replica generator based on the initial code phase `init_code_phase`, sampling frequency
`sampling_freq`, satellite space vehicle ID `svid`, and the two functions `gen_sampled_code` and `get_code_phase`
to generate the code signal and calculate the next code phase respectively. The code generators could refer to
any system e.g. L1, L5, etc.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function init_code_replica(init_code_phase, sampling_freq, sat_svid, gen_sampled_code, get_code_phase)
    early_prompt_late_phase = [-0.5, 0.0, 0.5]
    calc_signal(samples, freq, phase, sampling_freq) = begin 
      sample_shifts = round.(Int, early_prompt_late_phase .* sampling_freq / freq)
      sampled_code = gen_sampled_code(minimum(samples) + sample_shifts[1]:maximum(samples) + sample_shifts[3], freq, phase, sampling_freq, sat_svid)
      map(sample_shift -> sampled_code[(minimum(samples) - sample_shifts[1] + sample_shift):(maximum(samples) - sample_shifts[1] + sample_shift)], sample_shifts)
    end
    init_replica(init_code_phase, sampling_freq, calc_signal, get_code_phase)
end


"""
$(SIGNATURES)

Initializes the carrier replica generator based on the initial carrier phase `init_carrier_phase`, sampling frequency
`sampling_freq`.
A function is returned to calculate the replica with `num_samples` number of samples based on the current
frequency `freq`.
"""
function init_carrier_replica(init_carrier_phase, sampling_freq)
    init_replica(init_carrier_phase, sampling_freq, gen_carrier, get_carrier_phase)
end