function downconvert_and_correlate!(signal, output, start_sample, num_samples_to_integrate, gen_carrier_replica, gen_code_replica, sample_shift)
    shifts = SVector(1, sample_shift + 1, 2 * sample_shift + 1) # Currently only early prompt and late
    code_sampled_buffer = CircularBuffer{Int}(2 * sample_shift + 1)
    append!(code_sampled_buffer, gen_code_replica.(1 - sample_shift:1 + sample_shift))
    for sample_idx = 1:num_samples_to_integrate
        carrier = gen_carrier_replica(sample_idx)
        for (output_idx, code_shift_idx) in enumerate(shifts)
            code = code_sampled_buffer[code_shift_idx]
            for ant_idx = 1:size(signal, 2)
                @inbounds output[output_idx][ant_idx] += signal[sample_idx + start_sample - 1,ant_idx] * code * carrier
            end
        end
        push!(code_sampled_buffer, gen_code_replica(1 + sample_shift + sample_idx))
    end
    output
end

function gen_replica_downconvert_correlate!(corr_res, system, signal, total_integration_time, integrated_samples, start_sample, sample_freq, interm_freq, code_sample_shift, carrier_doppler, code_doppler, sat_prn, velocity_aiding)
    num_samples = size(signal, 1)
    total_integration_samples = floor(Int, total_integration_time * sample_freq)
    num_samples_to_integrate = min(total_integration_samples - integrated_samples, num_samples - (start_sample - 1))
    gen_carrier_replica, gen_code_replica = create_gen_replicas(system, interm_freq, sample_freq, carrier_doppler, code_doppler, corr_res.carrier_phase, corr_res.code_phase, sat_prn)
    correlated_signals = downconvert_and_correlate!(signal, corr_res.outputs, start_sample, num_samples_to_integrate, gen_carrier_replica, gen_code_replica, code_sample_shift)
    next_carrier_phase, next_code_phase = calc_phases(system, interm_freq, sample_freq, carrier_doppler, code_doppler, corr_res.carrier_phase, corr_res.code_phase, num_samples_to_integrate)
    next_start_sample = start_sample + num_samples_to_integrate
    next_integrated_samples = integrated_samples + num_samples_to_integrate
    CorrelatorResults(next_carrier_phase, next_code_phase, correlated_signals), next_start_sample, next_integrated_samples, total_integration_samples
end

function create_gen_replicas(system, interm_freq, sample_freq, carrier_doppler, code_doppler, carrier_phase, code_phase, sat_prn)
    gen_carrier_replica(x) = gen_carrier(x, -(interm_freq + carrier_doppler), -carrier_phase, sample_freq)
    gen_code_replica(x) = gen_code(system, x, system.code_freq + code_doppler, code_phase, sample_freq, sat_prn)
    gen_carrier_replica, gen_code_replica
end

function calc_phases(system, interm_freq, sample_freq, carrier_doppler, code_doppler, carrier_phase, code_phase, num_samples_to_integrate)
    next_carrier_phase = calc_carrier_phase(num_samples_to_integrate, interm_freq + carrier_doppler, carrier_phase, sample_freq)
    next_code_phase = calc_code_phase(num_samples_to_integrate, system.code_freq + code_doppler, code_phase, sample_freq, size(system.codes, 1))
    next_carrier_phase, next_code_phase
end

function calc_sample_shift(system, sample_freq, preferred_phase)
    sample_shift = round(preferred_phase * sample_freq / system.code_freq)
    actual_phase = sample_shift * system.code_freq / sample_freq
    Int(sample_shift), actual_phase
end
