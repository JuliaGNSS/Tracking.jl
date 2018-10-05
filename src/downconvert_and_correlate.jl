function downconvert_and_correlate!(signal, output, gen_carrier_replica, gen_code_replica, sample_shift)
    shifts = SVector(1, sample_shift + 1, 2 * sample_shift + 1) # Currently only early prompt and late
    code_sampled_buffer = CircularBuffer{Int}(2 * sample_shift + 1)
    append!(code_sampled_buffer, gen_code_replica.(1 - sample_shift:1 + sample_shift))
    for sample_idx = 1:size(signal, 1)
        carrier = gen_carrier_replica(i)
        for (output_idx, code_shift_idx) in enumerate(shifts)
            code = code_sampled_buffer[code_shift_idx]
            for ant_idx = 1:size(signal, 2)
                @inbounds output[output_idx][ant_idx] += signal[sample_idx,ant_idx] * code * carrier
            end
        end
        push!(code_sampled_buffer, gen_code_replica(1 + sample_shift + sample_idx))
    end
    output
end

function calc_sample_shift(code_freq, sample_freq, preferred_phase)
    sample_shift = round(preferred_phase * sample_freq / code_freq)
    actual_phase = sample_shift * code_freq / sample_freq
    Int(sample_shift), actual_phase
end
