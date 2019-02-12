struct CN0Buffer
    summed_total_power::Float64
    summed_abs_inphase_ampl::Float64
    last_valid_cn0::Float64
    update_time::typeof(1.0s)
    current_updated_time::typeof(1.0s)
end

function CN0Buffer(update_time)
    CN0Buffer(0.0, 0.0, NaN, update_time, 0.0s)
end

function estimate_CN0(cn0_buffer::CN0Buffer, integration_time, prompt_correlator_output)
    summed_total_power = cn0_buffer.summed_total_power + abs2(prompt_correlator_output)
    summed_abs_inphase_ampl = cn0_buffer.summed_abs_inphase_ampl + abs(real(prompt_correlator_output))
    current_updated_time = cn0_buffer.current_updated_time + integration_time
    last_valid_cn0 = cn0_buffer.last_valid_cn0
    if current_updated_time >= cn0_buffer.update_time
        inphase_power = (summed_abs_inphase_ampl / round(upreferred(current_updated_time / integration_time)))^2
        total_power = summed_total_power / round(upreferred(current_updated_time / integration_time))
        last_valid_cn0 = inphase_power / (total_power - inphase_power) / upreferred(integration_time / s)
        summed_total_power = 0.0
        summed_abs_inphase_ampl = 0.0
        current_updated_time = 0.0s
    end
    CN0Buffer(summed_total_power, summed_abs_inphase_ampl, last_valid_cn0, cn0_buffer.update_time, current_updated_time)
end
