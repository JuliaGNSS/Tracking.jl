function estimate_C╱N₀(corr_outs, corr_time_ms)
    inphase_power = mean(abs.(real.(corr_outs)))^2
    total_power = mean(abs2.(corr_outs))
    SNR = 10 * log10(inphase_power / (total_power - inphase_power))
    return SNR - 10 * log10(1e-3 * corr_time_ms)
end


function init_joined_tracking(
    l1_init_carrier_phase, l1_init_carrier_freq, l1_init_code_phase, l1_init_code_freq, l1_f_s, l1_beamform, l1_pll_disc_bandwidth, l1_dll_disc_bandwidth, l1_scale_factor,
    l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, l5_scale_factor,
    Δt, sat_prn)

    l1_beamformed_corr_outs = CircularBuffer{Complex{Float64}}(10)
    l5_beamformed_corr_outs = CircularBuffer{Complex{Float64}}(10)
    gpsl1_tracking_loop = init_tracking(l1_init_carrier_phase, l1_init_carrier_freq, l1_init_code_phase, l1_init_code_freq, Δt, l1_f_s, l1_beamform, l1_pll_disc_bandwidth, l1_dll_disc_bandwidth, sat_prn, l1_scale_factor)
    gpsl5_tracking_loop = init_tracking(l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, Δt, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, sat_prn, l5_scale_factor)
    (gpsl1_signal, gpsl5_signal, gpsl1_aiding = 0, gpsl5_aiding = 0) -> _joined_tracking(gpsl1_signal, gpsl5_signal, gpsl1_aiding, gpsl5_aiding, gpsl1_tracking_loop, gpsl5_tracking_loop, l1_beamform, l5_beamform, l1_beamformed_corr_outs, l5_beamformed_corr_outs)
end


function _joined_tracking(gpsl1_signal, gpsl5_signal, gpsl1_aiding, gpsl5_aiding, gpsl1_tracking_loop, gpsl5_tracking_loop, l1_beamform, l5_beamform, l1_beamformed_corr_outs, l5_beamformed_corr_outs)
    
    next_gpsl1_tracking_loop, gpsl1_code_phase, gpsl1_prompts_corr_signal, gpsl1_carrier_freq_update = gpsl1_tracking_loop(gpsl1_signal, gpsl1_aiding)
    next_gpsl5_tracking_loop, gpsl5_code_phase, gpsl5_prompts_corr_signal, gpsl5_carrier_freq_update = gpsl5_tracking_loop(gpsl5_signal, gpsl5_aiding)
    append!(l1_beamformed_corr_outs,l1_beamform(gpsl1_prompts_corr_signal))
    append!(l5_beamformed_corr_outs, l5_beamform(gpsl5_prompts_corr_signal))
    l1_SNR = estimate_C╱N₀(l1_beamformed_corr_outs, length(l1_beamformed_corr_outs))
    l5_SNR = estimate_C╱N₀(l5_beamformed_corr_outs, length(l5_beamformed_corr_outs))
    internal_gpsl1_aiding = 0
    internal_gpsl5_aiding = 0
#=     
    if l1_SNR < 10
        if l5_SNR > 10
            #internal_gpsl1_aiding = gpsl1_carrier_freq_update ToDO
        else
            #println("both signals lost")
        end
    else
        if l5_SNR < 10
            #internal_gpsl5_aiding = gpsl5_carrier_freq_update 
        end
    end  =#
    (next_gpsl1_signal, next_gpsl5_signal, next_gpsl1_aiding = 0, next_gpsl5_aiding = 0) -> _joined_tracking(next_gpsl1_signal, next_gpsl5_signal, next_gpsl1_aiding + internal_gpsl1_aiding, next_gpsl5_aiding + internal_gpsl5_aiding, next_gpsl1_tracking_loop, next_gpsl5_tracking_loop, l1_beamform, l5_beamform, l1_beamformed_corr_outs, l5_beamformed_corr_outs),
    gpsl1_code_phase, gpsl1_carrier_freq_update,
    gpsl5_code_phase, gpsl5_carrier_freq_update
end
