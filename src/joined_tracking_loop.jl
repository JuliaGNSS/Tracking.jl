function estimate_C╱N₀(corr_outs, corr_time_ms)
    inphase_power = mean(abs.(real.(corr_outs)))^2
    total_power = mean(abs2.(corr_outs))
    SNR = 10 * log10(inphase_power / (total_power - inphase_power))
    return SNR - 10 * log10(1e-3 * corr_time_ms)
end


function init_joined_tracking(
    l1_init_carrier_phase, l1_init_carrier_freq, l1_init_code_phase, l1_init_code_freq, l1_f_s, l1_beamform, l1_pll_disc_bandwidth, l1_dll_disc_bandwidth, l1_scale_factor,
    l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, l5_scale_factor,
    Δt, sat_prn
)
    gpsl1_tracking_loop = init_tracking(l1_init_carrier_phase, l1_init_carrier_freq, l1_init_code_phase, l1_init_code_freq, Δt, l1_f_s, l1_beamform, l1_pll_disc_bandwidth, l1_dll_disc_bandwidth, sat_prn, l1_scale_factor)
    gpsl5_tracking_loop = init_tracking(l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, Δt, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, sat_prn, l5_scale_factor)
    (gpsl1_signal, gpsl1_aiding, gpsl5_signal, gpsl5_aiding) -> _joined_tracking(gpsl1_signal, gpsl1_aiding, gpsl5_signal, gpsl5_aiding, gpsl1_tracking_loop, gpsl5_tracking_loop)
end


function _joined_tracking(gpsl1_signal, gpsl1_aiding, gpsl5_signal, gpsl5_aiding, gpsl1_tracking_loop, gpsl5_tracking_loop)

#=     if gpsl1_trouble()
        if !gpsl5_trouble()
            get_aiding and values from gpsl5
        else
            println("we lost the signals")
        end
    else
        if gpsl5_trouble()
            get aiding and values from gpsl1_aid
        end
    end =#
    next_gpsl1_tracking_loop, gpsl1_code_phase, gpsl1_carrier_freq_update = gpsl1_tracking_loop(gpsl1_signal, gpsl1_aiding)
    next_gpsl1_tracking_loop, gpsl1_code_phase, gpsl1_carrier_freq_update = gpsl1_tracking_loop(gpsl5_signal, gpsl5_aiding)
    
    (next_gpsl1_signal, next_gpsl1_aiding, next_gpsl5_signal, next_gpsl5_aiding) -> _joined_tracking(next_gpsl1_signal, next_gpsl1_aiding, next_gpsl5_signal, next_gpsl5_aiding, next_gpsl1_tracking_loop, next_gpsl5_tracking_loop)
end