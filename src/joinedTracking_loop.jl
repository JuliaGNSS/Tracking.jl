function estimate_C╱N₀(corr_outs, corr_time_ms)
    inphase_power = mean(abs.(real.(corr_outs)))^2
    total_power = mean(abs2.(corr_outs))
    SNR = 10 * log10(inphase_power / (total_power - inphase_power))
    return SNR - 10 * log10(1e-3 * corr_time_ms)
end


function init_joined_tracking(
    l1_init_carrier_phase, l1_init_carrier_freq, l1_init_code_phase, l1_init_code_freq, Δt, l1_f_s, beamform, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn, scale_factor
    init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn, scale_factor
)
    gpsl1_tracking_loop = init_tracking()
    gpsl5_tracking_loop = init_tracking()
    (gpsl1_signal, gpsl1_aiding, gpsl5_signal, gpsl5_aiding) -> _joined_tracking()
end


function _joined_tracking()
    next_gpsl1_tracking_loop, gpsl1_code_phase, gpsl1_carrier_freq_update = gpsl1_tracking_loop()
    next_gpsl1_tracking_loop, gpsl1_code_phase, gpsl1_carrier_freq_update = gpsl1_tracking_loop()
    if gpsl1_trouble()
        if !gpsl5_trouble()
            get_aiding and values from gpsl5
        else
            println("we lost the signals")
        end
    else
        if gpsl5_trouble()
            get aiding and values from gpsl1_aid
        end
    end
    
    (gpsl1_signal, gpsl1_aiding, gpsl5_signal, gpsl5_aiding) -> _joined_tracking()
end