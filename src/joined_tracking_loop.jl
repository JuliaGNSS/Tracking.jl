"""
$(SIGNATURES)

Calculate the  mean Signal to Noise ratio, for the `corr_outs` in the time frame `corr_time_ms`.

"""

function estimate_C╱N₀(corr_outs, corr_time_ms)
    inphase_power = mean(abs.(real.(corr_outs)))^2
    total_power = mean(abs2.(corr_outs))
    SNR = 10 * log10(inphase_power / (total_power - inphase_power))
    return SNR - 10 * log10(1e-3 * corr_time_ms)
end

"""
$(SIGNATURES)

Initialize the joined_tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a joined_trackin_loop function to track 2 signals at the same time.

# Arguments
- `l1_init_carrier_phase::Float`: the initial replicated carrier signal phase in rad
- `l1_init_carrier_freq::Float`: the initial replicated carrier signal frequency in Hz
- `l1_init_carrier_phase::Float`: the initial replicated PRN code offset in chip
- `l1_init_carrier_freq::Float`: the initial replicated PRN code signal frequency in Hz
- `l1_f_s::Float`: the signal sampling frequency in Hz
- `l1_beamform::Function`: a beamforming function to produce a new early_prompt_late tuple from multiple inputs
- `l1_pll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth for the PLL discriminator, in Hz
- `l1_dll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth for the DLL discriminator, in Hz
- `l1_scale_factor::Float`: scale filter that is apllied to the carrier outputs before adding them to the code loop
-  all the above for l5
- `sat_prn::Integer`: the satellite PRN code number, choose from 1-32
- `Δt::Float`: the loop update time intervall in seconds


# Examples
```julia-repl
    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end
    l1_scale_factor = 1.023e6/1575.43e6
    l5_scale_factor = 1.023e6/1176.45e6
    velocity_aiding = 0.0
    joined_tracking_loop = init_joined_tracking(
        #     l5_init_carrier_phase, l5_init_carrier_freq, l5_init_code_phase, l5_init_code_freq, l5_f_s, l5_beamform, l5_pll_disc_bandwidth, l5_dll_disc_bandwidth, l5_scale_factor
        1/3 * π, 50, 2.0, 1023e3, 4e6, beamform, 18.0, 1.0, l1_scale_factor,
        2/3 * π, 50, 1.0, 1023e4, 4e7, beamform, 18.0, 1.0, l5_scale_factor,
        1e-3, 1 # Δt, sat_prn
    )
    l1_test_signal = cis.(2 * π * 50 / 4e6 * (1:32000) + 1 / 3 * π)
    l1_samples_code = GNSSSignals.gen_sat_code(1:32000, 1.023e6, 2.0, 4e6, SATELLITE_1_CODE) * 10^(-20/20)
    l1_incoming_signals = [1, 1] .* (l1_test_signal[1:4000] .* l1_samples_code[1:4000])'

    l5_test_signal = cis.(2 * π * 50 / 4e7 * (1:102300) + 2 / 3 * π)
    l5_samples_code = GNSSSignals.gen_sat_code(1:102300, 1.023e7, 1.0, 4e7, L5_SAT1_CODE) * 10^(-20/20)
    l5_incoming_signals = [1, 1] .* (l5_test_signal[1:40000] .* l5_samples_code[1:40000])'

    next_joined_tracking_loop, l1_code_phase, l1_carrier_freq_update, l5_code_phase, l5_carrier_freq_update = joined_tracking_loop(l1_incoming_signals, l5_incoming_signals)
```
    """

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

"""
$(SIGNATURES)

Should be initialized by init_joined_tracking, uses the provided 'PLL', 'DLL' and 'beamform' functions together with the provided antenna 'signals', the provided 'scale_factor, the 'velocity_aiding', and the replicated samples/codes 'carrier_replica' and 'code_replicas' to calculate the functions and samples/codes for the next timestep for both L1 and L5.
Returns the _tracking function for the next time step together with the the code_phase, the carrier_frequency_update, and the prompt of the correlated signals for L1 and L5.

"""

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
