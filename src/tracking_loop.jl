# For now
function downconvert(x, replica)
    x .* conj(replica)
end

function correlate(x, replica)
    replica' * x
end


function init_tracking(init_PLL, init_DLL, init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, n_samples, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn)
    PLL, init_carrier_replica = init_PLL(init_carrier_phase, init_carrier_freq, n_samples,  f_s, pll_disc_bandwidth, Δt)
    DLL, init_code_replica = init_DLL(init_code_phase, init_code_freq, n_samples,  f_s, dll_disc_bandwidth, Δt, sat_prn)
    signal -> _tracking(signal, PLL, DLL, init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, n_samples, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn)
end

function _tracking(signal, PLL, DLL, carrier_replica, code_replica, beamform)
    downconverted_signal = downconvert(signal, carrier_replica)
    correlated_signal = correlate(downconverted_signal, code_replica)
    beamformed_signal = beamform(correlated_signal)
    next_PLL, next_carrier_replica = PLL(beamformed_signal)
    next_DLL, next_code_replicas, code_phase = DLL(beamformed_signal)
    next_signal -> _tracking(next_signal, next_PLL, next_DLL, next_carrier_replica, next_code_replicas, beamform), code_phase , prompt(correlated_signal), prompt(beamformed_signal) 
end
