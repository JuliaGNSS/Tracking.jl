function downconvert(x, replica)
    map(signal -> signal .* conj(replica), x)
end

function correlate(x, replicas)
    map(replica -> map(single_antenna_signal -> replica' * single_antenna_signal, x), replicas)
end


function init_tracking(init_PLL, init_DLL, init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, num_samples, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn)
    PLL, init_carrier_replica = init_PLL(init_carrier_phase, init_carrier_freq, num_samples,  f_s, pll_disc_bandwidth, Δt)
    DLL, init_code_replicas = init_DLL(init_code_phase, init_code_freq, num_samples,  f_s, dll_disc_bandwidth, Δt, sat_prn)
    signals -> _tracking(signals, PLL, DLL, init_carrier_replica, init_code_replicas, beamform)
end

function _tracking(signals, PLL, DLL, carrier_replica, code_replicas, beamform)
    downconverted_signals = downconvert(signals, carrier_replica)
    correlated_signals = correlate(downconverted_signals, code_replicas)
    beamformed_signal = beamform(correlated_signals)
    next_PLL, next_carrier_replica = PLL(beamformed_signal)
    next_DLL, next_code_replicas, code_phase = DLL(beamformed_signal)
    next_signal -> _tracking(next_signal, next_PLL, next_DLL, next_carrier_replica, next_code_replicas, beamform), code_phase, prompt(correlated_signals), prompt(beamformed_signal) 
end
