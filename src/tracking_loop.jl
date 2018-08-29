"""
$(SIGNATURES)

Take one or multiple antenna signals 'x' and downgrade them by componantwise multiplying with a complex conjugated carrier replication signal 'replica' and return the result.

"""
function downconvert(x, replica)
    x .* conj(replica)
end

"""
$(SIGNATURES)

Take one or multiple, allready downconverted antenna signals 'x' and an replicated satellite code , return the scalarproduct of each combination.

"""
function correlate(x, replica)
    replica.' * x
end

"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::AbstractGNSSSystem, inits::Initials, interm_freq, sampling_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    aiding_scale_factor = system.code_freq / system.center_freq
    gen_code_replica = init_code_replica(system, system.code_freq + inits.code_doppler, inits.code_phase, sampling_freq, sat_prn)
    gen_carrier_replica = init_carrier_replica(interm_freq + inits.carrier_doppler, inits.carrier_phase, sampling_freq)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    (signals, beamform, velocity_aiding = 0.0) -> _tracking(signals, beamform, sampling_freq, gen_carrier_replica, gen_code_replica, 0.0, 0.0, carrier_loop, code_loop, aiding_scale_factor, velocity_aiding)
end

"""
$(SIGNATURES)

Should be initialized by init_tracking, uses the provided 'PLL', 'DLL' and 'beamform' function together with the provided antenna 'signals', the provided 'aiding_scale_factor, the 'velocity_aiding', and the replicated samples/codes 'carrier_replica' and 'code_replicas' to calculate the functions and samples/codes for the next timestep.
Returns the _tracking function for the next time step together with the the code_phase, the carrier_frequency_update, and the prompt of the correlated signals.

"""
function _tracking(signals, beamform, sampling_freq, gen_carrier_replica, gen_code_replica, carrier_doppler, code_doppler, carrier_loop, code_loop, aiding_scale_factor, velocity_aiding)
    num_samples = size(signals, 1)
    Δt =  num_samples / sampling_freq
    next_gen_carrier_replica, carrier_replica, next_carrier_phase = gen_carrier_replica(num_samples, carrier_doppler)
    next_gen_code_replica, code_replicas, next_code_phase = gen_code_replica(num_samples, code_doppler)
    downconverted_signals = downconvert(signals, carrier_replica)
    correlated_signals = map(replica -> correlate(downconverted_signals, replica).', code_replicas)
    beamformed_signal = hcat(map(beamform, correlated_signals)...)
    next_carrier_loop, carrier_freq_update = carrier_loop(beamformed_signal, Δt)
    next_code_loop, code_freq_update = code_loop(beamformed_signal, Δt)
    next_carrier_doppler = carrier_freq_update + velocity_aiding
    next_code_doppler = code_freq_update + next_carrier_doppler * aiding_scale_factor
    (next_signal, beamform, next_velocity_aiding = 0.0) -> _tracking(next_signal, beamform, sampling_freq, next_gen_carrier_replica, next_gen_code_replica, next_carrier_doppler, next_code_doppler, next_carrier_loop, next_code_loop, aiding_scale_factor, velocity_aiding), TrackingResults(next_carrier_doppler, next_carrier_phase, next_code_doppler, next_code_phase, prompt(correlated_signals))
end

function init_carrier_loop(bandwidth)
    (correlator_output, Δt) -> _loop(correlator_output, pll_disc, init_3rd_order_bilinear_loop_filter(bandwidth), Δt)
end

function init_code_loop(bandwidth)
    (correlator_output, Δt) -> _loop(correlator_output, dll_disc, init_2nd_order_bilinear_loop_filter(bandwidth), Δt)
end

function _loop(correlator_output, disc, loop_filter, Δt)
    next_loop_filter, freq_update = loop_filter(disc(correlator_output), Δt)
    (next_correlator_output, next_Δt) -> _loop(next_correlator_output, disc, next_loop_filter, next_Δt), freq_update
end