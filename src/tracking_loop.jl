"""
$(SIGNATURES)

Take one or multiple antenna signals `x` and downgrade them by componantwise multiplying with a complex conjugated carrier replication signal `replica` and return the result.

"""
function downconvert(x, replica)
    x .* conj(replica)
end

"""
$(SIGNATURES)

Take one or multiple, allready downconverted antenna signals `x` and an replicated satellite code , return the scalarproduct of each combination.

"""
function correlate(x, replica)
    transpose(replica) * x / size(x, 1)
end

"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::AbstractGNSSSystem, inits::Initials, interm_freq, sample_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    gen_code_replica = init_code_replica(system, system.code_freq + inits.code_doppler, inits.code_phase, sample_freq, sat_prn)
    gen_carrier_replica = init_carrier_replica(interm_freq + inits.carrier_doppler, inits.carrier_phase, sample_freq)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    (signal, beamform, velocity_aiding = 0.0Hz) -> _tracking(system, signal, beamform, sample_freq, gen_carrier_replica, gen_code_replica, 0.0Hz, 0.0Hz, carrier_loop, code_loop, velocity_aiding)
end

"""
$(SIGNATURES)

Should be initialized by init_tracking, uses the provided `PLL`, `DLL` and `beamform` function together with the provided antenna `system`, provided antenna `signals`, the `velocity_aiding`, and the replicated samples/codes `carrier_replica` and `code_replicas` to calculate the functions and samples/codes for the next timestep.
Returns the _tracking function for the next time step together with the the code_phase, the carrier_frequency_update, and the prompt of the correlated signals.

"""
function _tracking(system, signal, beamform, sample_freq, gen_carrier_replica, gen_code_replica, carrier_freq_update, code_freq_update, carrier_loop, code_loop, velocity_aiding)
    num_samples = size(signal, 1)
    Δt =  num_samples / sample_freq

    tracking_result, correlated_signals, next_gen_carrier_replica, next_gen_code_replica =
        gen_replica_downconvert_correlate(system, signal, sample_freq, gen_carrier_replica, gen_code_replica, carrier_freq_update, code_freq_update, system.center_freq, system.code_freq, velocity_aiding)
    beamformed_signal = beamform(correlated_signals)
    next_carrier_loop, next_carrier_freq_update = carrier_loop(beamformed_signal, Δt)
    next_code_loop, next_code_freq_update = code_loop(beamformed_signal, Δt)
    (next_signal, next_beamform, next_velocity_aiding = 0.0Hz) -> _tracking(system, next_signal, next_beamform, sample_freq, next_gen_carrier_replica, next_gen_code_replica, next_carrier_freq_update, next_code_freq_update, next_carrier_loop, next_code_loop, velocity_aiding), tracking_result
end

function gen_replica_downconvert_correlate(system, signal, sample_freq, gen_carrier_replica, gen_code_replica, carrier_freq_update, code_freq_update, center_freq_geometric_mean, code_freq_geometric_mean, velocity_aiding)
    num_samples = size(signal, 1)
    carrier_doppler = carrier_freq_update * system.center_freq / center_freq_geometric_mean + velocity_aiding
    code_doppler = code_freq_update * system.code_freq / code_freq_geometric_mean + carrier_doppler * system.code_freq / system.center_freq

    next_gen_carrier_replica, carrier_replica, next_carrier_phase = gen_carrier_replica(num_samples, carrier_doppler)
    next_gen_code_replica, code_replicas, next_code_phase = gen_code_replica(num_samples, code_doppler)

    downconverted_signal = downconvert(signal, carrier_replica)
    correlated_signals = map(replica -> transpose(correlate(downconverted_signal, replica)), code_replicas)

    tracking_result = TrackingResults(carrier_doppler, next_carrier_phase, code_doppler, next_code_phase, prompt(correlated_signals))

    tracking_result, reduce(hcat, correlated_signals), next_gen_carrier_replica, next_gen_code_replica
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
