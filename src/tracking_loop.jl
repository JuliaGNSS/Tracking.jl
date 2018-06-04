"""
$(SIGNATURES)

Take multiple antenna signals 'x' and downgrade them by componantwise multiplying of each individual signal with a complex conjugated carrier replication signal 'replica' and return the result.

"""
function _downconvert(x, replica)
    map(signal -> signal .* conj(replica), x)
end

"""
$(SIGNATURES)

Take multiple, allready downconverted antenna signals 'x' and multiple replicated satellite codes (e.g. with different chip offsets), return the scalarproduct of each combination.

"""
function _correlate(x, replicas)
    map(replica -> map(single_antenna_signal -> replica' * single_antenna_signal, x), replicas)
end

"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

# Arguments
- `init_carrier_phase::Float`: the initial replicated carrier signal phase in rad
- `init_carrier_freq::Float`: the initial replicated carrier signal frequency in Hz
- `init_carrier_phase::Float`: the initial replicated PRN code offset in chip
- `init_carrier_freq::Float`: the initial replicated PRN code signal frequency in Hz
- `Δt::Float`: the loop update time intervall in seconds
- `f_s::Float`: the signal sampling frequency in Hz
- `beamform::Function`: a beamforming function to produce a new early_prompt_late tuple from multiple inputs
- `num_samples::Integer`: the amount of samples per carrier or satellite code 
- `pll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth  for the PLL discriminator, in Hz
- `dll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth  for the DLL discriminator, in Hz
- `sat_prn::Integer`:the satellite PRN code number, choose from 1-32

# Examples
```julia-repl
    function beamform(x)
        beamformed_x = map(x_per_phaseshift -> [0.5 0.5 0.5 0.5] * (x_per_phaseshift), x)
        hcat(beamformed_x...)
    end
    test_signal = cis.(2 * π * 10 / 120 * (1:12))
    incoming_signals = [test_signal, test_signal, test_signal, test_signal]
    tracking_loop = Tracking.init_tracking(Tracking.init_PLL, Tracking.init_DLL, 0, 50, 0, 1023e3, 1e-3, 4e6, beamform, 12, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals)
```
    """
function init_tracking(init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, num_samples, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn)
    PLL, init_carrier_replica = init_PLL(init_carrier_phase, init_carrier_freq, num_samples,  f_s, pll_disc_bandwidth, Δt)
    DLL, init_code_replicas = init_DLL(init_code_phase, init_code_freq, num_samples,  f_s, dll_disc_bandwidth, Δt, sat_prn)
    signals -> _tracking(signals, PLL, DLL, init_carrier_replica, init_code_replicas, beamform)
end

"""
$(SIGNATURES)

Should be initialized by init_tracking, uses the provided 'PLL', 'DLL' and 'beamform' function together with the provided antenna 'signals' and the replicated samples/codes 'carrier_replica' and 'code_replicas' to calculate the functions and samples/codes for the next timestep.
Returns the _tracking function for the next time step together with the the code_phase and the prompt of the correlated and beamformed signals.


# Examples
```julia-repl
    function beamform(x)
        beamformed_x = map(x_per_phaseshift -> [0.5 0.5 0.5 0.5] * (x_per_phaseshift), x)
        hcat(beamformed_x...)
    end
    test_signal = cis.(2 * π * 10 / 120 * (1:12))
    incoming_signals = [test_signal, test_signal, test_signal, test_signal]
    tracking_loop = Tracking.init_tracking(Tracking.init_PLL, Tracking.init_DLL, 0, 50, 0, 1023e3, 1e-3, 4e6, beamform, 12, 18.0, 1.0, 1)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals)
```
"""
function _tracking(signals, PLL, DLL, carrier_replica, code_replicas, beamform)
    downconverted_signals = _downconvert(signals, carrier_replica)
    correlated_signals = _correlate(downconverted_signals, code_replicas)
    beamformed_signal = beamform(correlated_signals)
    next_PLL, next_carrier_replica = PLL(beamformed_signal)
    next_DLL, next_code_replicas, code_phase = DLL(beamformed_signal)
    next_signal -> _tracking(next_signal, next_PLL, next_DLL, next_carrier_replica, next_code_replicas, beamform), code_phase, prompt(correlated_signals), prompt(beamformed_signal) 
end
