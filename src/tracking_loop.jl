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
    x * replica
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
- `pll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth for the PLL discriminator, in Hz
- `dll_disc_bandwidth::Float`: the signal aquivalent noise bandwidth for the DLL discriminator, in Hz
- `sat_prn::Integer`: the satellite PRN code number, choose from 1-32
- `scale_factor::Float`:  scale filter that is apllied to the carrier outputs before adding them to the code loop

# Examples
```julia-repl
    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end
    scale_factor = 1.023e6/1575.43e6
    velocity_aiding = 0
    test_signal = cis.(2 * π * 10 / 120 * (1:12))
    incoming_signals = [test_signal, test_signal, test_signal, test_signal]
    tracking_loop = Tracking.init_tracking(Tracking.init_PLL, Tracking.init_DLL, 0, 50, 0, 1023e3, 1e-3, 4e6, beamform, 12, 18.0, 1.0, 1, scale_factor)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals, velocity_aiding)
```
    """
function init_tracking(init_carrier_phase, init_carrier_freq, init_code_phase, init_code_freq, Δt, f_s, beamform, pll_disc_bandwidth, dll_disc_bandwidth, sat_prn, scale_factor)
    PLL, init_carrier_replica = init_PLL(init_carrier_phase, init_carrier_freq, f_s, pll_disc_bandwidth, Δt)
    DLL, init_code_replicas = init_DLL(init_code_phase, init_code_freq, f_s, dll_disc_bandwidth, Δt, sat_prn)
    (signals, velocity_aiding = 0) -> _tracking(signals, PLL, DLL, init_carrier_replica, init_code_replicas, beamform, scale_factor, velocity_aiding)
end

"""
$(SIGNATURES)

Should be initialized by init_tracking, uses the provided 'PLL', 'DLL' and 'beamform' function together with the provided antenna 'signals', the provided 'scale_factor, the 'velocity_aiding', and the replicated samples/codes 'carrier_replica' and 'code_replicas' to calculate the functions and samples/codes for the next timestep.
Returns the _tracking function for the next time step together with the the code_phase, the carrier_frequency_update, and the prompt of the correlated signals.

"""
function _tracking(signals, PLL, DLL, carrier_replica, code_replicas, beamform, scale_factor, velocity_aiding)
    downconverted_signals = downconvert(signals, carrier_replica')
    correlated_signals = map(replica -> correlate(downconverted_signals, replica), code_replicas)
    println(" \n correlated: ",correlated_signals)
    beamformed_signal = hcat(map(beamform, correlated_signals)...)
    println("beamformed: ",beamformed_signal,"  beamformed abs: ",abs(beamformed_signal))
    next_PLL, next_carrier_replica, carrier_phase, carrier_frequency_update = PLL(beamformed_signal, velocity_aiding)
    next_DLL, next_code_replicas, code_phase = DLL(beamformed_signal, carrier_frequency_update * scale_factor)
    (next_signal, next_velocity_aiding = 0) -> _tracking(next_signal, next_PLL, next_DLL, next_carrier_replica, next_code_replicas, beamform, scale_factor, next_velocity_aiding), code_phase, prompt(correlated_signals), carrier_frequency_update
end
