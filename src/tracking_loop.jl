struct TrackingResults
    carrier_doppler::Float64
    carrier_phase::Float64
    code_doppler::Float64
    code_phase::Float64
    prompts_correlated_signals::Array{Complex{Float64}, 1}
end

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
- `scale_factor::Float`: scale filter that is apllied to the carrier outputs before adding them to the code loop

# Examples
```julia-repl
    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end
    scale_factor = 1.023e6/1575.43e6
    velocity_aiding = 0.0
    test_signal = cis.(2 * π * 10 / 120 * (1:12))
    incoming_signals = [test_signal, test_signal, test_signal, test_signal]
    tracking_loop = Tracking.init_tracking(Tracking.init_PLL, Tracking.init_DLL, 0, 50, 0, 1023e3, 1e-3, 4e6, beamform, 12, 18.0, 1.0, 1, scale_factor)
    next_tracking_loop, code_phase, prompt_correlated_signal, prompt_beamformed_signal = tracking_loop(incoming_signals, velocity_aiding)
```
    """
function init_tracking(init_carrier_phase, carrier_freq, interm_freq, init_carrier_doppler, init_code_phase, code_freq, sampling_freq, pll_disc_bandwidth, dll_disc_bandwidth, sat_svid, gen_sampled_code, get_code_phase)
    gen_code_replica = init_code_replica(init_code_phase, sampling_freq, sat_svid, gen_sampled_code, get_code_phase)
    gen_carrier_replica = init_carrier_replica(init_carrier_phase, sampling_freq)
    carrier_loop = init_carrier_loop(pll_disc_bandwidth)
    code_loop = init_code_loop(dll_disc_bandwidth)
    aiding_scale_factor = code_freq / carrier_freq
    code_doppler = init_carrier_doppler * aiding_scale_factor
    (signals, beamform, velocity_aiding = 0.0) -> _tracking(signals, beamform, sampling_freq, gen_carrier_replica, gen_code_replica, interm_freq + init_carrier_doppler, 0.0, code_freq + code_doppler, 0.0, carrier_loop, code_loop, aiding_scale_factor, velocity_aiding)
end

"""
$(SIGNATURES)

Should be initialized by init_tracking, uses the provided 'PLL', 'DLL' and 'beamform' function together with the provided antenna 'signals', the provided 'aiding_scale_factor, the 'velocity_aiding', and the replicated samples/codes 'carrier_replica' and 'code_replicas' to calculate the functions and samples/codes for the next timestep.
Returns the _tracking function for the next time step together with the the code_phase, the carrier_frequency_update, and the prompt of the correlated signals.

"""
function _tracking(signals, beamform, sampling_freq, gen_carrier_replica, gen_code_replica, init_carrier_freq, carrier_doppler, init_code_freq, code_doppler, carrier_loop, code_loop, aiding_scale_factor, velocity_aiding)
    num_samples = size(signals, 2)
    Δt =  num_samples / sampling_freq
    next_gen_carrier_replica, carrier_replica, next_carrier_phase = gen_carrier_replica(num_samples, init_carrier_freq + carrier_doppler)
    next_gen_code_replica, code_replicas, next_code_phase = gen_code_replica(num_samples, init_code_freq + code_doppler)
    downconverted_signals = downconvert(signals, carrier_replica')
    correlated_signals = map(replica -> correlate(downconverted_signals, replica), code_replicas)
    beamformed_signal = hcat(map(beamform, correlated_signals)...)
    next_carrier_loop, carrier_freq_update = carrier_loop(beamformed_signal, Δt)
    next_code_loop, code_freq_update = code_loop(beamformed_signal, Δt)
    next_carrier_doppler = carrier_freq_update + velocity_aiding
    next_code_doppler = code_freq_update + next_carrier_doppler * aiding_scale_factor
    (next_signal, beamform, next_velocity_aiding = 0.0) -> _tracking(next_signal, beamform, sampling_freq, next_gen_carrier_replica, next_gen_code_replica, init_carrier_freq, next_carrier_doppler, init_code_freq, next_code_doppler, next_carrier_loop, next_code_loop, aiding_scale_factor, velocity_aiding), TrackingResults(next_carrier_doppler, next_carrier_phase, next_code_doppler, next_code_phase, prompt(correlated_signals))
end

function init_carrier_loop(bandwidth)
    (correlator_output, Δt) -> _loop(correlator_output, pll_disc, init_3rd_order_loop_filter(bandwidth), Δt)
end

function init_code_loop(bandwidth)
    (correlator_output, Δt) -> _loop(correlator_output, dll_disc, init_2nd_order_loop_filter(bandwidth), Δt)
end

function _loop(correlator_output, disc, loop_filter, Δt)
    next_loop_filter, freq_update = loop_filter(disc(correlator_output), Δt)
    (next_correlator_output, next_Δt) -> _loop(next_correlator_output, disc, next_loop_filter, next_Δt), freq_update
end