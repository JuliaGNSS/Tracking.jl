"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::Union{AbstractGNSSSystem, Vector{<:AbstractGNSSSystem}}, inits::Union{Initials, Vector{Initials}}, max_total_integration_time, sample_freq, interm_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    check_init_track_consistency(system, sample_freq)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    initial_dopplers = init_dopplers(inits)
    initial_phases = init_phases(inits)
    max_total_integration_samples = floor.(Int, max_total_integration_time .* sample_freq)
    return (signal, beamform, velocity_aiding = default_velocity_aiding(system)) -> begin
        initial_track_results = init_track_results(signal, 0, max_total_integration_samples)
        correlated_signals = init_correlated_signals(signal)
        _tracking!(initial_track_results, correlated_signals, system, signal, 1, max_total_integration_time, 0, beamform, sample_freq, interm_freq, inits, initial_dopplers, initial_phases, carrier_loop, code_loop, sat_prn, velocity_aiding)
    end
end

function _tracking!(track_results, correlated_signals, system, signal, start_sample, max_total_integration_time, integrated_samples, beamform, sample_freq, interm_freq, inits, dopplers, phases, carrier_loop, code_loop, sat_prn, velocity_aiding)
    code_sample_shift, actual_code_phase_shift = calc_sample_shift(system, sample_freq, 0.5)

    next_correlated_signals, next_phases, next_start_sample, next_integrated_samples, max_total_integration_samples =
        gen_replica_downconvert_correlate!(correlated_signals, system, signal, max_total_integration_time, integrated_samples, start_sample, sample_freq, interm_freq, code_sample_shift, dopplers, phases, sat_prn, velocity_aiding)

    if all(next_integrated_samples .== max_total_integration_samples)
        next_carrier_loop, next_code_loop, next_carrier_freq_update, next_code_freq_update =
            beamform_and_update_loops(next_correlated_signals, max_total_integration_time, beamform, carrier_loop, code_loop, actual_code_phase_shift)
        next_dopplers = aid_dopplers(system, inits, next_carrier_freq_update, velocity_aiding, next_code_freq_update)
        track_results = push_track_results!(track_results, next_correlated_signals, start_sample, max_total_integration_samples, next_dopplers, next_phases)
        next_integrated_samples = 0
        next_correlated_signals = init_correlated_signals(signal)
    else
        next_dopplers = dopplers
        next_carrier_loop, next_code_loop = carrier_loop, code_loop
    end
    if check_need_new_signal(next_start_sample, signal)
        return (next_signal, next_beamform, next_velocity_aiding = default_velocity_aiding(system)) -> begin
            initial_track_results = init_track_results(next_signal, next_integrated_samples, max_total_integration_samples)
            _tracking!(initial_track_results, next_correlated_signals, system, next_signal, 1, max_total_integration_time, next_integrated_samples, next_beamform, sample_freq, interm_freq, inits, next_dopplers, next_phases, next_carrier_loop, next_code_loop, sat_prn, next_velocity_aiding)
        end, track_results
    else
        _tracking!(track_results, next_correlated_signals, system, signal, next_start_sample, max_total_integration_time, next_integrated_samples, beamform, sample_freq, interm_freq, inits, next_dopplers, next_phases, next_carrier_loop, next_code_loop, sat_prn, velocity_aiding)
    end
end

function init_track_results(signal, integrated_samples, total_integration_samples)
    StructArray{TrackingResults}(undef, floor(Int, (size(signal, 1) + integrated_samples) / total_integration_samples))
end

function check_need_new_signal(next_start_sample, signal)
    next_start_sample == size(signal, 1) + 1
end

function aid_dopplers(system, inits, carrier_freq_update, code_freq_update, velocity_aiding)
    carrier_doppler = carrier_freq_update + velocity_aiding
    code_doppler = code_freq_update + carrier_doppler * system.code_freq / system.center_freq
    TrackingDopplers(inits.carrier_doppler + carrier_doppler, inits.code_doppler + code_doppler)
end

function beamform_and_update_loops(correlated_signals, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
    beamformed_signals = map(correlated_signal -> beamform(correlated_signal), correlated_signals)
    next_carrier_loop, carrier_freq_update = carrier_loop(pll_disc(beamformed_signals), Δt)
    next_code_loop, code_freq_update = code_loop(dll_disc(beamformed_signals, 2 * actual_code_phase_shift), Δt)
    next_carrier_loop, next_code_loop, carrier_freq_update, code_freq_update
end

function push_track_results!(track_results, correlated_signals, start_sample, max_total_integration_samples, dopplers, phases)
    track_result_idx = floor(Int, start_sample / max_total_integration_samples) + 1
    track_results[track_result_idx] = TrackingResults(dopplers.carrier, phases.carrier, dopplers.code, phases.code, prompt(correlated_signals))
    track_results
end

function init_correlated_signals(signal)
    [zeros(ComplexF64, size(signal, 2)) for i = 1:3]
end

function init_phases(initials::Initials)
    TrackingPhases(initials.carrier_phase, initials.code_phase)
end

function init_carrier_loop(bandwidth)
    init_3rd_order_bilinear_loop_filter(bandwidth)
end

function init_code_loop(bandwidth)
    init_2nd_order_bilinear_loop_filter(bandwidth)
end

function default_velocity_aiding(system)
    0.0Hz
end

function init_dopplers(inits)
    TrackingDopplers(inits.carrier_doppler, inits.code_doppler)
end

function check_init_track_consistency(system, sample_rate)
    #Everything is fine here
end
