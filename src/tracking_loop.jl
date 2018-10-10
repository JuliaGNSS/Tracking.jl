"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::AbstractGNSSSystem, inits::Initials, max_total_integration_time, sample_freq, interm_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    max_total_integration_samples = floor(Int, max_total_integration_time * sample_freq)
    return (signal, beamform, velocity_aiding = 0.0Hz) -> begin
        num_ants = size(signal, 2)
        correlated_signals = [zeros(ComplexF64, num_ants) for i = 1:3] # Early, prompt, late
        init_track_results = StructArray{TrackingResults}(undef, floor(Int, size(signal, 1) / max_total_integration_samples))
        _tracking!(init_track_results, correlated_signals, system, signal, 1, max_total_integration_time, 0, beamform, sample_freq, interm_freq, inits.carrier_phase, inits.code_phase, 0.0Hz, 0.0Hz, inits.carrier_doppler, inits.code_doppler, carrier_loop, code_loop, sat_prn, velocity_aiding)
    end
end

function _tracking!(track_results, correlated_signals, system, signal, start_sample, max_total_integration_time, integrated_samples, beamform, sample_freq, interm_freq, carrier_phase, code_phase, carrier_freq_update, code_freq_update, init_carrier_doppler, init_code_doppler, carrier_loop, code_loop, sat_prn, velocity_aiding)
    code_sample_shift, actual_code_phase_shift = calc_sample_shift(system, sample_freq, 0.5)

    correlated_signals, carrier_doppler, code_doppler, next_carrier_phase, next_code_phase, next_start_sample, next_integrated_samples, max_total_integration_samples =
        gen_replica_downconvert_correlate!(correlated_signals, system, signal, max_total_integration_time, integrated_samples, start_sample, sample_freq, interm_freq, code_sample_shift, init_carrier_doppler, init_code_doppler, carrier_freq_update, code_freq_update, carrier_phase, code_phase, sat_prn, velocity_aiding)

    if next_integrated_samples == max_total_integration_samples
        next_carrier_loop, next_code_loop, next_carrier_freq_update, next_code_freq_update =
            beamform_and_update_loops(correlated_signals, max_total_integration_time, beamform, carrier_loop, code_loop)
        track_results = push_track_results!(track_results, start_sample, max_total_integration_samples, carrier_doppler, carrier_phase, code_doppler, code_phase, correlated_signals)
        next_integrated_samples, correlated_signals = init_correlated_signals(signal)
    else
        next_carrier_loop, next_carrier_freq_update = carrier_loop, carrier_freq_update
        next_code_loop, next_code_freq_update = code_loop, code_freq_update
    end
    if next_start_sample < size(signal, 1)
        _tracking!(track_results, correlated_signals, system, signal, next_start_sample, max_total_integration_time, next_integrated_samples, beamform, sample_freq, interm_freq, next_carrier_phase, next_code_phase, next_carrier_freq_update, next_code_freq_update, init_carrier_doppler, init_code_doppler, next_carrier_loop, next_code_loop, sat_prn, velocity_aiding)
    else
        return (next_signal, next_beamform, next_velocity_aiding = 0.0Hz) -> begin
            init_track_results = StructArray{TrackingResults}(undef, floor(Int, (size(next_signal, 1) + next_integrated_samples) / max_total_integration_samples))
            _tracking!(init_track_results, correlated_signals, system, next_signal, 1, max_total_integration_time, next_integrated_samples, next_beamform, sample_freq, interm_freq, next_carrier_phase, next_code_phase, next_carrier_freq_update, next_code_freq_update, init_carrier_doppler, init_code_doppler, next_carrier_loop, next_code_loop, sat_prn, next_velocity_aiding)
        end, track_results
    end
end

function check_need_new_signal(next_start_sample, signal)
    next_start_sample == size(signal, 1)
end

function beamform_and_update_loops(correlated_signals, Δt, beamform, carrier_loop, code_loop)
    beamformed_signals = map(correlated_signal -> beamform(correlated_signal), correlated_signals)
    next_carrier_loop, carrier_freq_update = carrier_loop(beamformed_signals, Δt)
    next_code_loop, code_freq_update = code_loop(beamformed_signals, Δt)
    next_carrier_loop, next_code_loop, carrier_freq_update, code_freq_update
end

function push_track_results!(track_results, start_sample, max_total_integration_samples, carrier_doppler, carrier_phase, code_doppler, code_phase, correlated_signals)
    track_result_idx = floor(Int, start_sample / max_total_integration_samples) + 1
    track_results[track_result_idx] = TrackingResults(carrier_doppler, carrier_phase, code_doppler, code_phase, prompt(correlated_signals))
    track_results
end

function init_correlated_signals(signal)
    next_integrated_samples = 0
    correlated_signals = [zeros(ComplexF64, size(signal, 2)) for i = 1:3]
    next_integrated_samples, correlated_signals
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
