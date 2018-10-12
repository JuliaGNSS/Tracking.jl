"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::Union{AbstractGNSSSystem, Vector{<:AbstractGNSSSystem}}, inits::Union{Initials, Vector{Initials}}, max_total_integration_time, sample_freq, interm_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    check_init_track_consistency(system, sample_freq)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    init_carrier_doppler, init_code_doppler = init_dopplers(inits)
    max_total_integration_samples = floor.(Int, max_total_integration_time .* sample_freq)
    return (signal, beamform, velocity_aiding = default_velocity_aiding(system)) -> begin
        initial_track_results = init_track_results(signal, 0, max_total_integration_samples)
        integrated_samples, corr_res = init_correlated_signals(signal, inits)
        _tracking!(initial_track_results, corr_res, system, signal, 1, max_total_integration_time, integrated_samples, beamform, sample_freq, interm_freq, init_carrier_doppler, init_code_doppler, inits, carrier_loop, code_loop, sat_prn, velocity_aiding)
    end
end

function _tracking!(track_results, corr_res, system, signal, start_sample, max_total_integration_time, integrated_samples, beamform, sample_freq, interm_freq, carrier_doppler, code_doppler, inits, carrier_loop, code_loop, sat_prn, velocity_aiding)
    code_sample_shift, actual_code_phase_shift = calc_sample_shift(system, sample_freq, 0.5)

    next_corr_res, next_start_sample, next_integrated_samples, max_total_integration_samples =
        gen_replica_downconvert_correlate!(corr_res, system, signal, max_total_integration_time, integrated_samples, start_sample, sample_freq, interm_freq, code_sample_shift, carrier_doppler, code_doppler, sat_prn, velocity_aiding)

    if all(next_integrated_samples .== max_total_integration_samples)
        next_carrier_loop, next_code_loop, next_carrier_freq_update, next_code_freq_update =
            beamform_and_update_loops(next_corr_res, max_total_integration_time, beamform, carrier_loop, code_loop, actual_code_phase_shift)
        next_carrier_doppler, next_code_doppler = aid_dopplers(system, inits, next_carrier_freq_update, velocity_aiding, next_code_freq_update)
        track_results = push_track_results!(track_results, next_corr_res, start_sample, max_total_integration_samples, next_carrier_doppler, next_code_doppler)
        next_integrated_samples, next_corr_res = init_correlated_signals(signal, next_corr_res)
    else
        next_carrier_loop, next_carrier_doppler = carrier_loop, carrier_doppler
        next_code_loop, next_code_doppler = code_loop, code_doppler
    end
    if check_need_new_signal(next_start_sample, signal)
        return (next_signal, next_beamform, next_velocity_aiding = default_velocity_aiding(system)) -> begin
            initial_track_results = init_track_results(next_signal, next_integrated_samples, max_total_integration_samples)
            _tracking!(initial_track_results, next_corr_res, system, next_signal, 1, max_total_integration_time, next_integrated_samples, next_beamform, sample_freq, interm_freq, next_carrier_doppler, next_code_doppler, inits, next_carrier_loop, next_code_loop, sat_prn, next_velocity_aiding)
        end, track_results
    else
        _tracking!(track_results, next_corr_res, system, signal, next_start_sample, max_total_integration_time, next_integrated_samples, beamform, sample_freq, interm_freq, next_carrier_doppler, next_code_doppler, inits, next_carrier_loop, next_code_loop, sat_prn, velocity_aiding)
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
    inits.carrier_doppler + carrier_doppler, inits.code_doppler + code_doppler
end

function beamform_and_update_loops(correlated_signals, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
    beamformed_signals = map(correlated_signal -> beamform(correlated_signal), correlated_signals)
    next_carrier_loop, carrier_freq_update = carrier_loop(pll_disc(beamformed_signals), Δt)
    next_code_loop, code_freq_update = code_loop(dll_disc(beamformed_signals, 2 * actual_code_phase_shift), Δt)
    next_carrier_loop, next_code_loop, carrier_freq_update, code_freq_update
end

function beamform_and_update_loops(corr_res::CorrelatorResults, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
    beamform_and_update_loops(corr_res.outputs, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
end

function push_track_results!(track_results, corr_res, start_sample, max_total_integration_samples, carrier_doppler, code_doppler)
    track_result_idx = floor(Int, start_sample / max_total_integration_samples) + 1
    track_results[track_result_idx] = TrackingResults(carrier_doppler, corr_res.carrier_phase, code_doppler, corr_res.code_phase, prompt(corr_res.outputs))
    track_results
end

function init_correlated_signals(signal, carrier_phase, code_phase)
    integrated_samples = 0
    correlated_signals = [zeros(ComplexF64, size(signal, 2)) for i = 1:3]
    corr_res = CorrelatorResults(carrier_phase, code_phase, correlated_signals)
    integrated_samples, corr_res
end

function init_correlated_signals(signal, initials::Initials)
    init_correlated_signals(signal, initials.carrier_phase, initials.code_phase)
end

function init_correlated_signals(signal, corr_res::CorrelatorResults)
    init_correlated_signals(signal, corr_res.carrier_phase, corr_res.code_phase)
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
    inits.carrier_doppler, inits.code_doppler
end

function check_init_track_consistency(system, sample_rate)
    #Everything is fine here
end
