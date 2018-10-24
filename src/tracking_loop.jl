"""
$(SIGNATURES)

Initialize the tracking_loop by providing initial inputs to create the replicated carrier and satellite PRN code, the PLL, the DLL, and all the therfore needed componants; return a trackin_loop function.

"""
function init_tracking(system::Union{AbstractGNSSSystem, NTuple{N, AbstractGNSSSystem}}, inits::Union{Initials, NTuple{N, Initials}}, max_total_integration_time, sample_freq, interm_freq, pll_bandwidth, dll_bandwidth, sat_prn) where N
    check_init_track_consistency(system, sample_freq)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    initial_dopplers = init_dopplers(inits)
    initial_phases = init_phases(inits)
    beamformed_prompt_correlator_buffer = create_beamformed_correlator_buffer(system)
    max_total_integration_samples = floor.(Int, max_total_integration_time .* sample_freq)
    return (signal, beamform, velocity_aiding = default_velocity_aiding(system)) -> begin
        initial_track_results = init_track_results(signal, 0, max_total_integration_samples)
        correlated_signals = init_correlated_signals(signal)
        _tracking!(initial_track_results, correlated_signals, system, signal, 1, max_total_integration_time, 0, beamform, beamformed_prompt_correlator_buffer, sample_freq, interm_freq, inits, initial_dopplers, initial_phases, carrier_loop, code_loop, sat_prn, velocity_aiding)
    end
end

function _tracking!(track_results, correlated_signals, system, signal, start_sample, max_total_integration_time, integrated_samples, beamform, beamformed_prompt_correlator_buffer, sample_freq, interm_freq, inits, dopplers, phases, carrier_loop, code_loop, sat_prn, velocity_aiding)
    code_sample_shift, actual_code_phase_shift = calc_sample_shift(system, sample_freq, 0.5)

    phases, max_total_integration_time = check_found_neuman_hofman_code(system, phases, beamformed_prompt_correlator_buffer, max_total_integration_time)

    next_correlated_signals, next_phases, next_start_sample, next_integrated_samples, max_total_integration_samples =
        gen_replica_downconvert_correlate!(correlated_signals, system, signal, max_total_integration_time, integrated_samples, start_sample, sample_freq, interm_freq, code_sample_shift, dopplers, phases, sat_prn, velocity_aiding)

    if all(next_integrated_samples .== max_total_integration_samples)
        next_correlated_signals_scaled = next_correlated_signals ./ max_total_integration_samples
        beamformed_signals, next_carrier_loop, next_code_loop, next_carrier_freq_update, next_code_freq_update =
            beamform_and_update_loops(next_correlated_signals_scaled, max_total_integration_time, beamform, carrier_loop, code_loop, actual_code_phase_shift)
        next_dopplers = aid_dopplers(system, inits, next_carrier_freq_update, next_code_freq_update, velocity_aiding)
        next_phases = find_neuman_hofman_code(system, next_phases, beamformed_prompt_correlator_buffer, prompt(beamformed_signals))
        track_results = push_track_results!(track_results, next_correlated_signals_scaled, prompt(beamformed_signals), start_sample, max_total_integration_samples, next_dopplers, next_phases)
        next_integrated_samples = 0
        next_correlated_signals = init_correlated_signals(signal)
    else
        next_dopplers = dopplers
        next_carrier_loop, next_code_loop = carrier_loop, code_loop
    end
    if check_need_new_signal(next_start_sample, signal)
        return (next_signal, next_beamform, next_velocity_aiding = default_velocity_aiding(system)) -> begin
            initial_track_results = init_track_results(next_signal, next_integrated_samples, max_total_integration_samples)
            _tracking!(initial_track_results, next_correlated_signals, system, next_signal, 1, max_total_integration_time, next_integrated_samples, next_beamform, beamformed_prompt_correlator_buffer, sample_freq, interm_freq, inits, next_dopplers, next_phases, next_carrier_loop, next_code_loop, sat_prn, next_velocity_aiding)
        end, track_results
    else
        _tracking!(track_results, next_correlated_signals, system, signal, next_start_sample, max_total_integration_time, next_integrated_samples, beamform, beamformed_prompt_correlator_buffer, sample_freq, interm_freq, inits, next_dopplers, next_phases, next_carrier_loop, next_code_loop, sat_prn, velocity_aiding)
    end
end

function create_beamformed_correlator_buffer(system::GPSL5)
    CircularBuffer{Float64}(length(system.neuman_hofman_code))
end

function create_beamformed_correlator_buffer(system)
    CircularBuffer{Float64}(0)
end

function check_found_neuman_hofman_code(system::GPSL5, phases, beamformed_prompt_correlator_buffer, max_total_integration_time)
    if isfull(beamformed_prompt_correlator_buffer)
        return phases, max_total_integration_time
    else
        return TrackingPhases(phases.carrier, mod(phases.code, system.code_length / length(system.neuman_hofman_code))), 1ms
    end
end

function check_found_neuman_hofman_code(system, phases, beamformed_prompt_correlator_buffer, max_total_integration_time)
    phases, max_total_integration_time
end

function find_neuman_hofman_code(system::GPSL5, phases, beamformed_prompt_correlator_buffer, beamformed_prompt_signal)
    if !isfull(beamformed_prompt_correlator_buffer)
        push!(beamformed_prompt_correlator_buffer, real(beamformed_prompt_signal))
        neuman_hofman_code_length = length(system.neuman_hofman_code)
        base_code_length = system.code_length / neuman_hofman_code_length
        base_code_phase = mod(phases.code, system.code_length / neuman_hofman_code_length)
        if isfull(beamformed_prompt_correlator_buffer)
            code_replica = vcat(system.neuman_hofman_code, zeros(length(beamformed_prompt_correlator_buffer) - neuman_hofman_code_length))
            cross_corr = abs.(ifft(fft(code_replica) .* conj(fft(beamformed_prompt_correlator_buffer))))
            max_value, max_index = findmax(cross_corr)
            new_code_phase = mod((max_index - 1 - (base_code_phase > base_code_length / 2)) * base_code_length + base_code_phase, system.code_length)
            return TrackingPhases(phases.carrier, new_code_phase)
        else
            return TrackingPhases(phases.carrier, base_code_phase)
        end
    else
        return phases
    end
end

function find_neuman_hofman_code(system, phases, beamformed_prompt_correlator_buffer, beamformed_prompt_signal)
    phases
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
    #beamformed_signals = map(correlated_signal -> beamform(correlated_signal .* sign(real(mean(prompt(correlated_signal))))), correlated_signals)
    beamformed_signals = map(correlated_signal -> beamform(correlated_signal), correlated_signals)
    next_carrier_loop, carrier_freq_update = carrier_loop(pll_disc(beamformed_signals), Δt)
    next_code_loop, code_freq_update = code_loop(dll_disc(beamformed_signals, 2 * actual_code_phase_shift), Δt)
    beamformed_signals, next_carrier_loop, next_code_loop, carrier_freq_update, code_freq_update
end

function push_track_results!(track_results, correlated_signals, beamformed_prompt_signal, start_sample, max_total_integration_samples, dopplers, phases)
    track_result_idx = floor(Int, start_sample / max_total_integration_samples) + 1
    #track_results[track_result_idx] = TrackingResults(dopplers.carrier, phases.carrier, dopplers.code, phases.code, prompt(correlated_signals) .* sign(real(mean(prompt(correlated_signals)))), prompt(beamformed_signals))
    track_results[track_result_idx] = TrackingResults(dopplers.carrier, phases.carrier, dopplers.code, phases.code, prompt(correlated_signals), beamformed_prompt_signal)
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
