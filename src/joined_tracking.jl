#const FL1 = 1.57542e9
#const FL5 = 1.17645e9
#const CFL1 = 1023e3
#const CFL5 = 1023e4

function init_joined_tracking(systems::Vector{T}, tracking_results::Vector{TrackingResults}, sample_freqs, interm_freqs, pll_bandwidth, dll_bandwidth, sat_prn) where T <: AbstractGNSSSystem
    gen_code_replicas = map((system, track_result, sample_freq) -> init_code_replica(system, system.code_freq + track_result.code_doppler, track_result.code_phase, sample_freq, sat_prn), systems, tracking_results, sample_freqs)
    gen_carrier_replicas = map((track_result, interm_freq, sample_freq) -> init_carrier_replica(interm_freq + track_result.carrier_doppler, track_result.carrier_phase, sample_freq), tracking_results, interm_freqs, sample_freqs)
    carrier_loop = init_carrier_loop(pll_bandwidth)
    code_loop = init_code_loop(dll_bandwidth)
    (signals, beamform, velocity_aidings = zeros(length(systems))) -> _joined_tracking(systems, signals, sample_freqs, beamform, gen_carrier_replicas, gen_code_replicas, 0.0, 0.0, carrier_loop, code_loop, velocity_aidings)
end

function gen_replica_downconvert_correlate(system, signal, sample_freq, gen_carrier_replica, gen_code_replica, carrier_freq_update, code_freq_update, center_freq_geometric_mean, code_freq_geometric_mean, velocity_aiding)
    num_samples = size(signal, 1)
    carrier_doppler = carrier_freq_update * system.center_freq / center_freq_geometric_mean + velocity_aiding
    code_doppler = code_freq_update * system.code_freq / code_freq_geometric_mean + carrier_doppler * system.code_freq / system.center_freq

    next_gen_carrier_replica, carrier_replica, next_carrier_phase = gen_carrier_replica(num_samples, carrier_doppler)
    next_gen_code_replica, code_replicas, next_code_phase = gen_code_replica(num_samples, code_doppler)

    downconverted_signal = downconvert(signal, carrier_replica)
    correlated_signals = map(replica -> correlate(downconverted_signal, replica).', code_replicas)

    tracking_result = TrackingResults(carrier_doppler, next_carrier_phase, code_doppler, next_code_phase, prompt(correlated_signals))

    tracking_result, reduce(hcat, correlated_signals) ./ num_samples, next_gen_carrier_replica, next_gen_code_replica
end

function _joined_tracking(systems, signals, sample_freqs, beamform, gen_carrier_replicas, gen_code_replicas, carrier_freq_update, code_freq_update, carrier_loop, code_loop, velocity_aidings)
    num_systems = length(systems)
    num_samples = map(signal -> size(signal, 1), signals)
    num_ants = sum(map(signal -> size(signal, 2), signals))
    Δts = map((num_samples, sample_freq) -> num_samples / sample_freq, num_samples, sample_freqs)
    if !all(Δts .== Δts[1])
        error("Signal time interval must be identical for all given signals.")
    end
    Δt = Δts[1]

    center_freq_geometric_mean = prod([system.center_freq for system in systems])^(1 / num_systems)
    code_freq_geometric_mean = prod([system.code_freq for system in systems])^(1 / num_systems)
    
    intermediate_results = map(systems, signals, sample_freqs, gen_carrier_replicas, gen_code_replicas, velocity_aidings) do system, signal, sample_freq, gen_carrier_replica, gen_code_replica, velocity_aiding
        gen_replica_downconvert_correlate(system, signal, sample_freq, gen_carrier_replica, gen_code_replica, carrier_freq_update, code_freq_update, center_freq_geometric_mean, code_freq_geometric_mean, velocity_aiding)
    end
    tracking_results = [intermediate_result[1] for intermediate_result in intermediate_results]
    correlated_signals = [intermediate_result[2] for intermediate_result in intermediate_results]
    next_gen_carrier_replicas = [intermediate_result[3] for intermediate_result in intermediate_results]
    next_gen_code_replicas = [intermediate_result[4] for intermediate_result in intermediate_results]

    joined_correlated_signals = reduce(vcat, correlated_signals)
    beamformed_signal = beamform(joined_correlated_signals)

    next_carrier_loop, next_carrier_freq_update = carrier_loop(beamformed_signal, Δt)
    next_code_loop, next_code_freq_update = code_loop(beamformed_signal, Δt)

    (next_signals, next_beamform, next_velocity_aidings = zeros(length(systems))) -> _joined_tracking(systems, next_signals, sample_freqs, next_beamform, next_gen_carrier_replicas, next_gen_code_replicas, next_carrier_freq_update, next_code_freq_update, next_carrier_loop, next_code_loop, next_velocity_aidings), tracking_results
end

#= function init_joined_tracking(l1_system::GPSL1, l5_system::GPSL5, l1_results::TrackingResults, l5_results::TrackingResults, l1_sampling_freq, l5_sampling_freq, l1_interm_freq, l5_interm_freq, pll_bandwidth, dll_bandwidth, sat_prn)
    l1_gen_code_replica = init_code_replica(l1_system, l1_system.code_freq + l1_results.code_doppler, l1_results.code_phase, l1_sampling_freq, sat_prn)
    l5_gen_code_replica = init_code_replica(l5_system, l5_system.code_freq + l5_results.code_doppler, l5_results.code_phase, l5_sampling_freq, sat_prn)
    l1_gen_carrier_replica = init_carrier_replica(l1_interm_freq + l1_results.carrier_doppler, l1_results.carrier_phase, l1_sampling_freq)
    l5_gen_carrier_replica = init_carrier_replica(l5_interm_freq + l5_results.carrier_doppler, l5_results.carrier_phase, l5_sampling_freq)
    l1_carrier_loop = init_3rd_order_bilinear_loop_filter(pll_bandwidth)
    l5_carrier_loop = init_3rd_order_bilinear_loop_filter(pll_bandwidth)
    l1_code_loop = init_2nd_order_bilinear_loop_filter(dll_bandwidth)
    l5_code_loop = init_2nd_order_bilinear_loop_filter(dll_bandwidth)
    beamformer = [0.5, 0.5, 0.5, 0.5]
    correlator_buffer = zeros(Complex{Float64}, 4, 20)
    buffer_index = 0
    (l1_signals, l5_signals, velocity_aiding = 0.0) -> _joined_tracking(l1_signals, l5_signals, beamformer, correlator_buffer, buffer_index, l1_gen_carrier_replica, l5_gen_carrier_replica, l1_gen_code_replica, l5_gen_code_replica, 0.0, 0.0, 0.0, 0.0, l1_carrier_loop, l5_carrier_loop, l1_code_loop, l5_code_loop, l1_system.code_freq / l1_system.center_freq, l5_system.code_freq / l5_system.center_freq, velocity_aiding, l1_sampling_freq)
end

function _joined_tracking(l1_signals, l5_signals, beamformer, correlator_buffer, buffer_index, l1_gen_carrier_replica, l5_gen_carrier_replica, l1_gen_code_replica, l5_gen_code_replica, l1_carrier_freq_update, l5_carrier_freq_update, l1_code_freq_update, l5_code_freq_update, l1_carrier_loop, l5_carrier_loop, l1_code_loop, l5_code_loop, l1_aiding_scale_factor, l5_aiding_scale_factor, velocity_aiding, l1_sampling_freq)
    l1_num_samples = size(l1_signals, 1)
    l5_num_samples = size(l5_signals, 1)
    Δt = l1_num_samples / l1_sampling_freq
    l1_carrier_doppler = l1_carrier_freq_update + velocity_aiding
    l5_carrier_doppler = l5_carrier_freq_update + velocity_aiding
    l1_code_doppler = l1_code_freq_update + l1_carrier_doppler * l1_aiding_scale_factor
    l5_code_doppler = l5_code_freq_update + l5_carrier_doppler * l5_aiding_scale_factor

    next_l1_gen_carrier_replica, l1_carrier_replica, next_l1_carrier_phase = l1_gen_carrier_replica(l1_num_samples, l1_carrier_doppler)
    next_l5_gen_carrier_replica, l5_carrier_replica, next_l5_carrier_phase = l5_gen_carrier_replica(l5_num_samples, l5_carrier_doppler) #scale factor for carrier_doppler

    next_l1_gen_code_replica, l1_code_replicas, next_l1_code_phase = l1_gen_code_replica(l1_num_samples, l1_code_doppler)
    next_l5_gen_code_replica, l5_code_replicas, next_l5_code_phase = l5_gen_code_replica(l5_num_samples, l5_code_doppler)

    l1_downconverted_signals = downconvert(l1_signals, l1_carrier_replica)
    l5_downconverted_signals = downconvert(l5_signals, l5_carrier_replica)

    l1_correlated_signals = map(replica -> correlate(l1_downconverted_signals, replica).', l1_code_replicas)
    l5_correlated_signals = map(replica -> correlate(l5_downconverted_signals, replica).', l5_code_replicas)

    correlated_signals = vcat(hcat(l1_correlated_signals...), hcat(l5_correlated_signals...) ./ 10)
    beamformed_signal = beamformer' * correlated_signals

    correlator_buffer[:, mod(buffer_index, 20) + 1] = correlated_signals[:,2]
    if mod(buffer_index, 20) + 1 == 20
        u, U = eig(correlator_buffer * correlator_buffer')
        b = sortperm(abs.(u))
        U_sorted = U[:,b]
        next_beamformer = U_sorted[:,end]
        beamformer = next_beamformer' * beamformer * next_beamformer / abs(next_beamformer' * beamformer)
    end

    carrier_phase_error = pll_disc(beamformed_signal)
    code_phase_error = dll_disc(beamformed_signal)

    l1_carrier_phase_error = carrier_phase_error * 2 * FL1 / (FL1 + FL5)
    l5_carrier_phase_error = carrier_phase_error * 2 * FL5 / (FL1 + FL5)
    l1_code_phase_error = code_phase_error * 2 * CFL1 / (CFL1 + CFL5)
    l5_code_phase_error = code_phase_error * 2 * CFL5 / (CFL1 + CFL5)

    next_l1_carrier_loop, next_l1_carrier_freq_update = l1_carrier_loop(l1_carrier_phase_error, Δt)
    next_l5_carrier_loop, next_l5_carrier_freq_update = l5_carrier_loop(l5_carrier_phase_error, Δt)
    next_l1_code_loop, next_l1_code_freq_update = l1_code_loop(l1_code_phase_error, Δt)
    next_l5_code_loop, next_l5_code_freq_update = l5_code_loop(l5_code_phase_error, Δt)

    l1_results = TrackingResults(l1_carrier_doppler, next_l1_carrier_phase, l1_code_doppler, next_l1_code_phase, prompt(l1_correlated_signals))
    l5_results = TrackingResults(l5_carrier_doppler, next_l5_carrier_phase, l5_code_doppler, next_l5_code_phase, prompt(l5_correlated_signals))
    joined_results = JoinedTrackingResults(l1_results, l5_results)

    (next_l1_signals, next_l5_signals, next_velocity_aiding = 0.0) -> _joined_tracking(next_l1_signals, next_l5_signals, beamformer, correlator_buffer, buffer_index + 1, next_l1_gen_carrier_replica, next_l5_gen_carrier_replica, next_l1_gen_code_replica, next_l5_gen_code_replica, next_l1_carrier_freq_update, next_l5_carrier_freq_update, next_l1_code_freq_update, next_l5_code_freq_update, next_l1_carrier_loop, next_l5_carrier_loop, next_l1_code_loop, next_l5_code_loop, l1_aiding_scale_factor, l5_aiding_scale_factor, next_velocity_aiding, l1_sampling_freq), joined_results
end
 =#
