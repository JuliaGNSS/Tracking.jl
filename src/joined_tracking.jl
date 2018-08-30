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