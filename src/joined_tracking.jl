
struct GNSSSystem
    f₀::Float64 # Carrier center frequency
    code_f₀::Float64 # Code center frequency
    sampling_freq::Float64
    calc_next_code_phase::Function
    gen_sampled_code::Function
end

struct JoinedTrackingResults
    l1_results::TrackingResults
    l5_results::TrackingResults
end

function init_joined_tracking(l1_system::GNSSSystem, l5_system::GNSSSystem, l1_results::TrackingResults, l5_results::TrackingResults, l1_init_carrier_freq, l5_init_carrier_freq, pll_disc_bandwidth, dll_disc_bandwidth, sat_svid)
    #sampling_freq durch Δt ersetzen und Δt vorher ausrechnen
    l1_gen_code_replica = init_code_replica(l1_results.code_phase, l1_system.sampling_freq, sat_svid, l1_system.gen_sampled_code, l1_system.calc_next_code_phase)
    l5_gen_code_replica = init_code_replica(l5_results.code_phase, l5_system.sampling_freq, sat_svid, l5_system.gen_sampled_code, l5_system.calc_next_code_phase)
    l1_gen_carrier_replica = init_carrier_replica(l1_results.carrier_phase, l1_system.sampling_freq)
    l5_gen_carrier_replica = init_carrier_replica(l5_results.carrier_phase, l5_system.sampling_freq)
    carrier_loop = init_carrier_loop(pll_disc_bandwidth) # pll
    code_loop = init_code_loop(dll_disc_bandwidth) # dll
    (l1_signals, l5_signals, beamform, velocity_aiding = 0.0) -> _joined_tracking(l1_signals, l5_signals, beamform, l1_gen_carrier_replica, l5_gen_carrier_replica, l1_gen_code_replica, l5_gen_code_replica, l1_init_carrier_freq + l1_results.carrier_doppler, l5_init_carrier_freq + l5_results.carrier_doppler, l1_system.code_f₀ + l1_results.code_doppler, l5_system.code_f₀ + l5_results.code_doppler, #=carrier_freq_update=# 0.0, #=code_freq_update=# 0.0, carrier_loop, code_loop, #=l1_aiding_scale_factor=# l1_system.code_f₀ / l1_system.f₀, #=l5_aiding_scale_factor=# l5_system.code_f₀ / l5_system.f₀, velocity_aiding)
end

function _joined_tracking(l1_signals, l5_signals, beamform, l1_gen_carrier_replica, l5_gen_carrier_replica, l1_gen_code_replica, l5_gen_code_replica, l1_init_carrier_freq, l5_init_carrier_freq, l1_init_code_freq, l5_init_code_freq, carrier_freq_update, code_freq_update, carrier_loop, code_loop, l1_aiding_scale_factor, l5_aiding_scale_factor, velocity_aiding)
    l1_num_samples = size(l1_signals, 2)
    l5_num_samples = size(l5_signals, 2)
    Δt =  l1_num_samples / sampling_freq
    l1_carrier_doppler = (carrier_freq_update + velocity_aiding) / sqrt(115/154)
    l5_carrier_doppler = (carrier_freq_update + velocity_aiding) * sqrt(115/154)
    l1_code_doppler = code_freq_update / sqrt(10) + l1_carrier_doppler *  l1_aiding_scale_factor
    l5_code_doppler = code_freq_update * sqrt(10) + l5_carrier_doppler *  l5_aiding_scale_factor

    next_l1_gen_carrier_replica, l1_carrier_replica, l1_carrier_phase = l1_gen_carrier_replica(l1_num_samples, l1_init_carrier_freq + l1_carrier_doppler)
    next_l5_gen_carrier_replica, l5_carrier_replica, l5_carrier_phase = l5_gen_carrier_replica(l5_num_samples, l5_init_carrier_freq + l5_carrier_doppler) #scale factor for carrier_doppler

    next_l1_gen_code_replica, l1_code_replicas, l1_code_phase = l1_gen_code_replica(l1_num_samples, l1_init_code_freq + l1_code_doppler)
    next_l5_gen_code_replica, l5_code_replicas, l5_code_phase = l5_gen_code_replica(l5_num_samples, l5_init_code_freq + l5_code_doppler)

    l1_downconverted_signals = downconvert(l1_signals, l1_carrier_replica')
    l5_downconverted_signals = downconvert(l5_signals, l5_carrier_replica')

    l1_correlated_signals = map(replica -> correlate(l1_downconverted_signals, replica), l1_code_replicas)
    l5_correlated_signals = map(replica -> correlate(l5_downconverted_signals, replica), l5_code_replicas)

    #ToDo beamform with PLL to filter data, for not simulated signals
    correlated_signals = vcat(l1_correlated_signals, l5_correlated_signals)
    beamformed_signal = hcat(map(beamform, correlated_signals)...)
    next_carrier_loop, next_carrier_freq_update = carrier_loop(beamformed_signal, Δt) #pll 
    next_code_loop, next_code_freq_update = code_loop(beamformed_signal, Δt) #dll

    (next_l1_signals, next_l5_signals, beamform, next_velocity_aiding = 0.0) -> _joined_tracking(next_l1_signals, next_l5_signals, beamform, next_l1_gen_carrier_replica, next_l5_gen_carrier_replica, next_l1_gen_code_replica, next_l5_gen_code_replica, l1_init_carrier_freq, l5_init_carrier_freq, l1_init_code_freq, l5_init_code_freq, next_carrier_freq_update, next_code_freq_update, next_carrier_loop, next_code_loop, l1_aiding_scale_factor, l5_aiding_scale_factor, next_velocity_aiding), JoinedTrackingResults(TrackingResults(l1_carrier_doppler, l1_carrier_phase, l1_code_doppler, l1_code_phase, prompt(l1_correlated_signals)), TrackingResults(l5_carrier_doppler, l5_carrier_phase, l5_code_doppler, l5_code_phase, prompt(l5_correlated_signals)))
end
