function calc_sample_shift(systems::Vector{<:AbstractGNSSSystem}, sample_freqs, preferred_phase)
    calc_sample_shift(systems[1], sample_freqs[1], preferred_phase)
end

function gen_replica_downconvert_correlate!(correlated_signals, systems::Vector{<:AbstractGNSSSystem}, signals, total_integration_samples, integrated_samples, start_samples, sample_freqs, interm_freqs, code_sample_shift, dopplers, phases, sat_prn, velocity_aidings)
    destruct(gen_replica_downconvert_correlate!.(correlated_signals, systems, signals, total_integration_samples, integrated_samples, start_samples, sample_freqs, interm_freqs, code_sample_shift, dopplers, phases, sat_prn, velocity_aidings))
end

function init_track_results(signals::Vector{<:AbstractArray}, integrated_samples, total_integration_samples)
    init_track_results.(signals, integrated_samples, total_integration_samples)
end

function aid_dopplers(systems::Vector{<:AbstractGNSSSystem}, inits, carrier_freq_update, code_freq_update, velocity_aidings)
    center_freq_mean = mean([system.center_freq for system in systems])
    code_freq_mean = mean([system.code_freq for system in systems])
    matched_carrier_freq_update = carrier_freq_update .* getfield.(systems, :center_freq) ./ center_freq_mean
    matched_code_freq_update = code_freq_update .* getfield.(systems, :code_freq) ./ code_freq_mean
    aid_dopplers.(systems, inits, matched_carrier_freq_update, matched_code_freq_update, velocity_aidings)
end

function beamform_and_update_loops(correlated_signals::Vector{Vector{Vector{ComplexF64}}}, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
    joined_correlated_signals = foldl((corrs1, corrs2) -> vcat.(corrs1, corrs2), correlated_signals)
    beamform_and_update_loops(joined_correlated_signals, Δt, beamform, carrier_loop, code_loop, actual_code_phase_shift)
end

function push_track_results!(track_results, correlated_signals::Vector{Vector{Vector{ComplexF64}}}, start_samples, max_total_integration_samples, dopplers, phases)
    push_track_results!.(track_results, correlated_signals, start_samples, max_total_integration_samples, dopplers, phases)
end

function init_correlated_signals(signals::Vector{<:AbstractArray})
    init_correlated_signals.(signals)
end

function init_phases(initials::Vector{Initials})
    init_phases.(initials)
end

function check_init_track_consistency(systems::Vector{<:AbstractGNSSSystem}, sample_freqs)
    code_sample_freq_ratio = map((system, sample_freq) -> system.code_freq / sample_freq, systems, sample_freqs)
    if !all(code_sample_freq_ratio .== code_sample_freq_ratio[1])
        error("The code sample freq ratio must be identical for all given systems.")
    end
end

function default_velocity_aiding(systems::Vector{<:AbstractGNSSSystem})
    zeros(typeof(1.0Hz), length(systems))
end

function init_dopplers(inits::Vector{Initials})
    init_dopplers.(inits)
end

function check_need_new_signal(next_start_samples::Vector, signals)
    all(next_start_samples .== size.(signals, 1) .+ 1)
end
