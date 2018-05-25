

function init_locked_loop(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, Δt, sampling_freq)
    phase = calc_phase(Δt, init_freq, init_phase, sampling_freq)
    init_replica = calc_signal(1:Δt, init_freq, init_phase, sampling_freq)
    signal -> _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Δt, sampling_freq), init_replica, phase 
end

function _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Δt, sampling_freq)
    next_loop_filter, freq_update = loop_filter(disc(signal))
    replica = calc_signal(1:Δt, init_freq + freq_update, phase, sampling_freq)
    next_phase = calc_phase(Δt, init_freq + freq_update, phase, sampling_freq)
    next_signal -> _locked_loop(next_signal, disc, next_loop_filter, calc_phase, calc_signal, next_phase, init_freq, Δt, sampling_freq), replica, phase
end

function init_PLL(init_phase, init_freq, Δt, sampling_freq, bandwidth) 
  init_locked_loop(pll_disc, init_3rd_order_loop_filter(bandwidth ,Δt), GNSSSignals.get_carrier_phase, GNSSSignals.gen_carrier, init_phase, init_freq, Δt, sampling_freq)
end

function init_DLL(init_phase, init_freq, Δt, sampling_freq, sat_prn, bandwidth)
  early_prompt_late_phase = [-0.5, 0, 0.5]
  gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
  code = gen_sampled_code(1:Δt, init_freq, init_phase, sampling_freq, sat_prn)
  code_phase = get_code_phase(Δt, init_freq, init_phase, sampling_freq)
  calc_signal(t, f, phase, sampling_freq) = map(phase_shift -> GNSSSignals.gen_sat_code(t, f, code_phase + phase_shift, sampling_freq, code), early_prompt_late_phase)
  loop_filter = init_3rd_order_loop_filter(bandwidth ,Δt)
  init_locked_loop(dll_disc, loop_filter, get_code_phase, calc_signal, init_phase, init_freq, Δt, sampling_freq)
end
