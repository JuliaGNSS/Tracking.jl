

function init_PDLL(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, Delta_t, sampling_freq)
    phase = calc_phase(init_phase, init_freq, Delta_t, sampling_freq)
    init_replica = calc_signal(init_phase, init_freq, 1:Delta_t, sampling_freq)
    signal -> _PDLL(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Delta_t, sampling_freq), init_replica, phase 
end

function _PDLL(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Delta_t, sampling_freq)
    next_loop_filter, freq_update = loop_filter(disc(signal))
    replica = calc_signal(1:Delta_t, init_freq + freq_update, phase, sampling_freq)
    next_phase = calc_phase(Delta_t, init_freq + freq_update, phase, sampling_freq)
    next_signal -> _PDLL(next_signal, disc, next_loop_filter, calc_phase, calc_signal, next_phase, init_freq, Delta_t, sampling_freq), replica, phase
end

# This is just an example. Please, use GNSSSignals instead.
function gen_carrier(t, f, φ₀, f_s)
  cis((2 * π * f / f_s) .* t .+ φ₀)
end
# This is just an example. Please, use GNSSSignals instead.
function calc_carrier_phase(t, f, φ₀, f_s)
  mod2pi((2 * π * f / f_s) * t + φ₀)
end
# This is just an example. Please, use GNSSSignals instead.
function gen_sat_code(t, f_c, φ₀, f_s, code)
  c = floor.(Int, (f_c / f_s) .* t .+ φ₀)
  return code[1 .+ mod.(c, length(code))]
end
# This is just an example. Please, use GNSSSignals instead.
function get_sat_code_phase(t, f_c, φ₀, f_s, code_length)
  return mod(f_c / f_s * t + φ₀ + code_length / 2, code_length) - code_length / 2
end

function init_PLL(init_phase, init_freq, Delta_t, sampling_freq)

    gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
    sat_prn = 1
    code = gen_sampled_code(1:4000, 1023, 20, 4e6, sat_prn)
    code_phase = get_code_phase(4000, 1023, 20, 4e6)
    carrier = gen_carrier(1:4000, 1e3, 20 * pi / 180, 4e6)
    carrier_phase = get_carrier_phase(4000, 1e3, 20 * pi / 180, 4e6)

    # loop_filter = init_3rd_order_loop_filter(...)
    # disc = PLL_disc(...)
    # calc_phase = calc_carrier_phase
    # calc_signal = gen_carrier
    init_PDLL(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, Delta_t, sampling_freq)
end

function init_DLL()
    early_prompt_late_phase = [-0.5, 0, 0.5]
    # ....
    calc_signal(t, f, phase, sampling_freq) = map(phase_shift -> gen_sat_code(t, f, phase + phase_shift, sampling_freq, code), early_prompt_late_phase)
    # ...
    init_PDLL(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, Delta_t, sampling_freq)
end
