
"""
$(SIGNATURES)

init_locked_loop(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, sampling_freq, n_samples )

Initialize the 'locked_loop' function, calculate the replicated signal `init_replica`and its 'phase' and return all three values.

# Arguments
  - `disc::Function`: the DLL or PLL discriminator function 
  - `loop_filter::Function`: the loop filter function of 1st, 2nd or 3rd order 
  - `calc_phase::Function`: the phase calculation function 
  - `calc_signal::Function`: the signal replication function 
  - `init_phase::Float` the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `n_samples::Integer`: the amount of chips per carrier or satellite code 

```

"""
function init_locked_loop(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, sampling_freq, n_samples )
    phase = calc_phase(n_samples , init_freq, init_phase, sampling_freq)
    init_replica = calc_signal(1:n_samples , init_freq, init_phase, sampling_freq)
    signal -> _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, sampling_freq, n_samples ), init_replica, phase 
end


"""
$(SIGNATURES)
_locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, sampling_freq, n_samples )

Calculate the replication_signal `replica`, the replication signals phase `next_phase`, and the loop filter function `next_loop_filter` for the next timestep and return them.

# Arguments
  - `signal::Array{Float,n_samples}`: the replicated signal 
  - `disc::Function`: the DLL or PLL discriminator function 
  - `loop_filter::Function`: the loop filter function of 1st, 2nd or 3rd order 
  - `calc_phase::Function`: the phase calculation function 
  - `calc_signal::Function`: the signal replication function 
  - `init_phase::Float`: the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `n_samples::Integer`: the amount of chips per carrier or satellite code 


# Examples
"""
function _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, sampling_freq, n_samples )
    next_loop_filter, freq_update = loop_filter(disc(signal))
    replica = calc_signal(1:n_samples , init_freq + freq_update[1], phase, sampling_freq)
    next_phase = calc_phase(n_samples , init_freq + freq_update[1], phase, sampling_freq)
    next_signal -> _locked_loop(next_signal, disc, next_loop_filter, calc_phase, calc_signal, next_phase, init_freq, sampling_freq, n_samples ), replica, phase
end


"""
$(SIGNATURES)

init_PLL(init_phase, init_freq, n_samples, sampling_freq, bandwidth, Δt) 

Call the initialization of a PLL locked loop and return the `locked_loop` function, the replicated signal `init_replica`, and the calculated phase `phase`.

# Arguments
  - `init_phase::Float`: the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `Δt::Float`: the loop update time intervall in seconds
  - `sampling_freq::Float`: the signal samlping frequency in Hz
  - `bandwidth::Float`: the signal aquivalent noise bandwidth  for the loop_filter, in Hz
  - `n_samples::Integer`: the amount of chips per carrier or satellite code 

  # Examples
  locked_loop, sample_code, phase = Tracking.init_PLL(20.0, 50, 4000, 4e6, 18.0, 1e-3)
"""
function init_PLL(init_phase, init_freq, n_samples, sampling_freq, bandwidth, Δt) 
  init_locked_loop(pll_disc, init_3rd_order_loop_filter(bandwidth ,Δt), GNSSSignals.get_carrier_phase, GNSSSignals.gen_carrier, init_phase, init_freq, sampling_freq, n_samples )
end


"""
$(SIGNATURES)

init_DLL(init_phase, init_freq, n_samples, sampling_freq, bandwidth, sat_prn, Δt )

Initialize a dll_locked_loop function by calculating the needed parameters and return the locked_loop function, the replicated signal `init_replica`, and the calculated phase `phase`.
  
To initialize the locked_loop folling has to b calculated
  - generates a sampled satellite code `code`
  - calculates the codes phase `code_phase`
  - creates a `calc_signal` function
  - creates a 3rd order loop filter function loop_filter

# Arguments
  - `init_phase::Float` the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `Δt::Float`: the loop update time intervall in seconds
  - `sampling_freq::Float`: the signal samlping frequency in Hz
  - `bandwidth::Float`: the signal aquivalent noise bandwidth  for the loop_filter, in Hz
  - `n_samples::Integer`: the amount of chips per carrier or satellite code 
  - `sat_prn::Integer`:the satellite PRN code number, choose from 1-32

  # Examples
  
  locked_loop, sample_code, phase = Tracking.init_DLL(20, 1023e3, 4000, 4e6, 1, 1e-3, 1)

"""
function init_DLL(init_phase, init_freq, n_samples, sampling_freq, bandwidth,  Δt, sat_prn )
  early_prompt_late_phase = [-0.5, 0, 0.5]
  gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
  code = gen_sampled_code(1:n_samples, init_freq, init_phase, sampling_freq, sat_prn)
  code_phase = get_code_phase(n_samples, init_freq, init_phase, sampling_freq)
  calc_signal(t, f, phase, sampling_freq) = map(phase_shift -> GNSSSignals.gen_sat_code(t, f, code_phase + phase_shift, sampling_freq, code), early_prompt_late_phase)
  loop_filter = init_3rd_order_loop_filter(bandwidth ,Δt)
  init_locked_loop(dll_disc, loop_filter, get_code_phase, calc_signal, init_phase, init_freq, sampling_freq, n_samples)
end

