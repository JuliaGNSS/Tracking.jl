"""
$(SIGNATURES)

Initialize the 'locked_loop' function, calculate the replicated signal `init_replica`and its 'phase' and return the function and both values.

# Arguments
  - `disc::Function`: the DLL or PLL discriminator function 
  - `loop_filter::Function`: the loop filter function of 1st, 2nd or 3rd order 
  - `calc_phase::Function`: the phase calculation function 
  - `calc_signal::Function`: the signal replication function 
  - `init_phase::Float` the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `num_samples::Integer`: the amount of samples per carrier or satellite code 

"""
function init_locked_loop(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, sampling_freq, num_samples)
    phase = calc_phase(num_samples, init_freq, init_phase, sampling_freq)
    init_replica = calc_signal(1:num_samples, init_freq, init_phase, sampling_freq)
    (signal, aiding = 0) -> _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, sampling_freq, num_samples, aiding), init_replica, phase
end


"""
$(SIGNATURES)

Calculate the replication_signal `replica`, the replication signals phase `next_phase`, and the loop filter function `next_loop_filter` for the next timestep and return them.

# Arguments
  - `signal::Array{Float,num_samples}`: the replicated signal 
  - `disc::Function`: the DLL or PLL discriminator function 
  - `loop_filter::Function`: the loop filter function of 1st, 2nd or 3rd order 
  - `calc_phase::Function`: the phase calculation function 
  - `calc_signal::Function`: the signal replication function 
  - `init_phase::Float`: the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `num_samples::Integer`: the amount of samples per carrier or satellite code 
  - `aiding`: the aiding provided by the PLL, 0 if this is the locked_loop of an PLL

"""
function _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, sampling_freq, num_samples, aiding)
    next_loop_filter, freq_update = loop_filter(disc(signal))
    next_freq = init_freq + freq_update + aiding
    replica = calc_signal(1:num_samples, next_freq, phase, sampling_freq)
    next_phase = calc_phase(num_samples, next_freq, phase, sampling_freq)
    (next_signal, next_aiding = 0) -> _locked_loop(next_signal, disc, next_loop_filter, calc_phase, calc_signal, next_phase, init_freq, sampling_freq, num_samples, next_aiding), replica, phase, next_freq
end


"""
$(SIGNATURES)

Call the initialization of a PLL locked loop and return the replicated signal `init_replica`, the calculated phase `phase`, and a pll_locked_loop function which can be used to calculate the consecutive values.

# Arguments
  - `init_phase::Float`: the initial signal phase in rad
  - `init_freq::Float`: the initial signal frequency in Hz
  - `Δt::Float`: the loop update time intervall in seconds
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `bandwidth::Float`: the signal aquivalent noise bandwidth for the loop_filter, in Hz

  # Examples
  ```julia-repl
  julia> #get correlator_output from correlator
  julia> PLL, sample_code, phase = Tracking.init_PLL(20.0, 50, 4000, 4e6, 18.0, 1e-3)
  julia> next_PLL, next_sampled_carrier, next_phase = PLL(correlator_output)
  ```

"""
function init_PLL(init_phase, init_freq, sampling_freq, bandwidth, Δt) 
  num_samples = sampling_freq * Δt
  init_locked_loop(pll_disc, init_3rd_order_loop_filter(bandwidth, Δt), GNSSSignals.get_carrier_phase, GNSSSignals.gen_carrier, init_phase, init_freq, sampling_freq, num_samples)
end


"""
$(SIGNATURES)

Initialize a dll_locked_loop function by calculating the needed parameters and return the replicated signal `init_replica`, the calculated phase `phase`, and a dll_locked_loop function which can be used to calculate the consecutive values.

# Arguments
  - `init_phase::Float` the initial signal offset in chip
  - `init_freq::Float`: the initial signal frequency in Hz
  - `Δt::Float`: the loop update time intervall in seconds
  - `sampling_freq::Float`: the signal sampling frequency in Hz
  - `bandwidth::Float`: the signal aquivalent noise bandwidth for the loop_filter, in Hz
  - `sat_prn::Integer`:the satellite PRN code number, choose from 1-32

  # Examples
    ```julia-repl
    julia> #get correlator_output from correlator
    julia> DLL, sample_code, phase = Tracking.init_DLL(20, 1023e3, 4000, 4e6, 1, 1e-3, 1)
    julia> next_DLL, next_sampled_code, next_phase = DLL(correlator_output);
    ```

"""
function init_DLL(init_phase, init_freq, sampling_freq, bandwidth, Δt, sat_prn)
  early_prompt_late_phase = [-0.5, 0, 0.5]
  gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
  calc_signal(samples, f, phase, sampling_freq) = map(phase_shift -> gen_sampled_code(samples, f, phase + phase_shift, sampling_freq, sat_prn), early_prompt_late_phase)
  loop_filter = init_2nd_order_loop_filter(bandwidth, Δt)
  num_samples = sampling_freq * Δt
  init_locked_loop(dll_disc, loop_filter, get_code_phase, calc_signal, init_phase, init_freq, sampling_freq, num_samples)
end

