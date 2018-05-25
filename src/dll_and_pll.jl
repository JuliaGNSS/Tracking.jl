
"""
$(SIGNATURES)

Initialize the locked_loop function
Takes several inputs to initialize the _locked_loop function:
  - the DLL or PLL discriminator function `disc`
  - the loop filter function of 1st, 2nd or 3rd order `loop_filter`
  - the phase calculation function `calc_phase`
  - the signal replication function `calc_signal`
  - the initial signal phase `init_phase`
  - the initial signal frequency `init_freq`
  - the loop update time intervall`Δt`
  - the signal samlping frequency `sampling_freq`

With the provided inputs and functions the replicated signal and the phase are calculated.

Returns a locked_loop function, the replicated signal `init_replica`, and the calculated phase `phase` 

"""
function init_locked_loop(disc, loop_filter, calc_phase, calc_signal, init_phase, init_freq, Δt, sampling_freq)
    phase = calc_phase(Δt, init_freq, init_phase, sampling_freq)
    init_replica = calc_signal(1:Δt, init_freq, init_phase, sampling_freq)
    signal -> _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Δt, sampling_freq), init_replica, phase 
end


"""
$(SIGNATURES)
Locked Loop function
Takes several inputs :
  - the replicated signal `signal`
  - the DLL or PLL discriminator function `disc`
  - the loop filter function of 1st, 2nd or 3rd order `loop_filter`
  - the phase calculation function `calc_phase`
  - the signal replication function `calc_signal`
  - the calculated signal phase `phase`
  - the initial signal frequency `init_freq`
  - the loop update time intervall`Δt`
  - the signal samlping frequency `sampling_freq`

Updates its own inputs for the next time step:
  - the replication_signal `replica`
  - the replication signals phase `next_phase`
  - the loop filter `next_loop_filter`

Returns the updated _locked_loop function, the updated replicated signal `replica`, and the new calculated phase `phase` 
"""
function _locked_loop(signal, disc, loop_filter, calc_phase, calc_signal, phase, init_freq, Δt, sampling_freq)
    next_loop_filter, freq_update = loop_filter(disc(signal))
    replica = calc_signal(1:Δt, init_freq + freq_update, phase, sampling_freq)
    next_phase = calc_phase(Δt, init_freq + freq_update, phase, sampling_freq)
    next_signal -> _locked_loop(next_signal, disc, next_loop_filter, calc_phase, calc_signal, next_phase, init_freq, Δt, sampling_freq), replica, phase
end


"""
$(SIGNATURES)

initialize PLL
Takes several inputs:
  - the initial signal phase `init_phase`
  - the initial signal frequency `init_freq`
  - the loop update time intervall`Δt`
  - the signal samlping frequency `sampling_freq`
  - the signal aquivalent the noise bandwidth `bandwidth`

Calls the initialization of a PLL locked loop with an appropriate discriminator `disc`
By calling init_locked_loop the following values are returned:
 - a locked_loop function
 - the replicated signal `init_replica`
 - the calculated phase `phase` 

"""
function init_PLL(init_phase, init_freq, Δt, sampling_freq, bandwidth) 
  init_locked_loop(pll_disc, init_3rd_order_loop_filter(bandwidth ,Δt), GNSSSignals.get_carrier_phase, GNSSSignals.gen_carrier, init_phase, init_freq, Δt, sampling_freq)
end


"""
$(SIGNATURES)

initialize DLL
Takes several inputs:
  - the initial signal phase `init_phase`
  - the initial signal frequency `init_freq`
  - the loop update time intervall`Δt`
  - the signal samlping frequency `sampling_freq`
  - the signal aquivalent the noise bandwidth `bandwidth`
  - the satellite PRN code number `sat_prn`

Uses functions from the GNSSSignals module.

Calls the initialization of a DLL locked loop after calculating the needed parameters:
  - generates a sampled satellite code `code`
  - calculates the codes phase `code_phase`
  - creates a `calc_signal` function
  - creates a 3rd order loop filter function loop_filter
  - 

By calling init_locked_loop the following values are returned:
 - a locked_loop function
 - the replicated signal `init_replica`
 - the calculated phase `phase` 

"""
function init_DLL(init_phase, init_freq, Δt, sampling_freq, bandwidth, sat_prn)
  early_prompt_late_phase = [-0.5, 0, 0.5]
  gen_sampled_code, get_code_phase = GNSSSignals.init_gpsl1_codes()
  code = gen_sampled_code(1:Δt, init_freq, init_phase, sampling_freq, sat_prn)
  code_phase = get_code_phase(Δt, init_freq, init_phase, sampling_freq)
  calc_signal(t, f, phase, sampling_freq) = map(phase_shift -> GNSSSignals.gen_sat_code(t, f, code_phase + phase_shift, sampling_freq, code), early_prompt_late_phase)
  loop_filter = init_3rd_order_loop_filter(bandwidth ,Δt)
  init_locked_loop(dll_disc, loop_filter, get_code_phase, calc_signal, init_phase, init_freq, Δt, sampling_freq)
end
