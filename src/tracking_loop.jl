"""
$(SIGNATURES)

Initializes the tracking function. Returns a tracking function to track a given signal.
The initialization needs the GNSS system `system`, the doppler and phase of the
carrier and code `inits` of type `Initials`, the sample frequency `sample_freq`
of type `Unitful.Hz`, the intermediate frequency of type `Unitful.Hz` and the
satellite PRN number `prn`.
Optianally you can set the PLL bandwidth (default: 18Hz), the DLL bandwidth
(default: 1Hz), the minimal integration time which will lead to a valid
correlation results (default: 0.5ms), the maximal integration time (default: 1ms),
the carrier loop function (default: `init_3rd_order_bilinear_loop_filter`), the
the code loop function (default: `init_2nd_order_bilinear_loop_filter`) the
Carrier-to-Noise-Density-Ratio (CN0) update time `cn0_update_time`.
The tracking function will only integrate in the boundaries of a bit shifts.
Moreover, it will track the complete given signal. Adjust the number of samples
in the signal to your preferred output update rate. The data bits are bufferd in
an UInt variable. It can hold up to 64 or 32 bits depending on your system.
The number of detected bits will be returned, too.
The returned tracking function will return a new tracking function for the next
iteration and the tracking results `TrackingResults`.
# Examples
```julia-repl
julia> gpsl1 = GPSL1()
julia> inits = Initials(gpsl1, 1200Hz, 231) # 1200 Hz doppler, 231 code phase
julia> track = init_tracking(gpsl1, inits, 5Mhz, 0Hz, 1)
julia> next_track, track_results = track(signal)
```
"""
function init_tracking(
        system,
        inits,
        sample_freq,
        interm_freq,
        sat_prn;
        num_ants = NumAnts(1),
        pll_bandwidth = 18Hz,
        dll_bandwidth = 1Hz,
        min_integration_time = 0.5ms,
        max_integration_time = 1ms,
        carrier_loop_func = init_3rd_order_bilinear_loop_filter,
        code_loop_func = init_2nd_order_bilinear_loop_filter,
        cn0_update_time = 20ms
    )
    code_shift = CodeShift{3}(system, sample_freq, 0.5) # 3: Early, Prompt, Late; should later depend on system
    dopplers = Dopplers(inits)
    phases = Phases(inits)
    carrier_loop = carrier_loop_func(pll_bandwidth)
    code_loop = code_loop_func(dll_bandwidth)
    correlator_outputs = init_correlator_outputs(num_ants, code_shift)
    data_bits = DataBits(system)
    last_valid_correlator_outputs = copy(correlator_outputs)
    last_valid_filtered_correlator_outputs = init_correlator_outputs(NumAnts(1), code_shift)
    cn0_state = CN0State(cn0_update_time)
    req_signal_and_track(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, min_integration_time, max_integration_time, 0, data_bits, cn0_state)
end

"""
$(SIGNATURES)

Running tracking function. Will return a new tracking function for the next
iteration and the tracking results.
"""
function _tracking(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, signal, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, post_corr_filter, min_integration_time, max_integration_time, signal_idx, integrated_samples, num_integrated_prns, data_bits, cn0_state, velocity_aiding)
    preferred_integration_time = calc_integration_time(data_bits, max_integration_time)
    num_samples_left_to_integrate = calc_num_samples_left_to_integrate(system, sample_freq, phases, preferred_integration_time)
    num_samples_signal_bound = calc_num_samples_signal_bound(signal, signal_idx)
    num_samples_to_integrate = min(num_samples_left_to_integrate, num_samples_signal_bound)
    correlator_outputs = correlate_and_dump(correlator_outputs, signal, system, sample_freq, interm_freq, dopplers, phases, code_shift, signal_idx, num_samples_to_integrate, sat_prn)
    integrated_samples, signal_idx = (integrated_samples, signal_idx) .+ num_samples_to_integrate
    phases = calc_next_phases(system, interm_freq, sample_freq, dopplers, phases, num_samples_to_integrate, data_bits)
    actual_integration_time = calc_actual_integration_time(integrated_samples, sample_freq)
    if num_samples_to_integrate == num_samples_left_to_integrate
        num_integrated_prns += calc_integrated_prns(system, integrated_samples, sample_freq)
        if actual_integration_time >= min_integration_time
            last_valid_correlator_outputs = correlator_outputs
            filtered_correlator_outputs = post_corr_filter(correlator_outputs)
            last_valid_filtered_correlator_outputs = filtered_correlator_outputs
            carrier_loop, carrier_freq_update = carrier_loop(pll_disc(filtered_correlator_outputs), actual_integration_time)
            code_loop, code_freq_update = code_loop(dll_disc(filtered_correlator_outputs, 2 * code_shift.actual_shift), actual_integration_time)
            dopplers = aid_dopplers(system, inits, carrier_freq_update, code_freq_update, velocity_aiding)
            data_bits = buffer(data_bits, system, real(prompt(filtered_correlator_outputs)), num_integrated_prns)
            cn0_state = estimate_CN0(cn0_state, actual_integration_time, prompt(filtered_correlator_outputs))
        end
        correlator_outputs = zeros(typeof(correlator_outputs))
        integrated_samples = zero(integrated_samples)
    end
    if num_samples_to_integrate == num_samples_signal_bound
        track_results = TrackingResults(dopplers, phases, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, data_bits, num_integrated_prns, cn0_state.last_valid_cn0)
        return req_signal_and_track(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, min_integration_time, max_integration_time, integrated_samples, data_bits, cn0_state), track_results
    else
        _tracking(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, signal, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, post_corr_filter, min_integration_time, max_integration_time, signal_idx, integrated_samples, num_integrated_prns, data_bits, cn0_state, velocity_aiding)
    end
end

"""
$(SIGNATURES)

Requires a new signal and returns a new tracking function for next iteration.
"""
function req_signal_and_track(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, min_integration_time, max_integration_time, integrated_samples, data_bits, cn0_state)
    data_bits = DataBits(data_bits, 0, 0)
    (signal, post_corr_filter = x -> x, velocity_aiding = 0.0Hz) ->
        _tracking(correlator_outputs, last_valid_correlator_outputs, last_valid_filtered_correlator_outputs, signal, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, post_corr_filter, min_integration_time, max_integration_time, 1, integrated_samples, 0, data_bits, cn0_state, velocity_aiding)
end

"""
$(SIGNATURES)

Downconverts the signal and correlates with the code replica.
"""
function correlate_and_dump(correlator_outputs, signal, system, sample_freq, interm_freq, dopplers, phases, code_shift, start_sample, num_samples_to_integrate, sat_prn)
    gen_carrier_replica(x) = gen_carrier(x, -(interm_freq + dopplers.carrier), -phases.carrier, sample_freq)
    calc_code_replica_phase_unsafe(x) = calc_code_phase_unsafe(x, system.code_freq + dopplers.code, phases.code, sample_freq)
    downconvert_and_correlate(correlator_outputs, signal, system, start_sample, num_samples_to_integrate, gen_carrier_replica, calc_code_replica_phase_unsafe, code_shift, sat_prn)
end

"""
$(SIGNATURES)

Calculates the phases for the next iterations. The phase gets adjusted if no
data bit has been found yet.
"""
function calc_next_phases(system, interm_freq, sample_freq, dopplers, phases, num_samples_to_integrate, data_bits)
    next_carrier_phase = calc_carrier_phase(num_samples_to_integrate, interm_freq + dopplers.carrier, phases.carrier, sample_freq)
    next_code_phase = calc_code_phase(num_samples_to_integrate, system.code_freq + dopplers.code, phases.code, sample_freq, size(system.codes, 1))
    Phases(next_carrier_phase, adjust_code_phase(system, data_bits, next_code_phase))
end

"""
$(SIGNATURES)

Calculates the number of samples left to integrate until a bit shift might occur.
"""
function calc_num_samples_left_to_integrate(system, sample_freq, phases, integration_time)
    chips_left = integration_time * system.code_freq - mod(phases.code, convert(Int, integration_time * system.code_freq))
    ceil(Int, chips_left * sample_freq / system.code_freq)
end

"""
$(SIGNATURES)

Calculates the number of samples left to integrate until signal is running out
of samples.
"""
function calc_num_samples_signal_bound(signal, signal_idx)
    size(signal, 1) - signal_idx + 1
end

"""
$(SIGNATURES)

Calculates the number of PRNs that were integrated
"""
function calc_integrated_prns(system, integrated_samples, sample_freq)
    ceil(Int, integrated_samples / (sample_freq / system.code_freq * system.code_length))
end

"""
$(SIGNATURES)

Calculates the necessary integration time. If the data bit has not been found yet,
it will default to the hard coded 1ms. Otherwise it will return the specified
maximum integration time.
"""
function calc_integration_time(data_bits, max_integration_time)
    ifelse(found(data_bits), max_integration_time, 1ms)
end

"""
$(SIGNATURES)

Initializes the correlator outputs.
"""
function init_correlator_outputs(::NumAnts{1}, code_shift::CodeShift{N}) where N
    zeros(SVector{N, ComplexF64})
end

"""
$(SIGNATURES)

Aid dopplers. That is velocity aiding for the carrier doppler and carrier aiding
for the code doppler.
"""
function aid_dopplers(system, inits, carrier_freq_update, code_freq_update, velocity_aiding)
    carrier_doppler = carrier_freq_update + velocity_aiding
    code_doppler = code_freq_update + carrier_doppler * system.code_freq / system.center_freq
    Dopplers(inits.carrier_doppler + carrier_doppler, inits.code_doppler + code_doppler)
end

"""
$(SIGNATURES)

Calculates that actual integration time.
"""
function calc_actual_integration_time(integrated_samples, sample_freq)
    integrated_samples / sample_freq
end

"""
$(SIGNATURES)

Downconverts and correlates the signal.
"""
function downconvert_and_correlate(output, signal, system, start_sample, num_samples_to_integrate, gen_carrier_replica, calc_code_replica_phase_unsafe, code_shift, prn)
    mutual_output = MArray(output)
    shifts = init_shifts(code_shift)
    code_idx_wrap = 0
    @inbounds for sample_idx = 1:num_samples_to_integrate
        carrier = gen_carrier_replica(sample_idx)
        for (output_idx, code_sample_shift) in enumerate(shifts)
            code_idx = floor(Int, calc_code_replica_phase_unsafe(code_sample_shift + sample_idx)) + code_idx_wrap
            wrapped_code_idx, code_idx_wrap = wrap_code_idx(system, code_idx, code_idx_wrap)
            code_carrier = system.codes[wrapped_code_idx + 1,prn] * carrier
            dump!(mutual_output, signal, output_idx, sample_idx + start_sample - 1, code_carrier)
        end
    end
    SArray(mutual_output)
end

"""
$(SIGNATURES)

Wraps the code index with the maximal code index.
"""
function wrap_code_idx(system, code_idx, code_idx_wrap)
    if code_idx >= system.code_length
        code_idx_wrap -= system.code_length
        code_idx -= system.code_length
    elseif code_idx <= -1
        code_idx += system.code_length
    end
    code_idx, code_idx_wrap
end

"""
$(SIGNATURES)

Initializes the code shift register
"""
function init_shifts(code_shift::CodeShift{N}) where N
    min_max_code_shift_samples = (N - 1) / 2 * code_shift.samples
    SVector{N,Int}(-min_max_code_shift_samples:code_shift.samples:min_max_code_shift_samples)
end

"""
$(SIGNATURES)

Dumps the downconverted and correlated signal for the case that signal is a
vector.
"""
Base.@propagate_inbounds function dump!(output, signal::Vector, output_idx, sample, code_carrier)
    @fastmath output[output_idx] += signal[sample] * code_carrier
end

"""
$(SIGNATURES)

General code phase adjustment. Should be overriden by specified system, if
necessary.
"""
function adjust_code_phase(system, data_bits, phase)
    phase
end
