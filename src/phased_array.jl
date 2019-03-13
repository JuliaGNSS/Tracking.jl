"""
$(SIGNATURES)

Extends the tracking function to support phased antenna arrays.
"""
function init_tracking(
        ::Val{N},
        system,
        inits,
        sample_freq,
        interm_freq,
        sat_prn;
        pll_bandwidth = 18Hz,
        dll_bandwidth = 1Hz,
        min_integration_time = 0.5ms,
        max_integration_time = 1ms,
        carrier_loop_func = Tracking.init_3rd_order_bilinear_loop_filter,
        code_loop_func = Tracking.init_2nd_order_bilinear_loop_filter,
        cn0_update_time = 20ms
    ) where N
    code_shift = Tracking.CodeShift{3}(system, sample_freq, 0.5) # 3: Early, Prompt, Late
    dopplers = Tracking.Dopplers(inits)
    phases = Tracking.Phases(inits)
    carrier_loop = carrier_loop_func(pll_bandwidth)
    code_loop = code_loop_func(dll_bandwidth)
    correlator_outputs = init_correlator_outputs(Val(N), code_shift)
    data_bits = Tracking.DataBits(system)
    last_valid_correlator_outputs = copy(correlator_outputs)
    cn0_state = Tracking.CN0State(cn0_update_time)
    Tracking.req_signal_and_track(correlator_outputs, last_valid_correlator_outputs, system, sample_freq, interm_freq, inits, dopplers, phases, code_shift, carrier_loop, code_loop, sat_prn, min_integration_time, max_integration_time, 0, data_bits, cn0_state)
end

function init_correlator_outputs(::Val{A}, code_shift::Tracking.CodeShift{N}) where {N,A}
    zeros(SMatrix{N, A, ComplexF64})
end

Base.@propagate_inbounds function dump!(output, signal::Matrix, output_idx, sample, code_carrier)
    for ant_idx = 1:size(output, 2)
        @fastmath output[output_idx,ant_idx] += signal[sample,ant_idx] * code_carrier
    end
end

function veryearly(x::SMatrix{N,A,T}) where {N,A,T}
    x[(N - 1) >> 1 + 3,:]
end

function early(x::SMatrix{N,A,T}) where {N,A,T}
    x[(N - 1) >> 1 + 2,:]
end

function prompt(x::SMatrix{N,A,T}) where {N,A,T}
    x[(N - 1) >> 1 + 1,:]
end

function late(x::SMatrix{N,A,T}) where {N,A,T}
    x[(N - 1) >> 1,:]
end

function verylate(x::SMatrix{N,A,T}) where {N,A,T}
    x[(N - 1) >> 1 - 1,:]
end
