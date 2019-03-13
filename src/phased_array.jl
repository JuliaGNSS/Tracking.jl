function init_correlator_outputs(::NumAnts{A}, code_shift::Tracking.CodeShift{N}) where {N,A}
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
