"""
$(SIGNATURES)

Returns the very early
"""
function veryearly(x::SVector{N,T}) where {N,T}
    x[(N - 1) >> 1 + 3]
end

"""
$(SIGNATURES)

Returns the early
"""
function early(x::SVector{N,T}) where {N,T}
    x[(N - 1) >> 1 + 2]
end

"""
$(SIGNATURES)

Returns the prompt
"""
function prompt(x::SVector{N,T}) where {N,T}
    x[(N - 1) >> 1 + 1]
end

"""
$(SIGNATURES)

Returns the late
"""
function late(x::SVector{N,T}) where {N,T}
    x[(N - 1) >> 1]
end

"""
$(SIGNATURES)

Returns the very late
"""
function verylate(x::SVector{N,T}) where {N,T}
    x[(N - 1) >> 1 - 1]
end

"""
$(SIGNATURES)

Calculates the code offset in chips.
"""
function dll_disc(x, d = 1)
    E = abs(early(x))
    L = abs(late(x))
    (E - L) / (E + L) / 2 / (2 - d)
end


"""
$(SIGNATURES)

Calculates the phase error in radians.
"""
function pll_disc(x)
    p = prompt(x)
    atan(imag(p) / real(p))
end
