

early(x) = x[1]
prompt(x) = x[2]
late(x) = x[3]

"""
$(SIGNATURES)

DLL Diskriminator
Returns the calculated code offset in chips.
"""
function dll_disc(x)
    E = abs(early(x))
    L = abs(late(x))
    (E - L) / (E + L)
end


"""
$(SIGNATURES)

PLL Diskriminator
Returns the calculated Phase error in rad.
"""
function pll_disc(x)
    p = prompt(x)
    atan(imag(p) / real(p))
end

