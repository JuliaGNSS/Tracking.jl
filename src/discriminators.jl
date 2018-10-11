

early(x) = x[3]
prompt(x) = x[2]
late(x) = x[1]

"""
$(SIGNATURES)

DLL discriminator
Returns the calculated code offset in chips.
"""
function dll_disc(x, d = 1)
    E = abs(early(x))
    L = abs(late(x))
    (E - L) / (E + L) / 2 / (2 - d)
end


"""
$(SIGNATURES)

PLL discriminator
Returns the calculated Phase error in rad.
"""
function pll_disc(x)
    p = prompt(x)
    atan(imag(p) / real(p))
end
