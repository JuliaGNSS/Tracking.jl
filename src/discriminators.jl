

early(x) = x[1]
prompt(x) = x[2]
late(x) = x[3]

function dll_disc(x)
    E = abs(early(x))
    L = abs(late(x))
    (E - L) / (E + L)
end

function pll_disc(x)
    p = prompt(x)
    atan(imag(p) / real(p))
end

