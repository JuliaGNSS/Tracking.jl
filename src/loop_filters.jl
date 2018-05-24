function init_loop_filter(F, L, C, D)
    x = zeros(size(F,1))
    δΘ -> _loop_filter(x, δΘ, F, L, C, D)
end

function _loop_filter(x, δΘ, F, L, C, D)
    x_next = F * x + L * δΘ
    y = C * x + D * δΘ
    δΘ_next -> _loop_filter(x_next, δΘ_next, F, L, C, D), y
end

function init_1st_order_loop_filter(bandwidth)
    ω0 = bandwidth * 4
    F = 0
    L = ω0
    C = 1
    D = 0
    init_loop_filter(F, L, C, D)
end

function init_2nd_order_loop_filter(bandwidth,Δt)
    ω0 = bandwidth *1.89
    F = 1
    L = Δt*ω0^2
    C = 1
    D = sqrt(2)*ω0
    init_loop_filter(F, L, C, D)
end

function init_3rd_order_loop_filter(bandwidth,Δt)
    ω0 = bandwidth*1.2
    F = [1 Δt; 0 1]
    L = [Δt*1.1*ω0^2; Δt*ω0^3]
    C = [1 0]
    D = [2.4*ω0]
    init_loop_filter(F, L, C, D)
end


