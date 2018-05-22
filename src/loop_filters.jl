function init_loop_filter(F, L, C, D)
    x = zeros(size(F,1))
    delta_theta -> _loop_filter(x, delta_theta, F, L, C, D)
end

function _loop_filter(x, delta_theta, F, L, C, D)
    x_next = F * x + L * delta_theta
    y = C * x + D * delta_theta
    delta_theta_next -> _loop_filter(x_next, delta_theta_next, F, L, C, D), y
end

function init_1st_order_loop_filter(bandwidth)
    # F = ...
    # L = ...
    # ...
        init_loop_filter(F, L, C, D)
    end

function init_2nd_order_loop_filter(bandwidth)
# F = ...
# L = ...
# ...
    init_loop_filter(F, L, C, D)
end

function init_3rd_order_loop_filter(bandwidth)
# F = ...
# L = ...
# ...
    init_loop_filter(F, L, C, D)
end


