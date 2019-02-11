"""
$(SIGNATURES)

Initialize loop_filter
Uses the 4 matrices `F`,`L`,`C`,`D` (transition matrix, filter gain matrix and
the output matrices) to calculate the initial state vector `x` and returns a
loop function which takes the discriminator output `δΘ` and returns a new loop
function and the system output `y`
"""
function init_loop_filter(F, L, C, D)
    x = zero(L(0.0))
    req_error_and_filter(x, F, L, C, D)
end

function req_error_and_filter(x, F, L, C, D)
    (δΘ, Δt) -> _loop_filter(x, δΘ, Δt, F, L, C, D)
end


"""
$(SIGNATURES)

Internal loop filter to propagate the state and return the loop filter output.
"""
function _loop_filter(x, δΘ, Δt, F, L, C, D)
    Δt_in_sec = Float64(upreferred(Δt/s))
    next_x = F(Δt_in_sec) * x .+ L(Δt_in_sec) * δΘ
    y = dot(C(Δt_in_sec), x) + D(Δt_in_sec) * δΘ # next_x or x? Michael: next_x
    req_error_and_filter(next_x, F, L, C, D), y * Hz
end

"""
$(SIGNATURES)

Initialize a 1st order loop filter with noise bandwidth `bandwidth`. Returns a
loop filter which takes the discriminator output.
"""
function init_1st_order_loop_filter(bandwidth)
    ω0 = Float64(bandwidth/Hz) * 4.0
    F(Δt) = 0.0
    L(Δt) = ω0
    C(Δt) = 1.0
    D(Δt) = 0.0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order boxcar loop filter with noise bandwidth `bandwidth`.
Returns a loop filter which takes the discriminator output.
"""
function init_2nd_order_boxcar_loop_filter(bandwidth)
    ω0 = Float64(bandwidth/Hz) * 1.89
    F(Δt) = 1.0
    L(Δt) = Δt * ω0^2
    C(Δt) = 1.0
    D(Δt) = sqrt(2) * ω0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order bilinear loop filter with noise bandwidth `bandwidth`.
Returns a loop filter which takes the discriminator output.
"""
function init_2nd_order_bilinear_loop_filter(bandwidth)
    ω0 = Float64(bandwidth/Hz) * 1.89
    F(Δt) = 1.0
    L(Δt) = Δt * ω0^2
    C(Δt) = 1.0
    D(Δt) = sqrt(2) * ω0 + ω0^2 * Δt / 2
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 3rd order boxcar loop filter with noise bandwidth `bandwidth`.
Returns a loop filter which takes the discriminator output.
"""
function init_3rd_order_boxcar_loop_filter(bandwidth)
    ω0 = Float64(bandwidth/Hz) * 1.2
    F(Δt) = @SMatrix [1.0 Δt; 0.0 1.0]
    L(Δt) = @SVector [Δt * 1.1 * ω0^2, Δt * ω0^3]
    C(Δt) = @SVector [1.0, 0.0]
    D(Δt) = 2.4 * ω0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 3rd order bilinear loop filter with noise bandwidth `bandwidth`.
Returns a loop filter which takes the discriminator output.
"""
function init_3rd_order_bilinear_loop_filter(bandwidth)
    ω0 = Float64(bandwidth/Hz) * 1.2
    F(Δt) = @SMatrix [1.0 Δt; 0.0 1.0]
    L(Δt) = @SVector [Δt * 1.1 * ω0^2 + ω0^3 * Δt^2 / 2, Δt * ω0^3]
    C(Δt) = @SVector [1.0, Δt / 2]
    D(Δt) = 2.4 * ω0 + 1.1 * ω0^2 * Δt / 2 + ω0^3 * Δt^2 / 4
    init_loop_filter(F, L, C, D)
end
