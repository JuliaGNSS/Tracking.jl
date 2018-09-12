"""
$(SIGNATURES)

Initialize loop_filter
Uses the 4 matrices `F`,`L`,`C`,`D` (transition matrix, filter gain matrix and the output matrices)
to calculate the initial state vector `x` and provide a first loop function
which takes the discriminator output `δΘ` and returns a new loop function and the system output `y`
"""
function init_loop_filter(F, L, C, D)
    x = copy(L(0.0s)) .* 0.0
    (δΘ, Δt) -> _loop_filter(x, δΘ, Δt, F, L, C, D)
end


"""
$(SIGNATURES)

loop filter
This is the loop_filter function which will be used repeatedly after being initialized by init_loop_filter
it takes the momentary state vector `x`, the discriminator output `δΘ`, and the 4 matrices `F`,`L`,`C`,`D` (transition matrix, filter gain matrix and the output matrices)
it simulates another timestep and calculates the new state Vector `x_next` and the system output `y`
Returns a new loop_filter_function with updated parameters, and the system output `y`

"""
function _loop_filter(x, δΘ, Δt, F, L, C, D)
    next_x = F(Δt) * x .+ L(Δt) * δΘ
    y = dot(C(Δt), x) + D(Δt) * δΘ
    (next_δΘ, next_Δt) -> _loop_filter(next_x, next_δΘ, next_Δt, F, L, C, D), y
end

"""
$(SIGNATURES)

Initialize a 1st order loop_filter
Takes the noise `bandwidth` and calculates the appropriate matrices for an 1st order loop_filter and hands them to init_loop_filter
Returns a 1st order loop_filter function.
"""
function init_1st_order_loop_filter(bandwidth)
    ω0 = bandwidth * 4.0
    F(Δt) = 0.0
    L(Δt) = [ω0]
    C(Δt) = [1.0]
    D(Δt) = 0.0Hz
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise `bandwidth` and the loop update time `Δt`.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter function.
"""
function init_2nd_order_boxcar_loop_filter(bandwidth)
    ω0 = bandwidth * 1.89
    F(Δt) = 1.0
    L(Δt) = [Δt * ω0^2]
    C(Δt) = [1.0]
    D(Δt) = sqrt(2) * ω0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise `bandwidth` and the loop update time `Δt`.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter function.
"""
function init_2nd_order_bilinear_loop_filter(bandwidth)
    ω0 = bandwidth * 1.89
    F(Δt) = 1.0
    L(Δt) = [Δt * ω0^2]
    C(Δt) = [1.0]
    D(Δt) = sqrt(2) * ω0 + ω0^2 * Δt / 2
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise `bandwidth` and the loop update time `Δt`.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter function.
"""
function init_3rd_order_boxcar_loop_filter(bandwidth)
    ω0 = bandwidth * 1.2
    F(Δt) = [1.0 Δt; 0.0Hz 1.0]
    L(Δt) = [Δt * 1.1 * ω0^2; Δt * ω0^3]
    C(Δt) = [1.0, 0.0s]
    D(Δt) = 2.4 * ω0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise `bandwidth` and the loop update time `Δt`.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter function.
"""
function init_3rd_order_bilinear_loop_filter(bandwidth)
    ω0 = bandwidth * 1.2
    F(Δt) = [1.0 Δt; 0.0Hz 1.0]
    L(Δt) = [Δt * 1.1 * ω0^2 + ω0^3 * Δt^2 / 2; Δt * ω0^3]
    C(Δt) = [1.0, Δt / 2]
    D(Δt) = 2.4 * ω0 + 1.1 * ω0^2 * Δt / 2 + ω0^3 * Δt^2 / 4
    init_loop_filter(F, L, C, D)
end
