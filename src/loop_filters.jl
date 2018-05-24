"""
$(SIGNATURES)

initialize loop_filter
Uses the 4 matrices F,L,C,D (transition matrix, filter gain matrix and the output matrices )
to calculate the initial state vector x and provide a first loop function
which takes the discriminator output δΘ and returns a new loop function and the system output y
"""
function init_loop_filter(F, L, C, D)
    x = zeros(size(F,1))
    δΘ -> _loop_filter(x, δΘ, F, L, C, D)
end


"""
$(SIGNATURES)

loop filter
This is the loop_filter function which will be used repeatedly after beeing initialized by init_loop_filter
it takes the momentary state vector x, the discriminator output δΘ, and the 4 matrices F,L,C,D (transition matrix, filter gain matrix and the output matrices )
it simulates another timestep and calculates the new state Vector x_next and the system output y
Returns a new loop_filter_function with updated parameters, and the system output y 

"""
function _loop_filter(x, δΘ, F, L, C, D)
    x_next = F * x + L * δΘ
    y = C * x + D * δΘ
    δΘ_next -> _loop_filter(x_next, δΘ_next, F, L, C, D), y
end

"""
$(SIGNATURES)

Initialize a 1st order loop_filter
Takes the noise bandwidth and calculates the appropriate matrices for an 1st order loop_filter and hands them to init_loop_filter
Returns a 1st order loop_filter function.
"""
function init_1st_order_loop_filter(bandwidth)
    ω0 = bandwidth * 4
    F = 0
    L = ω0
    C = 1
    D = 0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise bandwidth and the loop update time Δt.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter fucntion.
"""
function init_2nd_order_loop_filter(bandwidth,Δt)
    ω0 = bandwidth *1.89
    F = 1
    L = Δt*ω0^2
    C = 1
    D = sqrt(2)*ω0
    init_loop_filter(F, L, C, D)
end

"""
$(SIGNATURES)

Initialize a 2nd order loop_filter
Takes the noise bandwidth and the loop update time Δt.
Calculates the appropriate matrices for an 2nd order loop_filter and hands them to init_loop_filter
Returns a 2nd order loop_filter fucntion.
"""
function init_3rd_order_loop_filter(bandwidth,Δt)
    ω0 = bandwidth*1.2
    F = [1 Δt; 0 1]
    L = [Δt*1.1*ω0^2; Δt*ω0^3]
    C = [1 0]
    D = [2.4*ω0]
    init_loop_filter(F, L, C, D)
end
