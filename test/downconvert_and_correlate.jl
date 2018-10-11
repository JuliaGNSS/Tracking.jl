@testset "Downconvert and correlate" begin

    num_ants = 4
    gpsl1 = GPSL1()
    carrier = gen_carrier.(1:4000, 50, 1.2, 4e6)
    code = gen_code.(Ref(gpsl1), 1:4000, 1023e3, 2, 4e6, 1)
    signal = [1, 1, 1, 1]' .* (carrier .* code)
    gen_carrier_replica(x) = gen_carrier(x, -50, -1.2, 4e6)
    gen_code_replica(x) = gen_code(gpsl1, x, 1023e3, 2, 4e6, 1)
    output = [zeros(ComplexF64, num_ants) for i = 1:3]
    Tracking.downconvert_and_correlate!(signal, output, 1, size(signal, 1), gen_carrier_replica, gen_code_replica, 2)
    @test output == [ones(num_ants) * 1952, ones(num_ants) * 4000, ones(num_ants) * 1952]
end
