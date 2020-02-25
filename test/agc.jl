@testset "AGC" begin
    @testset "Single signal" begin
        signal = rand(ComplexF64, 100)
        agc_signal = GainControlledSignal(signal, 5)
        real_max = maximum(real(x) for x in signal)
        imag_max = maximum(imag(x) for x in signal)

#        attenuation = sqrt(real_max^2 + imag_max^2)
        attenuation = max(real_max, imag_max)
        amplification = 1 << 5 / attenuation
        @test Tracking.get_attenuation(agc_signal) ≈ attenuation
        @test Tracking.get_amplitude_power(agc_signal) == 5
        @test real.(agc_signal.signal) ≈ floor.(Int16, real.(signal) .* amplification)
        @test imag.(agc_signal.signal) ≈ floor.(Int16, imag.(signal) .* amplification)

        signal = rand(Complex{Int16}, 100)
        agc_signal = GainControlledSignal(signal, 5)
        real_max = maximum(real(x) for x in signal)
        imag_max = maximum(imag(x) for x in signal)
#        attenuation = sqrt(float(real_max)^2 + float(imag_max)^2)
        attenuation = max(real_max, imag_max)
        amplification = 1 << 5 / attenuation
        @test Tracking.get_attenuation(agc_signal) ≈ attenuation
        @test Tracking.get_amplitude_power(agc_signal) == 5
        @test real.(agc_signal.signal) ≈ floor.(Int16, real.(signal) .* amplification)
        @test imag.(agc_signal.signal) ≈ floor.(Int16, imag.(signal) .* amplification)
    end

    @testset "Multiple signals" begin
        signal = randn(ComplexF64, 100, 4)
        agc_signal = GainControlledSignal(signal, 5)
        real_max = map(X -> maximum(real(x) for x in X), eachcol(signal))
        imag_max = map(X -> maximum(imag(x) for x in X), eachcol(signal))

#        attenuation = sqrt.(real_max.^2 .+ imag_max.^2)
        attenuation = max.(real_max, imag_max)
        amplification = (1 << 5) ./ attenuation
        @test Tracking.get_attenuation(agc_signal) ≈ attenuation
        @test Tracking.get_amplitude_power(agc_signal) == 5
        @test real.(agc_signal.signal) ≈ floor.(Int16, real.(signal) .* amplification')
        @test imag.(agc_signal.signal) ≈ floor.(Int16, imag.(signal) .* amplification')
    end
end
