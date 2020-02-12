@testset "Downconvert" begin

    downconverted_signal = StructArray{Complex{Int16}}(undef, 2500)

    phases = 2π * (1:2500) * 1000 / 2.5e6
    signal = StructArray{Complex{Int16}}((
        floor.(Int16, cos.(phases) * 1 << 7),
        floor.(Int16, sin.(phases) * 1 << 7)
    ))

    carrier_replica = copy(signal)

    Tracking.downconvert!(
        downconverted_signal,
        signal,
        carrier_replica,
        11,
        2490
    )

    @test signal[11:2500] .* conj(carrier_replica[11:2500]) == downconverted_signal[11:2500]
end

@testset "Downconvert multiple signals" begin

    downconverted_signal = StructArray{Complex{Int16}}(undef, 2500, 3)

    phases = 2π * (1:2500) * 1000 / 2.5e6
    signal = StructArray{Complex{Int16}}((
        floor.(Int16, cos.(phases) * 1 << 7),
        floor.(Int16, sin.(phases) * 1 << 7)
    ))
    signal_mat = repeat(signal, outer = (1,3))

    carrier_replica = copy(signal)

    Tracking.downconvert!(
        downconverted_signal,
        signal_mat,
        carrier_replica,
        11,
        2490
    )

    @test signal_mat[11:2500,:] .* conj(carrier_replica[11:2500]) ==
        downconverted_signal[11:2500,:]
end
