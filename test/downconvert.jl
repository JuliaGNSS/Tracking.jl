@testset "Downconvert" begin
    downconverted_signal = StructArray{Complex{Float32}}(undef, 2500)

    phases = 2π * (1:2500) * 1000 / 2.5e6
    signal_struct =
        StructArray{Complex{Float32}}((Float32.(cos.(phases)), Float32.(sin.(phases))))

    carrier_replica = copy(signal_struct)

    Tracking.downconvert!(downconverted_signal, signal_struct, carrier_replica, 11, 2490)

    @test downconverted_signal[11:2500] ≈
          signal_struct[11:2500] .* conj(carrier_replica[11:2500])
    downconverted_signal = StructArray{Complex{Float32}}(undef, 2500)
    signal = Array(signal_struct)
    Tracking.downconvert!(downconverted_signal, signal, carrier_replica, 11, 2490)
    @test downconverted_signal[11:2500] ≈ signal[11:2500] .* conj(carrier_replica[11:2500])
end

@testset "Downconvert multiple signals" begin
    downconverted_signal = StructArray{Complex{Float32}}(undef, 2500, 3)

    phases = 2π * (1:2500) * 1000 / 2.5e6
    signal = StructArray{Complex{Float32}}((Float32.(cos.(phases)), Float32.(sin.(phases))))
    signal_mat_struct = repeat(signal; outer = (1, 3))

    carrier_replica = copy(signal)

    Tracking.downconvert!(
        downconverted_signal,
        signal_mat_struct,
        carrier_replica,
        11,
        2490,
    )

    @test downconverted_signal[11:2500, :] ≈
          signal_mat_struct[11:2500, :] .* conj(carrier_replica[11:2500])
    downconverted_signal = StructArray{Complex{Float32}}(undef, 2500, 3)
    signal_mat = Array(signal_mat_struct)
    Tracking.downconvert!(downconverted_signal, signal_mat, carrier_replica, 11, 2490)

    @test downconverted_signal[11:2500, :] ≈
          signal_mat[11:2500, :] .* conj(carrier_replica[11:2500])
end
