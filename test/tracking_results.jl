@testset "Tracking results" begin
    results = Tracking.TrackingResults(
        TrackingState(GPSL1, 100Hz, 100),
        EarlyPromptLateCorrelator(NumAnts(2)),
        1.0Hz,
        1.0,
        true,
        Tracking.BitBuffer(),
        45dBHz
    )

    @test @inferred(get_carrier_doppler(results)) == 100Hz
    @test @inferred(get_carrier_phase(results)) == 0.0
    @test @inferred(get_code_doppler(results)) == 100Hz / 1540
    @test @inferred(get_code_phase(results)) == 100
    @test @inferred(get_correlator_carrier_phase(results)) ≈ 1.0 * 2π
    @test @inferred(get_correlator_carrier_frequency(results)) == 1.0Hz
    correlator = @inferred get_correlator(results)
    @test correlator == EarlyPromptLateCorrelator(NumAnts(2))
    @test @inferred(get_early(results)) == [0.0, 0.0]
    @test @inferred(get_prompt(results)) == [0.0, 0.0]
    @test @inferred(get_late(results)) == [0.0, 0.0]
    @test @inferred(get_bits(results)) == 0
    @test @inferred(get_num_bits(results)) == 0
    @test @inferred(get_secondary_code_or_bit_found(results)) == false
    @test @inferred(get_cn0(results)) == 45dBHz

    results = Tracking.TrackingResults(
        TrackingState(GPSL1, 100Hz, 100),
        EarlyPromptLateCorrelator(NumAnts(2)),
        1.0Hz,
        1.0,
        false,
        Tracking.BitBuffer(),
        45dBHz
    )
    correlator = @inferred get_correlator(results)
    @test get_early(correlator) == [0, 0]
    @test get_prompt(correlator) == [0, 0]
    @test get_late(correlator) == [0, 0]
end
