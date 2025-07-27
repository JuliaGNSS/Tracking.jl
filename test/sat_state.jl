@testset "Satellite state" begin
    gpsl1 = GPSL1()
    sat_state = @inferred SatState(gpsl1, 1, 5e6Hz, 10.0, 500.0Hz)
    @test get_prn(sat_state) == 1
    @test get_code_phase(sat_state) == 10.0
    @test get_code_doppler(sat_state) == 500.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 500.0Hz
    @test get_integrated_samples(sat_state) == 0
    @test get_signal_start_sample(sat_state) == 1
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test Tracking.found(get_secondary_code_or_bit_detector(sat_state)) == false
    @test length(get_bit_buffer(sat_state)) == 0

    acq = AcquisitionResults(
        gpsl1,
        5,
        5e6Hz,
        100.0Hz,
        524.6,
        45.0,
        1.0,
        randn(100, 100),
        -500:100.0:500,
    )
    sat_state = @inferred SatState(acq)
    @test get_prn(sat_state) == 5
    @test get_code_phase(sat_state) == 524.6
    @test get_code_doppler(sat_state) == 100.0Hz * get_code_center_frequency_ratio(gpsl1)
    @test get_carrier_phase(sat_state) == 0.0
    @test get_carrier_doppler(sat_state) == 100.0Hz
    @test get_integrated_samples(sat_state) == 0.0
    @test get_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_correlator(sat_state).accumulators == zeros(3)
    @test get_last_fully_integrated_filtered_prompt(sat_state) == complex(0.0, 0.0)
    @test Tracking.found(get_secondary_code_or_bit_detector(sat_state)) == false
    @test length(get_bit_buffer(sat_state)) == 0
end
