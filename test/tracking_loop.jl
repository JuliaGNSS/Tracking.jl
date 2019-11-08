@testset "Default post correlation Filter" begin
    post_corr_filter = Tracking.get_default_post_corr_filter(EarlyPromptLateCorrelator())
    @test @inferred(post_corr_filter(1.0 + 0.0im)) == 1.0 + 0.0im

    post_corr_filter = Tracking.get_default_post_corr_filter(
        EarlyPromptLateCorrelator(NumAnts(2))
    )
    @test @inferred(post_corr_filter(SVector(1.0 + 0.0im, 2.0 + 0.0im))) == 1.0 + 0.0im
end

@testset "Integration time" begin
    max_integration_time = 2ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(GPSL1, max_integration_time, bit_found)
    @test time == 1ms

    max_integration_time = 2ms
    bit_found = true
    time = @inferred Tracking.get_integration_time(GPSL1, max_integration_time, bit_found)
    @test time == 2ms

    max_integration_time = 2ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(
        GalileoE1B,
        max_integration_time,
        bit_found
    )
    @test time == 2ms

    max_integration_time = 10ms
    bit_found = false
    time = @inferred Tracking.get_integration_time(
        GalileoE1B,
        max_integration_time,
        bit_found
    )
    @test time == 4ms
end

@testset "Number of chips to integrate" begin
    max_integration_time = 1ms
    code_phase = 200
    bit_found = true
    chips = @inferred Tracking.get_num_chips_to_integrate(
        GPSL1,
        max_integration_time,
        code_phase,
        bit_found
    )
    @test chips == 1023 - 200

    code_phase = 1200
    chips = @inferred Tracking.get_num_chips_to_integrate(
        GPSL1,
        max_integration_time,
        code_phase,
        bit_found
    )
    @test chips == 2046 - 1200
end

@testset "Number of samples to integrate" begin
    max_integration_time = 1ms
    sample_frequency = 4e6Hz
    current_code_doppler = 0Hz
    current_code_phase = 0
    bit_found = true
    samples = @inferred Tracking.get_num_samples_left_to_integrate(
        GPSL1,
        max_integration_time,
        sample_frequency,
        current_code_doppler,
        current_code_phase,
        bit_found
    )
    @test samples == 4000
end

@testset "Carrier phase delta" begin
    intermediate_frequency = 100Hz
    carrier_doppler = 50Hz
    frequency = @inferred Tracking.get_current_carrier_frequency(
        intermediate_frequency,
        carrier_doppler
    )
    @test frequency == intermediate_frequency + carrier_doppler
end

@testset "Code phase delta" begin
    code_doppler = 1Hz
    frequency = @inferred Tracking.get_current_code_frequency(GPSL1, code_doppler)
    @test frequency == 1023e3Hz + 1Hz
end

@testset "Update phases" begin
    carrier_phase = π / 2
    carrier_frequency = 10Hz
    sample_frequency = 100Hz
    num_samples = 2000
    phase = @inferred Tracking.update_carrier_phase(
        num_samples,
        carrier_frequency,
        sample_frequency,
        carrier_phase
    )
    @test phase == mod(π / 2 + 0.1 * 2000 + 1, 2) - 1

    code_phase = 10
    code_frequency = 10Hz
    sample_frequency = 100Hz
    num_samples = 2000
    bit_found = true
    phase = @inferred Tracking.update_code_phase(
        GPSL1,
        num_samples,
        code_frequency,
        sample_frequency,
        code_phase,
        bit_found
    )
    @test phase == mod(10 + 0.1 * 2000, 1023)
end

@testset "Correlate" begin
    correlator = EarlyPromptLateCorrelator()
    early_late_sample_shift = 2
    start_sample = 1
    num_samples = 4000
    carrier_frequency = 800e3Hz
    code_frequency = 1023e3Hz
    sample_frequency = 4000e3Hz
    start_carrier_phase = 0.25
    @testset "Correlate for satellite $prn" for prn in 1:32
        range = 0:3999
        signal = cis.(2π .* carrier_frequency .* range ./ sample_frequency .+ π ./ 2) .*
            get_code.(GPSL1, code_frequency .* range ./ sample_frequency .+ 0, prn)
        start_code_phase = 0

        correlator_result = @inferred Tracking.correlate(
            GPSL1,
            correlator,
            signal,
            prn,
            early_late_sample_shift,
            start_sample,
            num_samples,
            carrier_frequency,
            code_frequency,
            sample_frequency,
            start_carrier_phase,
            start_code_phase
        )
        difference_early_late = abs(
            get_early(correlator_result) - get_late(correlator_result)
        )
        # WHY??
        @test get_early(correlator_result) ≈ 1952 || get_early(correlator_result) ≈ 2080.0 ||
            get_early(correlator_result) ≈ 1824.0
        @test get_prompt(correlator_result) ≈ 4000
        @test get_late(correlator_result) ≈ 1952 || get_late(correlator_result) ≈ 2080.0 ||
            get_late(correlator_result) ≈ 1824.0
        @test difference_early_late ≈ 0 atol = 4.5e-11
    end


    range = 0:7999
    num_samples = 8000
    prn = 1
    signal = cis.(2π .* carrier_frequency .* range ./ sample_frequency .+ π ./ 2) .*
        get_code.(GPSL1, code_frequency .* range ./ sample_frequency .+ 0, prn)
    start_code_phase = 0
    correlator_result = @inferred Tracking.correlate(
        GPSL1,
        correlator,
        signal,
        prn,
        early_late_sample_shift,
        start_sample,
        num_samples,
        carrier_frequency,
        code_frequency,
        sample_frequency,
        start_carrier_phase,
        start_code_phase
    )
    @test get_early(correlator_result) ≈ 3904
    @test get_prompt(correlator_result) ≈ 8000
    @test get_late(correlator_result) ≈ 3904

end

@testset "Doppler aiding" begin

    init_carrier_doppler = 10Hz
    init_code_doppler = 1Hz
    carrier_freq_update = 2Hz
    code_freq_update = -0.5Hz
    velocity_aiding = 3Hz

    carrier_freq, code_freq = @inferred Tracking.aid_dopplers(
        GPSL1,
        init_carrier_doppler,
        init_code_doppler,
        carrier_freq_update,
        code_freq_update,
        velocity_aiding
    )

    @test carrier_freq == 10Hz + 2Hz + 3Hz
    @test code_freq == 1Hz + (2Hz + 3Hz) / 1540 - 0.5Hz
end

@testset "Signal" begin
    correlator = EarlyPromptLateCorrelator()
    signal = collect(1:10)
    sample = 5
    signal_sample = @inferred Tracking.get_signal(correlator, signal, sample)
    @test signal_sample == 5
    @test @inferred(Tracking.get_num_samples(signal)) == 10

    correlator = EarlyPromptLateCorrelator(NumAnts(2))
    signal = [1, 2] * collect(1:10)'
    sample = 5
    signal_sample = @inferred Tracking.get_signal(correlator, signal, sample)
    @test signal_sample == [5, 10]
    @test typeof(signal_sample) <: SVector

    @test @inferred(Tracking.get_num_samples(signal)) == 10
end

@testset "Tracking" begin

    carrier_doppler = 200Hz
    start_code_phase = 100
    code_frequency = carrier_doppler / 1540 + 1023kHz
    sample_frequency = 4MHz
    prn = 1
    range = 0:3999
    start_carrier_phase = π / 2
    state = TrackingState(GPSL1, carrier_doppler, start_code_phase)


    signal = cis.(
            2π .* carrier_doppler .* range ./ sample_frequency .+ start_carrier_phase
        ) .*
        get_code.(
            GPSL1,
            code_frequency .* range ./ sample_frequency .+ start_code_phase,
            prn
        )

    track_result = @inferred track(signal, state, prn, sample_frequency)

    iterations = 1000
    code_phases = zeros(iterations)
    carrier_phases = zeros(iterations)
    tracked_code_phases = zeros(iterations)
    tracked_carrier_phases = zeros(iterations)
    tracked_code_dopplers = zeros(iterations)
    tracked_carrier_dopplers = zeros(iterations)
    for i = 1:iterations
        carrier_phase = mod2pi(2π * carrier_doppler * 4000 * i / sample_frequency +
            start_carrier_phase + π) - π
        code_phase = mod(
            code_frequency * 4000 * i / sample_frequency + start_code_phase,
            1023
        )
        signal = cis.(
                2π .* carrier_doppler .* range ./ sample_frequency .+ carrier_phase
            ) .*
            get_code.(
                GPSL1,
                code_frequency .* range ./ sample_frequency .+ code_phase,
                prn
            )
        track_result = @inferred track(
            signal,
            get_state(track_result),
            prn,
            sample_frequency
        )
        comp_carrier_phase = mod2pi(2π * carrier_doppler * 4000 * (i + 1) /
            sample_frequency + start_carrier_phase + π) - π
        comp_code_phase = mod(
            code_frequency * 4000 * (i + 1) / sample_frequency + start_code_phase,
            1023
        )
        tracked_code_phases[i] = get_code_phase(track_result)
        tracked_carrier_phases[i] = get_carrier_phase(track_result)
        tracked_carrier_dopplers[i] = get_carrier_doppler(track_result)/Hz
        tracked_code_dopplers[i] = get_code_doppler(track_result)/Hz
        code_phases[i] = comp_code_phase
        carrier_phases[i] = comp_carrier_phase
    end
    @test tracked_code_phases[end] ≈ code_phases[end] atol = 2e-4
    @test tracked_carrier_phases[end] ≈ carrier_phases[end] atol = 5e-2

#    using PyPlot
#    pygui(true)
#    figure("carrier_phases")
#    plot(tracked_carrier_phases - carrier_phases)
#    plot(carrier_phases)
#    figure("Code Phases")
#    plot(300 * (tracked_code_phases .- code_phases))
#    figure("Carrier Doppler")
#    plot(tracked_carrier_dopplers)
#    figure("Code Doppler")
#    plot(tracked_code_dopplers)

end
