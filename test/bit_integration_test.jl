module BitDetectionIntegrationTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals: GPSL1, gen_code, get_code_frequency
using Tracking:
    SatState,
    TrackState,
    track,
    get_sat_state,
    get_code_phase,
    get_bits,
    get_num_bits,
    has_bit_or_secondary_code_been_found

@testset "Bit detection integration test" begin
    system = GPSL1()
    sampling_frequency = 5e6Hz
    num_samples = 5000
    code_frequency = get_code_frequency(system)
    carrier_doppler = 0.0Hz

    track_state =
        TrackState(system, [SatState(system, 1, sampling_frequency, 0, carrier_doppler)];)

    bits = vcat(ones(20), zeros(20), ones(1))
    foreach(enumerate(bits)) do (index, bit)
        code_phase = (index - 1) * num_samples * code_frequency / sampling_frequency
        carrier_phase =
            2π * (index - 1) * num_samples * carrier_doppler / sampling_frequency
        signal =
            (bit * 2 - 1) .* cis.(
                2π * (0:(num_samples-1)) * carrier_doppler / sampling_frequency .+
                carrier_phase,
            ) .*
            gen_code(num_samples, system, 1, sampling_frequency, code_frequency, code_phase)
        track_state = track(signal, track_state, sampling_frequency)
        @test has_bit_or_secondary_code_been_found(track_state) == (index >= 40)
        @test get_bits(track_state) == (index == 40 ? 2 : 0)
        @test get_num_bits(track_state) == (index == 40 ? 2 : 0)
    end
end

end
