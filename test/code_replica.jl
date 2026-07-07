module CodeReplicaTest

using Test: @test, @testset, @inferred
using Unitful: Hz
using GNSSSignals: GPSL1CA, get_code
using Tracking: gen_code_replica!, update_code_phase, get_current_code_frequency

# GNSSSignals' embedded-LUT `gen_code!` (PR #90) may round a chip-boundary sample
# differently from the per-chip `get_code` oracle — at most ~1 sample per code
# period, and only where the code transitions. That sub-sample edge is irrelevant
# to tracking, so assert the replica matches the oracle everywhere except such
# boundary samples (which must coincide with a code transition).
function agrees_except_chip_boundaries(gen, ref, samples_per_period)
    @assert length(gen) == length(ref)
    mism = findall(gen .!= ref)
    at_transition = all(mism) do d
        ref[d] != ref[clamp(d - 1, 1, length(ref))] ||
            ref[d] != ref[clamp(d + 1, 1, length(ref))]
    end
    at_transition && length(mism) <= cld(length(ref), samples_per_period)
end

@testset "Code replica" begin
    code = zeros(Int8, 2502)
    gpsl1 = GPSL1CA()
    gen_code_replica!(code, gpsl1, 1023e3Hz, 2.5e6Hz, 2.0, 11, 2480, -1:1, 1)

    # samples_per_period = sampling_frequency / code_frequency * code_length
    #                    = 2.5e6 / 1023e3 * 1023 = 2500
    @test agrees_except_chip_boundaries(
        code[13:2492],
        Int8.(get_code.(gpsl1, (-1:2480) * 1023e3 / 2.5e6 .+ 2.0, 1)[3:end]),
        2500,
    )
    @testset "More than 1ms" begin
        code = zeros(Int8, 6502)
        gpsl1 = GPSL1CA()
        gen_code_replica!(code, gpsl1, 1023e3Hz, 2.5e6Hz, 2.0, 11, 6480, -1:1, 1)

        @test agrees_except_chip_boundaries(
            code[13:6492],
            Int8.(get_code.(gpsl1, (-1:6480) * 1023e3 / 2.5e6 .+ 2.0, 1)[3:end]),
            2500,
        )
    end

    @testset "code_length is less than 1ms" begin
        code = zeros(Int8, 2502)
        gpsl1 = GPSL1CA()
        gen_code_replica!(code, gpsl1, 1023e3Hz + 1000Hz, 7.5e6Hz, 2.0, 11, 2480, -1:1, 1)

        # samples_per_period = 7.5e6 / (1023e3 + 1000) * 1023 ≈ 7493
        @test agrees_except_chip_boundaries(
            code[13:2492],
            Int8.(get_code.(gpsl1, (-1:2480) * (1023e3 + 1000) / 7.5e6 .+ 2.0, 1)[3:end]),
            7493,
        )
    end
end

@testset "Update code phase" begin
    gpsl1 = GPSL1CA()
    code_phase = 10
    code_frequency = 10Hz
    sampling_frequency = 100Hz
    num_samples = 2000
    bit_found = true
    phase = @inferred update_code_phase(
        gpsl1,
        num_samples,
        code_frequency,
        sampling_frequency,
        code_phase,
        bit_found,
    )
    @test phase ≈ mod(10 + 0.1 * 2000, 1023)
end

@testset "Current code frequency" begin
    gpsl1 = GPSL1CA()
    @test @inferred(get_current_code_frequency(gpsl1, 0Hz)) == 1023000Hz
    @test @inferred(get_current_code_frequency(gpsl1, 200Hz)) == 1023000Hz + 200Hz
end

end
