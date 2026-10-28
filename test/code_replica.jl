module CodeReplicaTest

using Test: @test, @testset
using Unitful: Hz
using GNSSSignals: GPSL1CA, get_code
using Tracking: gen_code_replica!

# `gen_code!` (GNSSSignals #90) may round a chip-boundary sample differently from
# `get_code`, at most ~1 sample per period, so such samples must sit on a code
# transition and every other one must match.
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

end
