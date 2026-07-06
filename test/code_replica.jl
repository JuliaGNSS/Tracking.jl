module CodeReplicaTest

using Test: @test, @testset, @inferred
using Statistics: mean
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, GPSL5I, GalileoE1B, GalileoE1B_BOC11, get_code, gen_code, get_code_frequency
using Tracking:
    gen_code_replica!, update_code_phase, get_current_code_frequency, get_code_amplitude

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

# The code-amplitude factor `normalize` divides out so a unit-amplitude signal gives
# a unit prompt regardless of modulation. It is `1` for every ±1 embedded-LUT code
# (BPSK, BOC(1,1)) and the RMS `sqrt(a1^2 + a2^2)` of the multi-level CBOC integer
# approximation (Galileo E1B: (a1, a2) = (19, 6) → sqrt(397) ≈ 19.92).
@testset "get_code_amplitude" begin
    # ±1 codes → exactly 1 (LOC and non-CBOC BOC dispatch to the generic method).
    @test @inferred(get_code_amplitude(GPSL1CA())) == 1.0        # LOC (BPSK)
    @test @inferred(get_code_amplitude(GPSL5I())) == 1.0         # LOC (BPSK)
    @test @inferred(get_code_amplitude(GalileoE1B_BOC11())) == 1.0  # BOC(1,1), not CBOC

    # CBOC (Galileo E1B): the analytic RMS of the (19, 6) Int8 approximation.
    e1b = GalileoE1B()
    @test @inferred(get_code_amplitude(e1b)) ≈ sqrt(19.0^2 + 6.0^2)

    # Cross-check against the ACTUAL embedded-LUT code: generate a full E1B code
    # period and measure its per-sample RMS. This pins `get_code_amplitude` to what
    # GNSSSignals really bakes — if GNSSSignals ever changes its CBOC Int8 amplitudes
    # this fails loudly rather than silently mis-normalising.
    fs = 15e6Hz
    fc = get_code_frequency(e1b)
    nsamp = round(Int, (fs / 1Hz) * 4e-3)        # one 4 ms (4092-chip) code period
    code = gen_code(nsamp, e1b, 1, fs, fc, 0.0)
    @test eltype(code) == Int8                    # multi-level integer replica
    # Widen to Int before squaring — abs2 on Int8 overflows (25^2 = 625 wraps).
    @test sqrt(mean(abs2, Int.(code))) ≈ get_code_amplitude(e1b) rtol = 1e-2
end

end
