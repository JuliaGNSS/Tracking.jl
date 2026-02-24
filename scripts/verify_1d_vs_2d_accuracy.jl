using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Statistics: mean, std
using Random: Random
using GNSSSignals:
    GPSL1, gen_code, get_code_frequency, get_code_center_frequency_ratio, get_code_length
using Unitful: Hz
using Tracking:
    CPUDownconvertAndCorrelator,
    SystemSatsState,
    SatState,
    TrackState,
    track,
    get_code_phase,
    get_code_doppler,
    get_carrier_doppler,
    use_fast_code_replica,
    EarlyPromptLateCorrelator

function run_tracking(;
    sampling_frequency,
    max_relative_code_error,
    initial_code_phase_offset = 0.0,
    cn0_dBHz = Inf,  # Inf = noiseless
    iterations = 2000,
    warmup = 500,
    seed = 1234,
)
    gpsl1 = GPSL1()
    code_frequency_base = get_code_frequency(gpsl1)
    carrier_doppler = 200Hz
    code_frequency =
        carrier_doppler * get_code_center_frequency_ratio(gpsl1) + code_frequency_base
    start_code_phase = 100.0
    start_carrier_phase = π / 2
    prn = 1

    num_samples_per_ms = Int(ceil(Float64(sampling_frequency / 1000Hz)))
    range_samples = 0:(num_samples_per_ms - 1)
    downconvert_and_correlator = CPUDownconvertAndCorrelator(Val(sampling_frequency))

    add_noise = isfinite(cn0_dBHz)
    if add_noise
        Random.seed!(seed)
        signal_amplitude = 10^(cn0_dBHz / 20)
        noise_std = sqrt(Float64(sampling_frequency / Hz))
    end

    track_state = TrackState(
        gpsl1,
        [SatState(gpsl1, prn, start_code_phase + initial_code_phase_offset, carrier_doppler)],
    )

    signal =
        cis.(
            2π .* carrier_doppler .* range_samples ./ sampling_frequency .+
            start_carrier_phase,
        ) .*
        gen_code(
            num_samples_per_ms,
            gpsl1,
            prn,
            sampling_frequency,
            code_frequency,
            start_code_phase,
        )
    if add_noise
        signal = signal .* signal_amplitude .+ randn(ComplexF64, num_samples_per_ms) .* noise_std
    end

    track_state = track(
        signal,
        track_state,
        sampling_frequency;
        downconvert_and_correlator,
        max_relative_code_error,
    )

    code_phase_errors = zeros(iterations)

    for i = 1:iterations
        carrier_phase =
            mod2pi(
                2π * carrier_doppler * num_samples_per_ms * i / sampling_frequency +
                start_carrier_phase + π,
            ) - π
        code_phase = mod(
            code_frequency * num_samples_per_ms * i / sampling_frequency +
            start_code_phase,
            get_code_length(gpsl1),
        )

        signal =
            cis.(
                2π .* carrier_doppler .* range_samples ./ sampling_frequency .+
                carrier_phase,
            ) .*
            gen_code(
                num_samples_per_ms,
                gpsl1,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase,
            )
        if add_noise
            signal = signal .* signal_amplitude .+ randn(ComplexF64, num_samples_per_ms) .* noise_std
        end

        track_state = track(
            signal,
            track_state,
            sampling_frequency;
            downconvert_and_correlator,
            max_relative_code_error,
        )

        true_code_phase = mod(
            code_frequency * num_samples_per_ms * (i + 1) / sampling_frequency +
            start_code_phase,
            get_code_length(gpsl1),
        )
        code_phase_errors[i] = get_code_phase(track_state) - true_code_phase
    end

    steady = code_phase_errors[(warmup + 1):end]
    return (;
        bias_chips = mean(steady),
        std_chips = std(steady),
        max_chips = maximum(abs.(steady)),
        errors = code_phase_errors,
    )
end

# Chip-to-meter conversion factor (C/A code: ~293 m per chip)
const CHIP_TO_METER = 293.05

gpsl1 = GPSL1()
code_frequency_base = get_code_frequency(gpsl1)

println("Verification: 1D (sample-quantized) vs 2D (exact chip shift) code replica accuracy")
println("Tracking 2000 iterations (1 ms each), steady-state analysis on last 1500\n")

for sampling_frequency in [2.5e6Hz, 3.5e6Hz, 4.0e6Hz, 5.0e6Hz]
    oversampling = Float64(sampling_frequency / code_frequency_base)
    corr = EarlyPromptLateCorrelator()
    is_fast = use_fast_code_replica(corr, sampling_frequency, code_frequency_base)

    println("=" ^ 80)
    println("fs = $(Float64(sampling_frequency / 1e6Hz)) MHz  ($(round(oversampling, digits=2))x oversampling)")
    println("Default threshold (0.2) selects: $(is_fast ? "1D (fast)" : "2D (accurate)")")
    println("-" ^ 80)

    r2d = run_tracking(; sampling_frequency, max_relative_code_error = 0.0)
    r1d = run_tracking(; sampling_frequency, max_relative_code_error = 1.0)

    println("  Path  |  Bias (m)    |  Std (m)     |  Max error (m)")
    println("  ------+--------------+--------------+----------------")
    for (label, r) in [("2D", r2d), ("1D", r1d)]
        bias_m = round(r.bias_chips * CHIP_TO_METER, sigdigits = 4)
        std_m = round(r.std_chips * CHIP_TO_METER, sigdigits = 4)
        max_m = round(r.max_chips * CHIP_TO_METER, sigdigits = 4)
        println("  $label    |  $(lpad(bias_m, 10)) |  $(lpad(std_m, 10)) |  $(lpad(max_m, 10))")
    end
    println()
end

println("""
Conclusion:
  The 1D (sample-quantized) path produces symmetric early/late correlations,
  giving the DLL discriminator an unbiased error signal. The 2D path uses exact
  fractional chip offsets but at low oversampling this creates asymmetric
  early/late correlations that bias the DLL. The 1D path is therefore more
  accurate for standard code tracking at all tested sampling rates.

  The 2D path may still be useful for specialized applications that need exact
  correlator spacing (e.g. multipath analysis), which is why max_relative_code_error
  is exposed as a parameter on the track() function.
""")

# Part 2: Convergence from initial code phase offsets
println("\n")
println("=" ^ 80)
println("Part 2: Convergence from initial code phase offsets")
println("=" ^ 80)
println()
println("Testing whether 1D and 2D paths can track the true code phase")
println("when the initial code phase estimate is offset from the true value.")
println("Using fs = 4.0 MHz, 10000 iterations (10 s), perfect Doppler init.\n")

sampling_frequency_convergence = 4.0e6Hz
iterations_convergence = 10000
# Print error at specific time snapshots
snapshots = [1, 10, 50, 100, 500, 1000, 2000, 5000, 10000]

for offset in [0.2, 0.4, 0.6, 0.8]
    println("-" ^ 80)
    println("Initial code phase offset: $offset chips ($(round(offset * CHIP_TO_METER, digits=1)) m)")
    println("-" ^ 80)
    println("  Time (ms) |  2D error (m)  |  1D error (m)")
    println("  ----------+----------------+---------------")

    r2d = run_tracking(;
        sampling_frequency = sampling_frequency_convergence,
        max_relative_code_error = 0.0,
        initial_code_phase_offset = offset,
        iterations = iterations_convergence,
        warmup = iterations_convergence - 1,  # don't care about summary stats here
    )
    r1d = run_tracking(;
        sampling_frequency = sampling_frequency_convergence,
        max_relative_code_error = 1.0,
        initial_code_phase_offset = offset,
        iterations = iterations_convergence,
        warmup = iterations_convergence - 1,
    )

    for t in snapshots
        e2d = round(r2d.errors[t] * CHIP_TO_METER, sigdigits = 4)
        e1d = round(r1d.errors[t] * CHIP_TO_METER, sigdigits = 4)
        println("  $(lpad(t, 9)) |  $(lpad(e2d, 12)) |  $(lpad(e1d, 12))")
    end
    println()
end

# Part 3: Effect of noise (CN0 = 45 dBHz)
println("\n")
println("=" ^ 80)
println("Part 3: Steady-state accuracy with and without noise")
println("=" ^ 80)
println()
println("Comparing 1D vs 2D with CN0 = 45 dBHz (typical open-sky) and noiseless.")
println("No initial code phase offset, perfect Doppler init, 2000 iterations.\n")

for sampling_frequency in [2.5e6Hz, 4.0e6Hz, 5.0e6Hz]
    oversampling = Float64(sampling_frequency / code_frequency_base)
    println("-" ^ 80)
    println("fs = $(Float64(sampling_frequency / 1e6Hz)) MHz  ($(round(oversampling, digits=2))x oversampling)")
    println("-" ^ 80)
    println("  Noise     |  Path  |  Bias (m)    |  Std (m)     |  Max error (m)")
    println("  ----------+--------+--------------+--------------+----------------")

    for (noise_label, cn0) in [("noiseless", Inf), ("45 dBHz", 45.0)]
        for (label, max_err) in [("2D", 0.0), ("1D", 1.0)]
            r = run_tracking(;
                sampling_frequency,
                max_relative_code_error = max_err,
                cn0_dBHz = cn0,
            )
            bias_m = round(r.bias_chips * CHIP_TO_METER, sigdigits = 4)
            std_m = round(r.std_chips * CHIP_TO_METER, sigdigits = 4)
            max_m = round(r.max_chips * CHIP_TO_METER, sigdigits = 4)
            println("  $(rpad(noise_label, 9)) |  $label    |  $(lpad(bias_m, 10)) |  $(lpad(std_m, 10)) |  $(lpad(max_m, 10))")
        end
    end
    println()
end

# Part 4: Convergence with noise
println("=" ^ 80)
println("Part 4: Convergence from code phase offsets with CN0 = 45 dBHz")
println("=" ^ 80)
println()
println("Using fs = 4.0 MHz, 10000 iterations (10 s), perfect Doppler init.\n")

for offset in [0.2, 0.4, 0.6, 0.8]
    println("-" ^ 80)
    println("Initial code phase offset: $offset chips ($(round(offset * CHIP_TO_METER, digits=1)) m)")
    println("-" ^ 80)
    println("  Time (ms) |  2D error (m)  |  1D error (m)")
    println("  ----------+----------------+---------------")

    r2d = run_tracking(;
        sampling_frequency = sampling_frequency_convergence,
        max_relative_code_error = 0.0,
        initial_code_phase_offset = offset,
        cn0_dBHz = 45.0,
        iterations = iterations_convergence,
        warmup = iterations_convergence - 1,
    )
    r1d = run_tracking(;
        sampling_frequency = sampling_frequency_convergence,
        max_relative_code_error = 1.0,
        initial_code_phase_offset = offset,
        cn0_dBHz = 45.0,
        iterations = iterations_convergence,
        warmup = iterations_convergence - 1,
    )

    for t in snapshots
        e2d = round(r2d.errors[t] * CHIP_TO_METER, sigdigits = 4)
        e1d = round(r1d.errors[t] * CHIP_TO_METER, sigdigits = 4)
        println("  $(lpad(t, 9)) |  $(lpad(e2d, 12)) |  $(lpad(e1d, 12))")
    end
    println()
end
