"""
$(SIGNATURES)

Generate a code replica for the given `signal`. The
replica contains `num_samples` prompt samples as well as an additional
number of early and late samples specified by `correlator_sample_shifts`.
The codefrequency is specified by `code_frequency`, while the sampling
rate is given by `sampling_frequency`. The phase of the first prompt
sample is given by `start_code_phase`. The generated signal is returned
in the array `code_replica` with the first generated sample written to
index start_sample.
"""
function gen_code_replica!(
    code_replica,
    signal::AbstractGNSSSignal,
    code_frequency,
    sampling_frequency,
    start_code_phase::AbstractFloat,
    start_sample::Integer,
    num_samples::Integer,
    correlator_sample_shifts::AbstractVector,
    prn::Integer,
)
    earliest_sample_shift = maximum(correlator_sample_shifts)
    latest_sample_shift = minimum(correlator_sample_shifts)
    total_samples = num_samples + earliest_sample_shift - latest_sample_shift
    gen_code!(
        view(code_replica, start_sample:(start_sample+total_samples-1)),
        signal,
        prn,
        sampling_frequency,
        code_frequency,
        start_code_phase,
        latest_sample_shift,
    )
    code_replica
end

# Best integer amplitude pair (a1, a2) GNSSSignals bakes into the CBOC Int8 embedded
# LUT: the continued-fraction convergent of the target amplitude ratio
# √(boc1_power) : √(1-boc1_power) whose sum stays within the Int8 table budget
# a1 + a2 ≤ 127. This mirrors GNSSSignals' (private) `_cboc_int_amplitudes`; the E1
# default (boc1_power = 10/11) yields (19, 6). Kept in lock-step with GNSSSignals by
# the `get_code_amplitude` regression test, which measures the RMS of a freshly
# generated E1B code and asserts it matches `sqrt(a1^2 + a2^2)`.
function _cboc_int_amplitudes(boc1_power::Real)
    ratio = sqrt(boc1_power / (1 - boc1_power))          # target a1 : a2
    h2, h1 = 0, 1                                         # convergent numerators   (h_{-2}, h_{-1})
    k2, k1 = 1, 0                                         # convergent denominators (k_{-2}, k_{-1})
    x = float(ratio)
    best = (1, 1)
    for _ = 1:64                                          # depth cap; the budget breaks far sooner
        a = floor(Int, x)
        h0 = a * h1 + h2
        k0 = a * k1 + k2
        h0 + k0 <= typemax(Int8) || break
        (h0 >= 1 && k0 >= 1) && (best = (h0, k0))
        h2, h1 = h1, h0
        k2, k1 = k1, k0
        frac = x - a
        frac > 1e-9 || break                             # exhausted the expansion
        x = 1 / frac
    end
    best
end

"""
$(SIGNATURES)

The RMS amplitude of the sampled code replica for `signal`'s modulation, used to
normalise the correlator so a unit-amplitude signal yields a unit prompt
regardless of modulation.

GNSSSignals v3's embedded Int8 LUT emits a ±1 replica for BPSK/BOC/TMBOC
(amplitude `1`), but a multi-level integer approximation for CBOC (Galileo E1B/E1C:
sub-carrier levels ±13/±25). That approximation bakes the sqrt-power split as an
integer pair `(a1, a2)` scaled to the Int8 budget, so its per-sample RMS is
`sqrt(a1^2 + a2^2)` (the ±1 sub-carrier cross term averages to zero over a code
period) — ≈ 19.92 for E1B, not `1`. Dividing the correlator by this in `normalize`
undoes the integer scale, restoring the unit-power convention of the old
sqrt-power float codes (so CBOC and BPSK stay at comparable power downstream).
"""
get_code_amplitude(signal::AbstractGNSSSignal) = get_code_amplitude(get_modulation(signal))
get_code_amplitude(::GNSSSignals.Modulation) = 1.0
function get_code_amplitude(modulation::GNSSSignals.CBOC)
    # `boc2_sign` (CBOC(+) for E1B, CBOC(−) for E1C) flips a2's sign but not a2², so
    # the RMS is identical for both — derive it from the unsigned amplitude pair.
    a1, a2 = _cboc_int_amplitudes(modulation.boc1_power)
    sqrt(Float64(a1)^2 + Float64(a2)^2)
end

"""
$(SIGNATURES)

Updates the code phase.
"""
function update_code_phase(
    signal::AbstractGNSSSignal,
    num_samples,
    code_frequency,
    sampling_frequency,
    start_code_phase,
    secondary_code_or_bit_found,
)
    secondary_code_or_bit_length =
        get_data_frequency(signal) == 0Hz ? get_secondary_code_length(signal) :
        Int(
            get_code_frequency(signal) /
            (get_data_frequency(signal) * get_code_length(signal)),
        )

    code_length =
        get_code_length(signal) *
        (secondary_code_or_bit_found ? secondary_code_or_bit_length : 1)
    mod(code_frequency * num_samples / sampling_frequency + start_code_phase, code_length)
    #    fixed_point = sizeof(Int) * 8 - 1 - min_bits_for_code_length(S)
    #    delta = floor(Int, code_frequency * 1 << fixed_point / sampling_frequency)
    #    fixed_point_start_phase = floor(Int, start_code_phase * 1 << fixed_point)
    #    phase_fixed_point = delta * num_samples + fixed_point_start_phase
    #    mod(phase_fixed_point / 1 << fixed_point, code_length)
end

"""
$(SIGNATURES)

Calculates the current code frequency.
"""
function get_current_code_frequency(signal::AbstractGNSSSignal, code_doppler)
    code_doppler + get_code_frequency(signal)
end
