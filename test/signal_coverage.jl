module SignalCoverageTest

# Coverage guard: every concrete GNSSSignals signal must have the full per-signal
# API, survive a `track` pass and have a row in `docs/src/signals.md` — or be
# listed in `UNSUPPORTED` with a reason. Signals are discovered from the type
# tree, so a GNSSSignals release that adds one fails here.

using Test: @test, @testset, @inferred
using Unitful: Hz, upreferred
using InteractiveUtils: subtypes
using GNSSSignals
using GNSSSignals:
    AbstractGNSSSignal,
    LOC,
    gen_code,
    get_code_center_frequency_ratio,
    get_code_frequency,
    get_code_length,
    get_data_frequency,
    get_modulation,
    get_secondary_code_length,
    get_signal_name
import Tracking
using Tracking:
    AbstractCorrelator,
    NumAnts,
    TrackState,
    TrackedSat,
    TrackedSignal,
    default_carrier_loop_filter_bandwidth,
    effective_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    default_num_code_blocks_to_integrate,
    detect_bit_or_secondary_code_sync,
    get_code_block_buffer_type,
    get_code_phase,
    get_default_correlator,
    get_num_bits,
    get_preferred_num_code_blocks_to_integrate,
    get_sat_state,
    max_num_code_blocks_to_integrate,
    track

# Signals deliberately not supported yet (type name => reason). Empty since
# issue #236 added Galileo E5a-QP.
const UNSUPPORTED = Dict{Symbol,String}()

# Signals excluded from the end-to-end `track` pass only (their per-signal API is
# still checked), because one code block is too many samples for a unit test.
const NO_END_TO_END =
    Dict(:GPSL2CL => "1.5 s primary code period — one code block is millions of samples")

# Concrete `AbstractGNSSSignal` leaves defined by GNSSSignals itself; the module
# filter excludes fake signal types that other test files define.
function _signal_types()
    function leaves(T)
        subs = subtypes(T)
        isempty(subs) ? [T] : reduce(vcat, leaves.(subs))
    end
    types = filter(T -> Base.typename(T).module === GNSSSignals, leaves(AbstractGNSSSignal))
    sort(types; by = T -> string(nameof(T)))
end

# 5 samples per chip (at least 25 MHz for BOC to resolve the sub-carrier), so a
# primary code period is an integer number of samples.
function _sampling_frequency(signal)
    fs = 5 * get_code_frequency(signal)
    get_modulation(signal) isa LOC ? fs : max(fs, 25.0e6Hz)
end

_is_unsupported(T) = haskey(UNSUPPORTED, Symbol(string(nameof(T))))
_skips_end_to_end(T) = haskey(NO_END_TO_END, Symbol(string(nameof(T))))

const SUPPORTED =
    [Base.typename(T).wrapper() for T in _signal_types() if !_is_unsupported(T)]
const END_TO_END = [
    Base.typename(T).wrapper() for
    T in _signal_types() if !_is_unsupported(T) && !_skips_end_to_end(T)
]

@testset "Every GNSSSignals signal type is accounted for" begin
    known = Set(Symbol(string(nameof(T))) for T in _signal_types())
    for T in _signal_types()
        # Supported, or listed with a reason.
        supported =
            hasmethod(get_default_correlator, Tuple{Base.typename(T).wrapper,NumAnts})
        @test _is_unsupported(T) || supported
        # A supported signal must leave `UNSUPPORTED`, or it skips every check.
        @test !(_is_unsupported(T) && supported)
    end
    # No stale entries.
    for name in keys(UNSUPPORTED)
        @test name in known
    end
    for name in keys(NO_END_TO_END)
        @test name in known
    end
end

@testset "Per-signal API — $(get_signal_name(signal))" for signal in SUPPORTED
    for num_ants in (1, 3)
        correlator = @inferred get_default_correlator(signal, NumAnts(num_ants))
        @test correlator isa AbstractCorrelator
        @test Tracking.get_num_ants(correlator) == num_ants
    end

    # Sync buffer must hold one whole secondary-code period.
    B = @inferred get_code_block_buffer_type(signal)
    @test B <: Unsigned
    @test sizeof(B) * 8 >= get_secondary_code_length(signal)

    # The default must divide the ceiling, or `calc_num_code_blocks_to_integrate`
    # clamps it to one (issue #128).
    default_blocks = @inferred default_num_code_blocks_to_integrate(signal)
    max_blocks = @inferred max_num_code_blocks_to_integrate(signal)
    @test 1 <= default_blocks <= max_blocks
    @test max_blocks % default_blocks == 0
    @test get_preferred_num_code_blocks_to_integrate(TrackedSignal(signal)) ==
          default_blocks

    # The loop must be sane at the integration the signal starts at (E5a-QP: a
    # 31-block cycle, not one 64.5 µs block): carrier below 100 Hz, update rate
    # below 2 kHz.
    integration_time =
        upreferred(get_code_length(signal) * default_blocks / get_code_frequency(signal))
    carrier_bandwidth = @inferred effective_carrier_loop_filter_bandwidth(
        default_carrier_loop_filter_bandwidth(signal),
        integration_time,
    )
    code_bandwidth = @inferred default_code_loop_filter_bandwidth(signal)
    @test 0.0Hz < carrier_bandwidth < 100.0Hz
    @test carrier_bandwidth * integration_time <= 0.09 + 1e-12
    @test 0.0Hz < code_bandwidth < 100.0Hz
    @test get_code_frequency(signal) / (get_code_length(signal) * default_blocks) < 2000.0Hz

    # The detector must accept the signal's buffer width. Signals without a
    # method use the soft CFAR path, exercised by the end-to-end pass below.
    if hasmethod(detect_bit_or_secondary_code_sync, Tuple{typeof(signal),Int,B,Int})
        @test detect_bit_or_secondary_code_sync(signal, 6, zero(B), 0) isa
              Tracking.SyncResult
    end
end

# Smoke test of the whole chain over four default integrations (so E5a-QP gets
# whole 31-block cycles). This mostly stops before sync; sync is pinned in the
# per-signal test files.
@testset "Runs a clean signal through track — $(get_signal_name(signal))" for signal in
                                                                              END_TO_END
    # PRN 6 exists everywhere (first B2b-I PRN) and is MEO/IGSO on BeiDou, so
    # B1I/B3I carry their NH20 overlay.
    prn = 6
    start_code_phase = 100
    sampling_frequency = _sampling_frequency(signal)
    code_length = get_code_length(signal)
    code_frequency = get_code_frequency(signal)

    # Zero Doppler, so the code phase should stay where it was seeded.
    samples_per_period = code_length * sampling_frequency / code_frequency
    @test isinteger(upreferred(samples_per_period))
    num_samples =
        4 *
        default_num_code_blocks_to_integrate(signal) *
        round(Int, upreferred(samples_per_period))

    code = gen_code(
        num_samples,
        signal,
        prn,
        sampling_frequency,
        code_frequency,
        start_code_phase,
    )
    # CBOC / TMBOC replicas have integer amplitudes (±13/±25); normalise to unit
    # peak.
    samples = ComplexF32.(code ./ maximum(abs, code))

    track_state = TrackState(signal, [TrackedSat(signal, prn, start_code_phase, 0.0Hz)])
    track_state = track(samples, track_state, sampling_frequency)
    sat_state = get_sat_state(track_state, prn)

    # Sync may snap the code phase by whole primary periods; compare modulo.
    @test mod(get_code_phase(sat_state), code_length) ≈ start_code_phase atol = 0.1

    # Dataless signals track without ever emitting a bit.
    if iszero(get_data_frequency(signal))
        @test Tracking._calc_num_code_blocks_that_form_a_bit(signal) == 0
        @test get_num_bits(sat_state) == 0
    else
        # One bit must span a whole number of primary code blocks.
        blocks_per_bit = Tracking._calc_num_code_blocks_that_form_a_bit(signal)
        @test blocks_per_bit >= 1
        @test upreferred(
            blocks_per_bit * code_length / code_frequency * get_data_frequency(signal),
        ) ≈ 1
    end
end

# Every signal has exactly one row in the `docs/src/signals.md` matrix, and its
# machine-readable integration-policy cell matches the traits. Prose is left to
# review.
@testset "docs/src/signals.md capability matrix is complete and current" begin
    path = joinpath(@__DIR__, "..", "docs", "src", "signals.md")
    @test isfile(path)
    rows = Dict{Symbol,Tuple{Int,Int}}()
    for line in eachline(path)
        m = match(r"^\|\s*`(\w+)`\s*\|.*\|\s*(\d+)\s*\(max\s*(\d+)\)\s*\|", line)
        isnothing(m) && continue
        rows[Symbol(m[1])] = (parse(Int, m[2]), parse(Int, m[3]))
    end
    documented = Set(keys(rows))
    known = Set(Symbol(string(nameof(T))) for T in _signal_types())
    @test setdiff(known, documented) == Set{Symbol}()
    @test setdiff(documented, known) == Set{Symbol}()
    for signal in SUPPORTED
        name = Symbol(string(nameof(typeof(signal))))
        haskey(rows, name) || continue
        documented_default, documented_max = rows[name]
        @test documented_default == default_num_code_blocks_to_integrate(signal)
        @test documented_max == max_num_code_blocks_to_integrate(signal)
    end
end

end
