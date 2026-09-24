module SignalCoverageTest

# Cross-cutting coverage guard: every concrete signal type GNSSSignals defines
# must either have a complete Tracking.jl per-signal API (correlator default,
# sync-search buffer width, loop-bandwidth defaults, integration policy, sync
# detector), actually survive a pass through `track`, and be documented in the
# capability matrix of `docs/src/signals.md` — or be listed in `UNSUPPORTED`
# with a reason.
#
# The signal list is discovered from the `AbstractGNSSSignal` type tree rather
# than hard-coded, so a GNSSSignals release that adds a signal fails here — a
# loud, one-place reminder — instead of only failing later at the first
# `get_default_correlator` MethodError in user code.

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
    TrackState,
    TrackedSat,
    TrackedSignal,
    get_code_phase,
    get_num_bits,
    get_preferred_num_code_blocks_to_integrate,
    get_sat_state,
    track
import TrackingLoops
using TrackingLoops:
    AbstractCorrelator,
    NumAnts,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    default_num_code_blocks_to_integrate,
    detect_bit_or_secondary_code_sync,
    get_code_block_buffer_type,
    get_default_correlator,
    max_num_code_blocks_to_integrate

# Signals GNSSSignals defines that Tracking.jl deliberately does not support
# yet, keyed by type name, with the reason. Each would need a design decision,
# not just a dispatch method — so they fail loudly here rather than silently
# tracking with a nonsensical default.
#
# Empty as of issue #236: every concrete signal GNSSSignals 4 defines is
# tracked. Galileo E5a-QP was the last entry — it is an acquisition aid, but
# that is a statement about how a receiver uses it, not about whether this
# package can track it, and it now integrates whole 31-block (2 ms) code cycles
# like any other pilot (`src/galileo/e5a_qp.jl`).
const UNSUPPORTED = Dict{Symbol,String}()

# Signals excluded from the end-to-end `track` pass only (their per-signal API is
# still checked), because one code block is too many samples for a unit test.
const NO_END_TO_END = Dict(
    # 767250 chips at 511.5 kcps = a 1.5 s primary code period; the tracker
    # integrates in whole code blocks, so the shortest possible pass is ~3.8 M
    # samples at the minimum sampling rate.
    :GPSL2CL => "1.5 s primary code period — one code block is millions of samples",
)

# Concrete leaves of the `AbstractGNSSSignal` tree that GNSSSignals itself
# defines, instantiated below via their zero-argument constructors. The module
# filter matters: other test files define their own fake signal types (to probe
# band/rate edge cases), and those are subtypes too but have no such
# constructor — and are not signals this package is expected to support.
function _signal_types()
    function leaves(T)
        subs = subtypes(T)
        isempty(subs) ? [T] : reduce(vcat, leaves.(subs))
    end
    types = filter(T -> Base.typename(T).module === GNSSSignals, leaves(AbstractGNSSSignal))
    sort(types; by = T -> string(nameof(T)))
end

# 5 samples per chip, raised to 25 MHz for the BOC-class signals so their
# sub-carrier is resolved too. Chosen so a whole primary code period is always an
# exact integer number of samples, so the expected code phase after a whole
# number of periods is exactly the one that was seeded.
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
        # Either Tracking supports it, or it is listed with a reason. A new
        # GNSSSignals signal fails here first.
        supported =
            hasmethod(get_default_correlator, Tuple{Base.typename(T).wrapper,NumAnts})
        @test _is_unsupported(T) || supported
        # And the two exception lists must stay honest: a signal that gains
        # support has to leave `UNSUPPORTED`, or the entry silently exempts it
        # from every check below.
        @test !(_is_unsupported(T) && supported)
    end
    # No stale entries either — a name GNSSSignals no longer defines (renamed,
    # removed) is an exemption for a signal that cannot fail anything.
    for name in keys(UNSUPPORTED)
        @test name in known
    end
    for name in keys(NO_END_TO_END)
        @test name in known
    end
end

@testset "Per-signal API — $(get_signal_name(signal))" for signal in SUPPORTED
    # Correlator default, for one and for several antennas.
    for num_ants in (1, 3)
        correlator = @inferred get_default_correlator(signal, NumAnts(num_ants))
        @test correlator isa AbstractCorrelator
        @test TrackingLoops.get_num_ants(correlator) == num_ants
    end

    # Sync-search buffer must be an unsigned integer wide enough to hold one
    # whole secondary-code period — the hard rotation sweep's search horizon.
    B = @inferred get_code_block_buffer_type(signal)
    @test B <: Unsigned
    @test sizeof(B) * 8 >= get_secondary_code_length(signal)

    # Integration policy: the length a fresh `TrackedSignal` starts at, and the
    # structural ceiling it may be raised to. The default has to be a divisor of
    # the ceiling, or `calc_num_code_blocks_to_integrate` would silently clamp
    # it down to one (issue #128's rule, applied to the default itself).
    default_blocks = @inferred default_num_code_blocks_to_integrate(signal)
    max_blocks = @inferred max_num_code_blocks_to_integrate(signal)
    @test 1 <= default_blocks <= max_blocks
    @test max_blocks % default_blocks == 0
    @test get_preferred_num_code_blocks_to_integrate(TrackedSignal(signal)) ==
          default_blocks

    # Loop-bandwidth defaults are derived from the primary code period, so the
    # question is not whether they are positive — `BL = 0.018 / T` always is —
    # but whether the loop that results is a sane one to run. That is a question
    # about the integration the signal actually *starts* at, not about one
    # primary block: the estimator scales the carrier bandwidth by `1/N` for an
    # N-block integration, so both bounds below are taken across `N`. The bounds
    # are absolute: 100 Hz is already a very wide carrier loop (reference
    # receivers sit at 5-25 Hz), and an update rate above 2 kHz means the
    # tracker would be filtering faster than any real loop needs to. Galileo
    # E5a-QP is what makes the distinction load-bearing — 279 Hz per 64.5 µs
    # block, but 9 Hz across the 31-block (2 ms) cycle it integrates.
    carrier_bandwidth = @inferred default_carrier_loop_filter_bandwidth(signal)
    code_bandwidth = @inferred default_code_loop_filter_bandwidth(signal)
    @test 0.0Hz < carrier_bandwidth / default_blocks < 100.0Hz
    @test 0.0Hz < code_bandwidth < 100.0Hz
    @test get_code_frequency(signal) / (get_code_length(signal) * default_blocks) < 2000.0Hz

    # The sync detector must accept the signal's own buffer width. Signals
    # without a bespoke method fall through to the soft CFAR path inside
    # `_buffer_find_bit`, which is exercised by the end-to-end pass below.
    if hasmethod(detect_bit_or_secondary_code_sync, Tuple{typeof(signal),Int,B,Int})
        @test detect_bit_or_secondary_code_sync(signal, 6, zero(B), 0) isa
              TrackingLoops.SyncResult
    end
end

# Four default integrations is far short of any detector's horizon (the soft
# CFAR ones need 2 × the secondary-code period), so this deliberately stops
# before sync for the signals that have something to sync on: it is a smoke test
# that the whole chain — code replica, downconvert/correlate, discriminators,
# bit buffer, C/N₀ — runs for the signal and leaves the code phase where it was
# seeded. Sync behaviour is pinned per signal in the
# `test/<constellation>_<signal>.jl` files.
#
# The length is counted in whole *default integrations* rather than code blocks
# so Galileo E5a-QP, whose default is a 31-block cycle, gets a pass long enough
# to complete one; for every other signal the default is one block and this is
# the same four code periods as before.
@testset "Runs a clean signal through track — $(get_signal_name(signal))" for signal in
                                                                              END_TO_END
    # PRN 6 is defined for every constellation here: GPS (1-32/63), Galileo
    # (1-50), and BeiDou — where it is also the first PRN BeiDou B2b-I defines
    # and a MEO/IGSO satellite, so B1I/B3I carry their NH20 overlay.
    prn = 6
    start_code_phase = 100
    sampling_frequency = _sampling_frequency(signal)
    code_length = get_code_length(signal)
    code_frequency = get_code_frequency(signal)

    # Four whole default integrations, at zero Doppler so the expected code
    # phase after the pass is exactly the one we seeded; loop dynamics are
    # `track.jl`'s job.
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
    # CBOC / TMBOC replicas carry integer sub-carrier amplitudes (±13/±25) while
    # BPSK ones are ±1; normalise so every signal enters `track` at unit peak.
    samples = ComplexF32.(code ./ maximum(abs, code))

    track_state = TrackState(signal, [TrackedSat(signal, prn, start_code_phase, 0.0Hz)])
    track_state = track(samples, track_state, sampling_frequency)
    sat_state = get_sat_state(track_state, prn)

    # The code phase may have been snapped forward by whole primary-code periods
    # when the secondary code locked (that is what seeds the position in the
    # tiered-code cycle), so compare modulo the primary code length.
    @test mod(get_code_phase(sat_state), code_length) ≈ start_code_phase atol = 0.1

    # Correlation and tracking are not navigation decoding: a pilot or
    # acquisition aid must come through the same chain without ever being asked
    # for a bit. `_calc_num_code_blocks_that_form_a_bit` is 0 for those signals
    # and the post-sync `buffer` path returns early on it, so `soft_bits` stays
    # empty no matter how long they track.
    if iszero(get_data_frequency(signal))
        @test TrackingLoops._calc_num_code_blocks_that_form_a_bit(signal) == 0
        @test get_num_bits(sat_state) == 0
    else
        # Data-bearing: one bit must be a whole number of primary code blocks,
        # or every integration length would straddle a bit boundary.
        blocks_per_bit = TrackingLoops._calc_num_code_blocks_that_form_a_bit(signal)
        @test blocks_per_bit >= 1
        @test upreferred(
            blocks_per_bit * code_length / code_frequency * get_data_frequency(signal),
        ) ≈ 1
    end
end

# The capability matrix in `docs/src/signals.md` is the user-facing answer to
# "what does Tracking.jl do with this signal", so it has to stay in step with
# the code. Two things are checked: that every signal GNSSSignals defines has a
# row (and no row names a signal it no longer defines), and that the one
# machine-readable cell — the integration policy, the number most likely to
# drift — matches the traits. The prose cells are left to review.
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
    # A newly added GNSSSignals signal has no row; a removed one leaves a stale
    # row behind. Both fail here.
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
