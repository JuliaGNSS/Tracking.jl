module SignalCoverageTest

# Cross-cutting coverage guard: every concrete signal type GNSSSignals defines
# must either have a complete Tracking.jl per-signal API (correlator default,
# sync-search buffer width, loop-bandwidth defaults, sync detector) and actually
# survive a pass through `track`, or be listed in `UNSUPPORTED` with a reason.
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
    get_modulation,
    get_secondary_code_length,
    get_signal_name
import Tracking
using Tracking:
    AbstractCorrelator,
    NumAnts,
    TrackState,
    TrackedSat,
    default_carrier_loop_filter_bandwidth,
    default_code_loop_filter_bandwidth,
    detect_bit_or_secondary_code_sync,
    get_code_block_buffer_type,
    get_code_phase,
    get_default_correlator,
    get_sat_state,
    track

# Signals GNSSSignals defines that Tracking.jl deliberately does not support
# yet, keyed by type name, with the reason. Each needs a design decision, not
# just a dispatch method — so they fail loudly here rather than silently
# tracking with a nonsensical default.
const UNSUPPORTED = Dict(
    # Not a tracking signal: E5a-QP is the quick-acquisition aid, a dataless
    # BPSK(5) quasi-pilot whose 330-chip code repeats every 64.5 µs so a
    # receiver can find E5a cheaply and hand the satellite over to E5a-I /
    # E5a-Q (OS SIS ICD v2.2 §2.3.1.4). It carries neither data nor a secondary
    # code, so there is nothing to sync on, and its 64.5 µs primary code period
    # is two orders of magnitude shorter than any signal here: one code block —
    # Tracking's integration and loop-bandwidth unit — would run the loop at
    # ~15.5 kHz with a 279 Hz carrier bandwidth. Tracking it is not a matter of
    # adding a dispatch method.
    :GalileoE5aQP => "acquisition aid, not a tracking signal: no data, no secondary code, 64.5 µs primary code period",
)

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
        @test Tracking.get_num_ants(correlator) == num_ants
    end

    # Sync-search buffer must be an unsigned integer wide enough to hold one
    # whole secondary-code period — the hard rotation sweep's search horizon.
    B = @inferred get_code_block_buffer_type(signal)
    @test B <: Unsigned
    @test sizeof(B) * 8 >= get_secondary_code_length(signal)

    # Loop-bandwidth defaults are derived from the primary code period, so the
    # question is not whether they are positive — `BL = 0.018 / T` always is —
    # but whether the period they come from is a sane one to run a loop on. The
    # bounds below are absolute: 100 Hz is already a very wide carrier loop
    # (reference receivers sit at 5-25 Hz), and a pre-sync update rate above
    # 2 kHz means the tracker would be filtering faster than any real loop
    # needs to. Galileo E5a-QP fails both (279 Hz at 15.5 kHz) — which is why
    # it is in `UNSUPPORTED` rather than here.
    carrier_bandwidth = @inferred default_carrier_loop_filter_bandwidth(signal)
    code_bandwidth = @inferred default_code_loop_filter_bandwidth(signal)
    @test 0.0Hz < carrier_bandwidth < 100.0Hz
    @test 0.0Hz < code_bandwidth < 100.0Hz
    # The pre-sync integration is one primary code period, so its reciprocal is
    # the fastest rate the loop can be asked to update at.
    @test get_code_frequency(signal) / get_code_length(signal) < 2000.0Hz

    # The sync detector must accept the signal's own buffer width. Signals
    # without a bespoke method fall through to the soft CFAR path inside
    # `_buffer_find_bit`, which is exercised by the end-to-end pass below.
    if hasmethod(detect_bit_or_secondary_code_sync, Tuple{typeof(signal),Int,B,Int})
        @test detect_bit_or_secondary_code_sync(signal, 6, zero(B), 0) isa
              Tracking.SyncResult
    end
end

# Four code periods is far short of any detector's horizon (the soft CFAR ones
# need 2 × the secondary-code period), so this deliberately stops before sync:
# it is a smoke test that the whole chain — code replica, downconvert/correlate,
# discriminators, bit buffer, C/N₀ — runs for the signal and leaves the code
# phase where it was seeded. Sync behaviour is pinned per signal in the
# `test/<constellation>_<signal>.jl` files.
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

    # Four whole primary code periods, at zero Doppler so the expected code
    # phase after the pass is exactly the one we seeded; loop dynamics are
    # `track.jl`'s job.
    samples_per_period = code_length * sampling_frequency / code_frequency
    @test isinteger(upreferred(samples_per_period))
    num_samples = 4 * round(Int, upreferred(samples_per_period))

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
end

end
