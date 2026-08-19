module NoiseEstimatorTest

using Test: @test, @testset, @inferred, @test_throws
using Random: Xoshiro, randn
using StaticArrays: SVector, SMatrix
using Unitful: Hz, s, ms, ustrip, uconvert
using GNSSSignals: GPSL1CA, GPSL1C_P, GPSL5I, GalileoE1B
import Tracking
using Tracking:
    AbstractNoiseEstimator,
    CorrelatorNoiseEstimator,
    NoiseObservation,
    NoiseRefCN0Estimator,
    NWPRCN0Estimator,
    NumAnts,
    TrackState,
    TrackedSat,
    append_noise_observation!,
    get_noise_density,
    noise_observation,
    noise_observation_from_correlator,
    noise_observation_from_samples,
    update_noise!

# One dump of `n` samples of CN(0, σ²) noise despread by a ±1 code, as the three
# shapes a producer can report it in. `σ²` is the per-sample variance, so the
# density every one of them must land on is `σ²/f_s`.
function _white_dump(σ², n, fs, seed)
    rng = Xoshiro(seed)
    accumulation = complex(0.0, 0.0)
    sample_power = 0.0
    for _ = 1:n
        x = sqrt(σ² / 2) * complex(randn(rng), randn(rng))
        code = rand(rng, (-1.0, 1.0))
        accumulation += x * code
        sample_power += abs2(x)
    end
    (; accumulation, sample_power, n, fs)
end

@testset "the three builders land on the same N₀ scale" begin
    fs = 4e6Hz
    σ² = 3.0
    n = 4000

    # `noise_observation_from_samples` is the low-variance one — `M = n` looks —
    # so it pins the scale essentially exactly on a single dump.
    d = _white_dump(σ², n, fs, 20260806)
    from_samples = noise_observation_from_samples(d.sample_power, n, fs)
    @test ustrip(Hz^-1, from_samples.noise_density) ≈ ustrip(Hz^-1, σ² / fs) rtol = 0.05
    @test from_samples.num_sub_integrations == n
    @test from_samples.duration ≈ uconvert(s, n / fs)

    # A single despread dump is one look, so it carries 100 % relative error —
    # average many before comparing, which is the whole reason `M` is the weight.
    estimator = CorrelatorNoiseEstimator(; window_duration = 10.0s)
    for seed = 1:400
        dump = _white_dump(σ², 400, fs, seed)
        append_noise_observation!(estimator, noise_observation(dump.accumulation, 400, fs))
    end
    @test ustrip(Hz^-1, get_noise_density(estimator)) ≈ ustrip(Hz^-1, σ² / fs) rtol = 0.15

    # Pre-summing `M` dumps incoherently is the same density with `M` times the
    # weight, and reduces to the single-dump builder at `M = 1`.
    one = noise_observation(d.accumulation, n, fs)
    presummed = noise_observation_from_correlator(abs2(d.accumulation), 1, n, fs)
    @test presummed.noise_density == one.noise_density
    @test presummed.num_sub_integrations == one.num_sub_integrations
    @test presummed.duration == one.duration

    # `code_amplitude` undoes a multi-level code's integer scale, so a replica
    # scaled by `a` reports the same density.
    scaled = noise_observation(3.0 * d.accumulation, n, fs; code_amplitude = 3)
    @test ustrip(Hz^-1, scaled.noise_density) ≈ ustrip(Hz^-1, one.noise_density)

    # The sampling frequency may be spelled in any unit; the window's element
    # type must not depend on that.
    @test typeof(noise_observation(d.accumulation, n, 4.0e3Hz * 1000)) ===
          typeof(noise_observation(d.accumulation, n, fs))
end

@testset "simultaneous looks report their own span, not M times it" begin
    # The software source's `M` looks are the reference correlator's taps, which
    # all integrate the *same* N samples. Their density divides by `M·N` — three
    # taps are three looks at one scalar — but their wall-clock span is `N/f_s`,
    # and the window is bounded in time. Defaulting `duration` here would shrink
    # the window threefold.
    fs = 4e6Hz
    obs = noise_observation_from_correlator(3.0, 3, 3 * 4000, fs; duration = 4000 / fs)
    @test obs.num_sub_integrations == 3
    @test obs.duration ≈ uconvert(s, 4000 / fs)
    @test ustrip(Hz^-1, obs.noise_density) ≈ ustrip(Hz^-1, 3.0 / (3 * 4000 * fs))
end

@testset "the window is bounded in time and weighted by M" begin
    fs = 4e6Hz
    estimator = CorrelatorNoiseEstimator(; window_duration = 10.0ms)
    # 1 ms observations: the window settles at ten of them and keeps sliding.
    for k = 1:50
        append_noise_observation!(
            estimator,
            noise_observation_from_samples(4000.0 * k, 4000, fs),
        )
    end
    @test Base.length(estimator) == 10
    # Every entry carries the same `M`, so the mean is the plain mean of the last
    # ten densities — `k = 41 … 50`, i.e. `k̄ = 45.5` times `1/f_s`.
    @test ustrip(Hz^-1, get_noise_density(estimator)) ≈ ustrip(Hz^-1, 45.5 / fs)

    # An observation longer than the whole window is still kept: the window is a
    # lower bound on what it spans, never a reason to report nothing.
    long = CorrelatorNoiseEstimator(; window_duration = 1.0ms)
    append_noise_observation!(long, noise_observation_from_samples(1.0e5, 100_000, fs))
    @test Base.length(long) == 1
    @test !isnothing(get_noise_density(long))

    # `M` is what makes observations from different producers combinable: one
    # 4000-sample power-monitor look (M = 4000) against one despread dump
    # (M = 1) must come out essentially at the power monitor's value.
    mixed = CorrelatorNoiseEstimator(; window_duration = 10.0s)
    append_noise_observation!(mixed, noise_observation_from_samples(4000.0, 4000, fs))
    append_noise_observation!(
        mixed,
        noise_observation_from_correlator(4000.0 * 1000, 1, 4000, fs),
    )
    weighted = ustrip(Hz^-1, get_noise_density(mixed))
    @test weighted ≈ ustrip(Hz^-1, (4000 * 1.0 + 1 * 1000.0) / 4001 / fs)
end

@testset "an empty window has no density" begin
    estimator = CorrelatorNoiseEstimator()
    @test isnothing(get_noise_density(estimator))
    @test Base.length(estimator) == 0
    append_noise_observation!(
        estimator,
        noise_observation_from_samples(4000.0, 4000, 4e6Hz),
    )
    @test !isnothing(get_noise_density(estimator))

    # The public reader carries the `nothing` in its return type, which is what
    # makes it comfortable to call. The fold must never see that union, so it
    # goes through the internal splitter instead — and *that* has to be
    # concretely inferred, in both states, or the density threaded down to
    # `_apply_correlator_output` would box.
    D = Tracking.noise_density_type(estimator)
    @test Base.return_types(get_noise_density, (typeof(estimator),)) == [Union{Nothing,D}]
    empty_one = CorrelatorNoiseEstimator()
    @test @inferred(Tracking._noise_density_and_ready(empty_one)) == (zero(D), false)
    density, ready = @inferred Tracking._noise_density_and_ready(estimator)
    @test ready
    @test density == get_noise_density(estimator)
end

@testset "a zero measured floor is not a floor to divide by" begin
    # `get_noise_density` keeps its documented contract — the window is not empty,
    # so it reports what it holds. The fold's splitter is what decides whether the
    # figure is usable, and a zero density is not: `|P|²/0` is `Inf`, or `NaN` for
    # the zero prompt the same all-zero buffer produces, and a `NaN dB-Hz` passes
    # every lock threshold. A front-end dropout or buffer underrun reaches this
    # with no misuse at all.
    dead = CorrelatorNoiseEstimator()
    append_noise_observation!(dead, noise_observation_from_samples(0.0, 4000, 4e6Hz))
    D = Tracking.noise_density_type(dead)
    @test get_noise_density(dead) == zero(D)
    @test @inferred(Tracking._noise_density_and_ready(dead)) == (zero(D), false)

    # A producer's own arithmetic can hand over a non-finite one through the public
    # append path; same answer, and still inferred.
    for bad in (NaN, Inf)
        e = CorrelatorNoiseEstimator()
        append_noise_observation!(e, NoiseObservation(bad / 1.0Hz, 1, 1.0e-3s, Int16(1)))
        @test @inferred(Tracking._noise_density_and_ready(e)) == (zero(D), false)
    end
end

@testset "an observation lands whatever number type the producer used" begin
    # The builders' shared core `convert`s to the canonical `{NoiseDensity,
    # typeof(1.0s)}` rather than only `uconvert`ing the units, because a window
    # matches on the *type*: a producer spelling their sampling frequency `4f6Hz`
    # built a `Float32` observation that matched no method for the `Float64` window,
    # fell through to the abstract no-op and was dropped SILENTLY — the documented
    # hardware fill path, leaving the window empty forever and C/N₀ at -Inf.
    f32 = noise_observation(complex(1.0f0, 0.0f0), 4000, 4.0f6Hz)
    @test f32 isa NoiseObservation{Tracking.NoiseDensity,typeof(1.0s)}
    e32 = CorrelatorNoiseEstimator()
    append_noise_observation!(e32, f32)
    @test Base.length(e32) == 1
    @test ustrip(Hz^-1, get_noise_density(e32)) ≈ 6.25e-11 rtol = 1e-6

    # An `Int`-spelled one, and a duration given in `ms`, land on the same type.
    mixed = noise_observation_from_correlator(1.0, 1, 4000, 4_000_000Hz; duration = 1ms)
    @test mixed isa NoiseObservation{Tracking.NoiseDensity,typeof(1.0s)}

    # And an observation assembled by hand — never through a builder — is retyped
    # on append instead of being dropped.
    hand = NoiseObservation(1.0f-10 / 1.0f0Hz, 1, 1.0f-3s, Int16(3))
    @test !(hand isa NoiseObservation{Tracking.NoiseDensity,typeof(1.0s)})
    ehand = CorrelatorNoiseEstimator()
    append_noise_observation!(ehand, hand)
    @test Base.length(ehand) == 1
    @test ustrip(Hz^-1, get_noise_density(ehand)) ≈ 1.0e-10 rtol = 1e-6
end

@testset "the estimator is mutated in place, never rebuilt" begin
    # This is what lets per-signal state live in an immutable `TrackState`: the
    # window is a `Vector` field written in place, and the struct that comes back
    # is the identical object.
    estimator = CorrelatorNoiseEstimator()
    obs = noise_observation_from_samples(4000.0, 4000, 4e6Hz)
    @test append_noise_observation!(estimator, obs) === estimator
end

@testset "appending is allocation-free in steady state" begin
    # The four-times `sizehint!` headroom is what makes the FIFO's
    # `push!`/`popfirst!` pair measure exactly zero: at 1× or 2× Julia
    # periodically shifts the front offset back and reallocates. Wrapped in a
    # function because `@allocated` at module scope picks up boxing from untyped
    # global lookups.
    function push_many(estimator, obs, n)
        for _ = 1:n
            append_noise_observation!(estimator, obs)
        end
        estimator
    end
    estimator = CorrelatorNoiseEstimator()
    obs = noise_observation_from_samples(4000.0, 4000, 4e6Hz)
    push_many(estimator, obs, 5_000)          # warm up + fill the window

    # The contract is that nothing scales with the number of pushes. Asserting a
    # bare `== 0` would instead assert something about the compiler: on Julia
    # 1.10 the measurement carries a fixed ~16 B per *call* to the harness, which
    # does not move between ten thousand pushes and a million. Two very different
    # counts separate a per-push allocation — which would grow a hundredfold —
    # from a per-call one, which does not move at all.
    few = @allocated push_many(estimator, obs, 10_000)
    many = @allocated push_many(estimator, obs, 1_000_000)
    @test many == few
    @test few <= 128

    read_many(estimator, n) = (d = get_noise_density(estimator); for _ = 1:n
        d = get_noise_density(estimator)
    end; d)
    read_many(estimator, 100)
    few_reads = @allocated read_many(estimator, 1_000)
    many_reads = @allocated read_many(estimator, 10_000)
    @test many_reads == few_reads
    @test few_reads <= 128
end

@testset "update_noise! defaults to a no-op for a source fed from outside" begin
    # The two fill paths live on disjoint call graphs, so a source that is filled
    # by `append_noise_observation!` simply never reaches `update_noise!` — but
    # the default has to be a no-op rather than a `MethodError` so a user-defined
    # `AbstractNoiseEstimator` need only implement what it uses.
    struct ExternalOnlyNoiseEstimator <: AbstractNoiseEstimator end
    estimator = ExternalOnlyNoiseEstimator()
    @test update_noise!(estimator, nothing, 1, 10, nothing) === estimator
    @test isnothing(get_noise_density(estimator))
end

@testset "append_noise_observation! addresses a signal on TrackState" begin
    gpsl1 = GPSL1CA()
    obs = noise_observation_from_samples(4000.0, 4000, 4e6Hz)

    # Provisioned automatically, because the sat's estimator reads a density.
    single = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NoiseRefCN0Estimator())],
    )
    @test keys(single.noise_estimators) == (:GPSL1CA,)
    @test append_noise_observation!(single, obs) === single
    @test !isnothing(get_noise_density(single.noise_estimators.GPSL1CA))
    @test append_noise_observation!(single, obs, :GPSL1CA) === single
    # The signal type and an instance of it address the same window — a producer
    # usually holds one of those rather than the bare symbol.
    @test append_noise_observation!(single, obs, GPSL1CA) === single
    @test append_noise_observation!(single, obs, gpsl1) === single

    # A signal with no consumer has no estimator, and saying so beats a bare
    # NamedTuple `KeyError`. Reached by configuring an estimator that reads no
    # density, since the library default does read one.
    nwpr_only = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator(gpsl1))],
    )
    @test nwpr_only.noise_estimators === NamedTuple()
    @test_throws ArgumentError append_noise_observation!(nwpr_only, obs, :GPSL1CA)
    @test_throws ArgumentError append_noise_observation!(nwpr_only, obs)

    # Multi-signal: the bare form refuses to guess, and the two windows stay
    # apart. Across bands here, but see below for the case that per-band keying
    # could not express at all.
    multi = TrackState(;
        signals = (l1 = (gpsl1,), l5 = (GPSL5I(),)),
        noise_estimators = (
            GPSL1CA = CorrelatorNoiseEstimator(),
            GPSL5I = CorrelatorNoiseEstimator(),
        ),
    )
    @test_throws ArgumentError append_noise_observation!(multi, obs)
    append_noise_observation!(multi, obs, :GPSL1CA)
    @test !isnothing(get_noise_density(multi.noise_estimators.GPSL1CA))
    @test isnothing(get_noise_density(multi.noise_estimators.GPSL5I))
    # ... including at different sampling rates, which is the point of storing a
    # density rather than a power.
    append_noise_observation!(
        multi,
        noise_observation_from_samples(20_000.0, 20_000, 20e6Hz),
        :GPSL5I,
    )
    @test get_noise_density(multi.noise_estimators.GPSL1CA) !=
          get_noise_density(multi.noise_estimators.GPSL5I)
end

@testset "two signals on one band keep separate windows" begin
    # The case the per-band key could not express: GPS L1 C/A and Galileo E1B
    # share the L1 band, one antenna and one front end, but not a post-correlation
    # noise floor — BPSK(1) peaks where BOC(1,1) nulls. A producer must be able to
    # report a different density for each, and each consumer must read its own.
    gpsl1 = GPSL1CA()
    e1b = GalileoE1B()
    ts = TrackState(; signals = (gps = (gpsl1,), galileo = (e1b,)))
    @test keys(ts.noise_estimators) == (:GPSL1CA, :GalileoE1B)

    append_noise_observation!(
        ts,
        noise_observation_from_samples(4000.0, 4000, 4e6Hz),
        gpsl1,
    )
    append_noise_observation!(ts, noise_observation_from_samples(8000.0, 4000, 4e6Hz), e1b)
    l1ca_density = get_noise_density(ts.noise_estimators.GPSL1CA)
    e1b_density = get_noise_density(ts.noise_estimators.GalileoE1B)
    @test !isnothing(l1ca_density)
    @test !isnothing(e1b_density)
    # Twice the accumulated power over the same samples ⇒ twice the density, and
    # the two do not leak into one another.
    @test e1b_density ≈ 2 * l1ca_density
end

@testset "the per-signal density tuple folds and pairs positionally" begin
    # Replacing the fold's single `(density, ready)` pair with one pair per signal
    # is the only thing per-signal keying costs on the hot path, so it has to stay
    # a compile-time-shaped tuple: a `Union` here would box the density on its way
    # down to `_apply_correlator_output` and undo the allocation contracts in
    # `test/track_in_place.jl`.
    ts = TrackState(;
        signals = (modern = (GPSL1C_P(), GPSL1CA()), galileo = (GalileoE1B(),)),
    )
    @test keys(ts.noise_estimators) == (:GPSL1C_P, :GPSL1CA, :GalileoE1B)
    for group in Tuple(ts.groups)
        slot_type = eltype(group.satellites)
        returned = only(
            Base.return_types(
                Tracking._signal_noise_densities,
                (typeof(ts.noise_estimators), Type{slot_type}),
            ),
        )
        @test isconcretetype(returned)
        @test @inferred(
            Tracking._signal_noise_densities(ts.noise_estimators, slot_type)
        ) isa returned
    end

    # And the pairing is positional against the slot type, so a satellite carrying
    # one requiring and one non-requiring signal gets a *different* entry for each
    # — `nothing` (statically no source, surface it at the first record) beside a
    # real slot. Reading one signal's floor for the other is what this prevents.
    gpsl1, e1b = GPSL1CA(), GalileoE1B()
    mixed = TrackState(
        gpsl1,
        TrackedSat(
            (gpsl1, e1b),
            1,
            0.0,
            0.0Hz;
            cn0_estimator = (NoiseRefCN0Estimator(), NWPRCN0Estimator(e1b)),
        ),
    )
    @test keys(mixed.noise_estimators) == (:GPSL1CA,)
    noise = Tracking._signal_noise_densities(
        mixed.noise_estimators,
        eltype(only(Tuple(mixed.groups)).satellites),
    )
    @test length(noise) == 2
    @test last(noise[1]) == false           # provisioned, window still empty
    @test isnothing(first(noise[2]))        # no source, and none is wanted
    @test last(noise[2]) == true
end

@testset "nothing in this design is a mutable struct" begin
    # The whole per-signal design is shaped by this constraint — the
    # length-managed FIFO, the PRN carried in the observation, averaging in the
    # window rather than in a scalar — and nothing else enforces it.
    src = joinpath(dirname(dirname(@__DIR__)), "src")
    mutable_lines = String[]
    for (root, _, files) in walkdir(src), file in files
        endswith(file, ".jl") || continue
        for line in eachline(joinpath(root, file))
            startswith(line, "mutable struct") && push!(mutable_lines, line)
        end
    end
    @test isempty(mutable_lines)
    for T in (
        CorrelatorNoiseEstimator,
        NoiseObservation,
        NoiseRefCN0Estimator,
        Tracking.NoiseUpdateContext,
    )
        @test !ismutabletype(T)
    end
end

@testset "the builders take per-antenna input" begin
    # A producer with `M` elements reports all `M`: the same three shapes as the
    # scalar case, one dimension up. Collapsing them to a scalar first is what
    # the covariance exists to avoid — a beamformer's floor is `wᴴR̂w`, which a
    # single number cannot answer.
    fs = 4e6Hz
    n = 4000
    M = 3

    # A raw complex accumulation per antenna → the outer product `b·bᴴ`.
    b = SVector{M,ComplexF64}(2.0 + 0.0im, 0.0 + 1.0im, 1.0 - 1.0im)
    obs = @inferred noise_observation(b, n, fs)
    R = obs.noise_density
    @test R isa SMatrix{M,M}
    @test obs.num_sub_integrations == 1
    # Every element is the scalar builder's answer for that pair of antennas.
    for m = 1:M, k = 1:M
        @test R[m, k] ≈ (b[m] * conj(b[k])) / (n * fs)
    end
    # The diagonal is exactly what each antenna alone would have reported.
    for m = 1:M
        @test real(R[m, m]) ≈ noise_observation(b[m], n, fs).noise_density
    end

    # A pre-summed covariance from a correlator, and from a sample power monitor.
    pre_summed = b * b' + (2 * b) * (2 * b)'
    from_corr = @inferred noise_observation_from_correlator(pre_summed, 2, 2n, fs)
    @test from_corr.noise_density isa SMatrix{M,M}
    @test from_corr.noise_density ≈ pre_summed / (2n * fs)
    @test from_corr.num_sub_integrations == 2

    from_samples = @inferred noise_observation_from_samples(pre_summed, n, fs)
    @test from_samples.noise_density isa SMatrix{M,M}
    @test from_samples.noise_density ≈ pre_summed / (n * fs)

    # The window takes them, and averages them elementwise.
    estimator = CorrelatorNoiseEstimator(; num_ants = NumAnts(M))
    append_noise_observation!(estimator, obs)
    @test get_noise_density(estimator) ≈ R

    # A hand-built observation whose number type is not the canonical one is
    # still converted onto the window's, exactly as in the scalar case.
    hand_built = NoiseObservation(
        SMatrix{M,M,ComplexF32,M * M}(b * b') / 4.0f6Hz,
        1,
        1.0ms,
        Int16(3),
    )
    append_noise_observation!(estimator, hand_built)
    @test Base.length(estimator) == 2
end

end
