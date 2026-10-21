module NoiseEstimatorTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals: GPSL1CA, GPSL1C_P, GPSL5I, GalileoE1B
import Tracking
using Tracking: TrackState, TrackedSat
import TrackingLoops
using TrackingLoops:
    CorrelatorNoiseEstimator,
    NoiseObservation,
    NoiseRefCN0Estimator,
    NWPRCN0Estimator,
    append_noise_observation!,
    get_noise_density,
    noise_observation_from_samples

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
    # The signal type and an instance address the same window as the symbol.
    @test append_noise_observation!(single, obs, GPSL1CA) === single
    @test append_noise_observation!(single, obs, gpsl1) === single

    # A signal with no consumer has no estimator: `ArgumentError`, not a bare
    # `KeyError`.
    nwpr_only = TrackState(
        gpsl1,
        [TrackedSat(gpsl1, 1, 0.0, 0.0Hz; cn0_estimator = NWPRCN0Estimator(gpsl1))],
    )
    @test nwpr_only.noise_estimators === NamedTuple()
    @test_throws ArgumentError append_noise_observation!(nwpr_only, obs, :GPSL1CA)
    @test_throws ArgumentError append_noise_observation!(nwpr_only, obs)

    # Multi-signal: the bare form refuses to guess, and the windows stay apart.
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
    # ... including at different sampling rates (a density, not a power).
    append_noise_observation!(
        multi,
        noise_observation_from_samples(20_000.0, 20_000, 20e6Hz),
        :GPSL5I,
    )
    @test get_noise_density(multi.noise_estimators.GPSL1CA) !=
          get_noise_density(multi.noise_estimators.GPSL5I)
end

@testset "two signals on one band keep separate windows" begin
    # Same band, different post-correlation floors (see `AbstractNoiseEstimator`),
    # so a producer reports one density per signal.
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
    # Twice the power ⇒ twice the density, with no leak between the windows.
    @test e1b_density ≈ 2 * l1ca_density
end

@testset "the per-signal density tuple folds and pairs positionally" begin
    # `_signal_noise_densities` must return a concrete tuple: a `Union` would box
    # the density and break the allocation contracts in `test/track_in_place.jl`.
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

    # Pairing is positional against the slot type: a non-requiring signal gets
    # `nothing` beside the requiring one's real slot, never the other's floor.
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
    # The per-signal design depends on this and nothing else enforces it (see
    # `NoiseWindowTotals`).
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
        TrackingLoops.NoiseUpdateContext,
    )
        @test !ismutabletype(T)
    end
end

end
