module VectorTrackingTest

# Vector tracking end to end through `track!`: IF samples of GPS L1 C/A satellites
# broadcasting their (fixture) ephemerides as real LNAV bits, tracked with
# TrackingLoops' `VectorPLLAndDLL` and nothing vector-specific in between. The
# estimator syncs to the bits, decodes them, solves the scalar PVT and takes the
# satellites over by itself; this test only reads its solution.
#
# The true signal is synthesised from PositionVelocityTime's own models at a static
# receiver with a perfect clock: each satellite's transmit time is solved from the
# range at the receiver, and the code, carrier and bits follow from it. It takes a
# while — a first fix needs the whole clock and ephemeris, 18 s of bits after the
# start of a subframe 1 — so it only runs with
# `ENV["TRACKING_RUN_INTEGRATION_TEST"] = "true"`, like the Flexiband test.

using Test: @test, @testset, @test_skip
using Random: Xoshiro, randn
using Unitful: Hz, ustrip, @u_str
using LinearAlgebra: norm
using StaticArrays: SVector
using GNSSSignals:
    GPSL1CA,
    gen_code,
    get_code_frequency,
    get_center_frequency,
    get_code_center_frequency_ratio
using GNSSDecoder: GNSSDecoderState, get_time_of_week
using PositionVelocityTime:
    _precompile_states,
    _PRECOMPILE_GPS_L1CA_STATES,
    SatelliteState,
    SignalGroup,
    calc_pvt,
    calc_corrected_time,
    calc_satellite_position_and_velocity,
    calc_ρ_hat!,
    BiasColumns
using Tracking: TrackState, track!, add_satellite!
using TrackingLoops:
    VectorPLLAndDLL, navigation_solution, navigation_status, VTStatus, SPEED_OF_LIGHT

include("lnav_encoder.jl")

const FS = 2.048e6
const CHUNK = 8192 # 4 ms
const SIGNAL = GPSL1CA()
const CODE_FREQUENCY = ustrip(Hz, get_code_frequency(SIGNAL))
const CENTER_FREQUENCY = ustrip(Hz, get_center_frequency(SIGNAL))

mutable struct TrueSat
    decoder::GNSSDecoderState
    stream::LNAVStream
    clock_offset::Float64 # corrected − uncorrected transmit time
    transmit_time::Float64 # warm start of the range iteration
    carrier_phase::Float64 # cycles
end

# The corrected transmit time of a signal received at time of week `t` at `position`.
function transmit_time!(sat::TrueSat, position, t)
    columns = BiasColumns([1], 1, [0], 0)
    ξ = [position[1], position[2], position[3], 0.0]
    ρ = [0.0]
    t_t = sat.transmit_time
    for _ = 1:4
        orbit = calc_satellite_position_and_velocity(sat.decoder, t_t)
        calc_ρ_hat!(ρ, [SVector{3,Float64}(orbit.position)], ξ, columns)
        t_t = t - ρ[1] / SPEED_OF_LIGHT
    end
    sat.transmit_time = t_t
end

# The uncorrected transmit time (what the code and the bits follow) at receive time `t`.
uncorrected_transmit_time(sat::TrueSat, position, t) =
    transmit_time!(sat, position, t) - sat.clock_offset

function TrueSat(decoder, position, t)
    sat = TrueSat(decoder, LNAVStream(decoder.data), 0.0, t - 0.075, 0.0)
    t_t = transmit_time!(sat, position, t)
    # Invert the satellite clock correction around the transmit time.
    tow = get_time_of_week(decoder)
    u = t_t
    for _ = 1:3
        rate = 50.0
        num_bits = floor(Int, (u - tow) * rate)
        state = SatelliteState(;
            decoder = GNSSDecoderState(
                decoder;
                num_bits_after_valid_syncro_sequence = num_bits,
            ),
            system = SIGNAL,
            code_phase = (u - tow - num_bits / rate) * CODE_FREQUENCY,
            carrier_doppler = 0.0Hz,
        )
        u += t_t - calc_corrected_time(state)
    end
    sat.clock_offset = t_t - u
    sat
end

# Add one 4 ms chunk of satellite `sat`, received from time `t` on, to `samples`.
function add_signal!(samples, sat::TrueSat, prn, position, t)
    dt = CHUNK / FS
    u0 = uncorrected_transmit_time(sat, position, t)
    u1 = uncorrected_transmit_time(sat, position, t + dt)
    rate = (u1 - u0) / dt # transmit seconds per receive second
    doppler = (rate - 1) * CENTER_FREQUENCY
    code = gen_code(
        CHUNK,
        SIGNAL,
        prn,
        FS * Hz,
        CODE_FREQUENCY * rate * Hz,
        mod(u0 * CODE_FREQUENCY, 1023),
    )
    phase = sat.carrier_phase
    # At most one bit edge in 4 ms: the sample it falls on.
    k = floor(Int, u0 * 50)
    edge = floor(Int, ((k + 1) / 50 - u0) / rate * FS) + 1
    bits = (
        lnav_bit(sat.stream, k) ? 1.0f0 : -1.0f0,
        lnav_bit(sat.stream, k + 1) ? 1.0f0 : -1.0f0,
    )
    @inbounds for n = 1:CHUNK
        bit = n <= edge ? bits[1] : bits[2]
        samples[n] += ComplexF32(bit * code[n] * cis(2π * (phase + doppler * (n - 1) / FS)))
    end
    sat.carrier_phase = mod(phase + doppler * dt, 1.0)
    sat
end

@testset "Vector tracking from IF samples alone" begin
    if get(ENV, "TRACKING_RUN_INTEGRATION_TEST", "false") != "true"
        @info "Skipping the vector-tracking integration test. " *
              "Set ENV[\"TRACKING_RUN_INTEGRATION_TEST\"] = \"true\" to run it."
        @test_skip false
        return
    end
    states = _precompile_states(SIGNAL, _PRECOMPILE_GPS_L1CA_STATES, identity, SIGNAL)
    fix = calc_pvt(SignalGroup(SIGNAL, states); approximate_year = 2021)
    position = SVector(fix.position.x, fix.position.y, fix.position.z)
    t0 = next_subframe1_start(maximum(calc_corrected_time, states)) - 8.0
    decoders = [state.decoder for state in states[1:6]]
    sats = [TrueSat(decoder, position, t0) for decoder in decoders]

    estimator = VectorPLLAndDLL(
        SIGNAL;
        approximate_year = 2021,
        enable_ionospheric_correction = false,
        enable_tropospheric_correction = false,
    )
    track_state = TrackState(; signal = SIGNAL, doppler_estimator = estimator)
    # Acquired on the truth: the code phase and Doppler of the first sample.
    for (sat, decoder) in zip(sats, decoders)
        u0 = uncorrected_transmit_time(sat, position, t0)
        u1 = uncorrected_transmit_time(sat, position, t0 + 1e-3)
        doppler = ((u1 - u0) / 1e-3 - 1) * CENTER_FREQUENCY
        track_state = add_satellite!(
            track_state;
            prn = decoder.prn,
            code_phase = mod(u0 * CODE_FREQUENCY, 1023),
            carrier_doppler = doppler * Hz,
        )
    end

    # 45 dB-Hz: unit-amplitude signals in complex noise of N₀·fs.
    σ = Float32(sqrt(FS / (2 * 10^(45 / 10))))
    rng = Xoshiro(1)
    samples = Vector{ComplexF32}(undef, CHUNK)
    errors = Float64[]
    seeded_at = nothing
    num_chunks = round(Int, 30.0 / (CHUNK / FS))
    for k = 0:(num_chunks-1)
        t = t0 + k * CHUNK / FS
        for n in eachindex(samples)
            samples[n] = σ * complex(randn(rng, Float32), randn(rng, Float32))
        end
        for (sat, decoder) in zip(sats, decoders)
            add_signal!(samples, sat, decoder.prn, position, t)
        end
        track!(samples, track_state, FS * Hz)
        status = navigation_status(estimator)
        if status.running
            seeded_at = something(seeded_at, t - t0)
            pvt = navigation_solution(estimator)
            push!(
                errors,
                norm(SVector(pvt.position.x, pvt.position.y, pvt.position.z) - position),
            )
        end
    end
    status = navigation_status(estimator)
    @info "Vector tracking from IF samples" seeded_at max_tail_error =
        maximum(errors[(end-250):end]) position_std = status.position_std
    # The first fix comes once every satellite has decoded subframe 3, 8 + 18 s in.
    @test seeded_at !== nothing && 26.0 < seeded_at < 28.0
    @test status.running
    @test status.num_members == length(sats)
    @test length(navigation_solution(estimator).sats) == length(sats)
    @test maximum(errors[(end-250):end]) < 10.0
end

end
