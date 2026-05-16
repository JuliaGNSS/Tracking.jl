module MultiBandTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, GPSL5I, GalileoE1B, L1, L5,
    gen_code, get_code_frequency, get_code_center_frequency_ratio

using Tracking:
    TrackState, Measurement, NumAnts, SignalGroup,
    track!, add_satellite!,
    band_key, band_keys,
    get_carrier_doppler, get_code_doppler, get_num_ants

# Build a synthetic complex sample buffer for one PRN at one band, with a
# given Doppler offset and starting code phase.
function _make_signal(
    ::Type{S}, prn, carrier_doppler, start_code_phase,
    num_samples, sampling_frequency,
) where {S}
    s = S()
    code_freq = carrier_doppler * get_code_center_frequency_ratio(s) +
        get_code_frequency(s)
    range = 0:(num_samples - 1)
    cis.(2π .* carrier_doppler .* range ./ sampling_frequency) .*
        gen_code(num_samples, s, prn, sampling_frequency, code_freq, start_code_phase)
end

@testset "band_keys for single-band TrackState" begin
    ts = TrackState(; signal = GPSL1CA())
    @test band_keys(ts) == (:l1,)
end

@testset "band_keys for multi-band TrackState" begin
    ts = TrackState(; signals = (
        legacy_gps_l1 = (GPSL1CA(),),
        gps_l5        = (GPSL5I(),),
    ))
    # Both bands appear in first-encounter order.
    @test band_keys(ts) == (:l1, :l5)
end

@testset "Two groups sharing a band yield one band key" begin
    ts = TrackState(; signals = (
        legacy_gps = (GPSL1CA(),),
        # GalileoE1B also lives on L1 — but using GPSL1CA twice is enough
        # to prove dedup, and avoids pulling in another signal.
        another_l1 = (GPSL1CA(),),
    ))
    @test band_keys(ts) == (:l1,)
end

@testset "Single-band track! still accepts bare buffer" begin
    sampling_frequency = 4e6Hz
    num_samples = 4000
    prn = 1
    # Integer-typed Hz must work end-to-end alongside Float64-typed Hz.
    carrier_doppler = 200Hz
    start_code_phase = 100.0

    track_state = TrackState(; signal = GPSL1CA())
    add_satellite!(track_state;
        prn, code_phase = start_code_phase,
        carrier_doppler = carrier_doppler - 20Hz,
    )

    signal = _make_signal(
        GPSL1CA, prn, carrier_doppler, start_code_phase,
        num_samples, sampling_frequency,
    )

    # Bare-buffer shortcut survives the multi-band refactor.
    track!(signal, track_state, sampling_frequency)
    # Doppler converged toward the injected value.
    @test abs(get_carrier_doppler(track_state, :default, prn) - carrier_doppler) < 50Hz
end

@testset "Multi-band track! with two synthetic bands" begin
    # GPS L1 C/A at 4 MHz, GPS L5I at 25 MHz, same observation duration (1 ms).
    fs_l1 = 4e6Hz
    fs_l5 = 25e6Hz
    n_l1 = 4000   # 4 MHz × 1 ms
    n_l5 = 25000  # 25 MHz × 1 ms

    prn_l1 = 1
    prn_l5 = 1
    # Mix int- and float-typed Hz to prove both forms work end-to-end.
    cd_l1 = 200Hz
    cd_l5 = -150.0Hz
    cp_l1 = 100.0
    cp_l5 = 250.0

    signal_l1 = _make_signal(GPSL1CA, prn_l1, cd_l1, cp_l1, n_l1, fs_l1)
    signal_l5 = _make_signal(GPSL5I, prn_l5, cd_l5, cp_l5, n_l5, fs_l5)

    track_state = TrackState(; signals = (
        legacy_gps_l1 = (GPSL1CA(),),
        gps_l5        = (GPSL5I(),),
    ))
    add_satellite!(track_state;
        prn = prn_l1, group = :legacy_gps_l1,
        code_phase = cp_l1, carrier_doppler = cd_l1 - 20Hz,
    )
    add_satellite!(track_state;
        prn = prn_l5, group = :gps_l5,
        code_phase = cp_l5, carrier_doppler = cd_l5 - 20.0Hz,
    )

    measurements = (
        l1 = Measurement(signal_l1, fs_l1),
        l5 = Measurement(signal_l5, fs_l5),
    )

    track!(measurements, track_state)

    # Each group's tracked Doppler should have moved toward its injected value.
    cd_l1_tracked = get_carrier_doppler(track_state, :legacy_gps_l1, prn_l1)
    cd_l5_tracked = get_carrier_doppler(track_state, :gps_l5, prn_l5)
    @test abs(cd_l1_tracked - cd_l1) < 50Hz
    @test abs(cd_l5_tracked - cd_l5) < 50.0Hz
end

@testset "Measurement keys mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA())
    # Wrong key: TrackState has band L1 (key :l1) but we pass :l5.
    m = Measurement(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((l5 = m,), ts)
end

@testset "Mismatched durations error (no tolerance)" begin
    ts = TrackState(; signals = (
        legacy_gps_l1 = (GPSL1CA(),),
        gps_l5        = (GPSL5I(),),
    ))
    add_satellite!(ts; prn = 1, group = :legacy_gps_l1, carrier_doppler = 0Hz)
    add_satellite!(ts; prn = 1, group = :gps_l5,        carrier_doppler = 0Hz)

    # L1 has 4000 samples @ 4 MHz = 1.000 ms; L5 has 25001 samples @ 25 MHz
    # ≈ 1.00004 ms. Off by one sample at L5's rate — must reject.
    m_l1 = Measurement(zeros(ComplexF64, 4000),  4e6Hz)
    m_l5 = Measurement(zeros(ComplexF64, 25001), 25e6Hz)
    @test_throws ArgumentError track!((l1 = m_l1, l5 = m_l5), ts)
end

@testset "Antenna shape mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA(), num_ants = NumAnts(2))
    add_satellite!(ts; prn = 1, carrier_doppler = 0Hz)
    # 2-antenna group expects a Matrix with 2 columns; pass a Vector.
    m = Measurement(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((l1 = m,), ts)
end

@testset "Per-band antenna counts via SignalGroup entries" begin
    # L1 group with 2 antennas, L5 group with 1. Constructed by passing
    # `SignalGroup` instances as `signals` entries.
    ts = TrackState(; signals = (
        legacy_gps_l1 = SignalGroup((GPSL1CA(),); num_ants = NumAnts(2)),
        gps_l5        = SignalGroup((GPSL5I(),);  num_ants = NumAnts(1)),
    ))
    add_satellite!(ts; prn = 1, group = :legacy_gps_l1, carrier_doppler = 0Hz)
    add_satellite!(ts; prn = 1, group = :gps_l5,        carrier_doppler = 0Hz)

    @test get_num_ants(ts, :legacy_gps_l1, 1) == 2
    @test get_num_ants(ts, :gps_l5, 1) == 1
end

@testset "Same band, mismatched NumAnts is rejected" begin
    # Two groups both on L1 (GPSL1CA + GalileoE1B share L1) but declaring
    # different antenna counts. Must error at TrackState construction.
    @test_throws ArgumentError TrackState(; signals = (
        gps     = SignalGroup((GPSL1CA(),);    num_ants = NumAnts(2)),
        galileo = SignalGroup((GalileoE1B(),); num_ants = NumAnts(1)),
    ))
end

@testset "Two groups on the same band with matching NumAnts is fine" begin
    ts = TrackState(; signals = (
        gps     = SignalGroup((GPSL1CA(),);    num_ants = NumAnts(2)),
        galileo = SignalGroup((GalileoE1B(),); num_ants = NumAnts(2)),
    ))
    @test band_keys(ts) == (:l1,)
end

end # module
