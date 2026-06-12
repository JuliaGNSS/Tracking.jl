module MultiBandTest

using Test: @test, @testset, @inferred, @test_throws
using Unitful: Hz
using GNSSSignals:
    GPSL1CA, GPSL5I, GalileoE1B, L1, L5,
    gen_code, get_code_frequency, get_code_center_frequency_ratio

using Tracking:
    TrackState, BandMeasurement, NumAnts, SignalGroup,
    track, track!, add_satellite!,
    band_key, band_keys,
    get_samples, get_sampling_frequency, get_intermediate_frequency,
    get_carrier_doppler, get_code_doppler, get_num_ants,
    _validate_equal_durations

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

@testset "BandMeasurement kwarg constructor and accessors" begin
    buf = zeros(ComplexF64, 100)
    m_kw = BandMeasurement(buf; sampling_frequency = 4e6Hz, intermediate_frequency = 1.5e6Hz)
    @test get_samples(m_kw) === buf
    @test get_sampling_frequency(m_kw) == 4e6Hz
    @test get_intermediate_frequency(m_kw) == 1.5e6Hz

    # Kwarg form with default intermediate_frequency = zero(sampling_frequency).
    m_kw_default = BandMeasurement(buf; sampling_frequency = 4e6Hz)
    @test get_intermediate_frequency(m_kw_default) == 0.0Hz
end

@testset "BandMeasurement promotes mixed frequency types" begin
    buf = zeros(ComplexF64, 100)
    # Integer-typed IF alongside a Float64 sampling frequency — the most
    # natural spelling — must promote instead of throwing a MethodError.
    m = @inferred BandMeasurement(buf, 4e6Hz, 100Hz)
    @test get_sampling_frequency(m) === 4e6Hz
    @test get_intermediate_frequency(m) === 100.0Hz

    m_kw = BandMeasurement(buf; sampling_frequency = 4e6Hz, intermediate_frequency = 100Hz)
    @test get_intermediate_frequency(m_kw) === 100.0Hz
end

@testset "track with integer-typed intermediate_frequency" begin
    fs = 4e6Hz
    buf = _make_signal(GPSL1CA, 1, 200Hz, 100.0, 4000, fs)
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, code_phase = 100.0, carrier_doppler = 180Hz)
    # Integer IF must promote inside the bare-buffer wrapper, not throw.
    new_ts = track(buf, ts, fs; intermediate_frequency = 0Hz)
    @test new_ts isa TrackState
end

@testset "BandMeasurement rejects non-dense sample buffers (issue #134)" begin
    # The SIMD kernels read `samples` through raw pointers with dense
    # column-stride math; a non-contiguous strided view would silently
    # correlate the wrong data, so the constructor must reject it.
    buf = zeros(ComplexF64, 100, 4)
    vec_buf = zeros(ComplexF64, 100)

    # Dense arrays and contiguous views are accepted.
    @test BandMeasurement(vec_buf, 4e6Hz) isa BandMeasurement
    @test BandMeasurement(buf, 4e6Hz) isa BandMeasurement
    @test BandMeasurement(view(vec_buf, 1:50), 4e6Hz) isa BandMeasurement
    @test BandMeasurement(view(buf, :, 2), 4e6Hz) isa BandMeasurement
    @test BandMeasurement(view(buf, :, 2:3), 4e6Hz) isa BandMeasurement

    # Non-unit row stride / non-packed columns are rejected.
    @test_throws ArgumentError BandMeasurement(view(vec_buf, 1:2:99), 4e6Hz)
    @test_throws ArgumentError BandMeasurement(view(buf, 1:50, 2:3), 4e6Hz)
    @test_throws ArgumentError BandMeasurement(transpose(buf), 4e6Hz)
end

@testset "band_key for L1 and L5" begin
    @test @inferred(band_key(L1())) === :l1
    @test @inferred(band_key(L5())) === :l5
end

@testset "Single-BandMeasurement track/track! shortcut on single-band TrackState" begin
    fs = 4e6Hz
    num_samples = 4000
    prn = 1
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    buf = _make_signal(GPSL1CA, prn, carrier_doppler, start_code_phase, num_samples, fs)
    measurement = BandMeasurement(buf, fs)

    # Immutable variant.
    ts_imm = TrackState(; signal = GPSL1CA())
    ts_imm = add_satellite!(ts_imm; prn, code_phase = start_code_phase,
                   carrier_doppler = carrier_doppler - 20Hz)
    new_ts = track(measurement, ts_imm)
    @test abs(get_carrier_doppler(new_ts, :default, prn) - carrier_doppler) < 50Hz

    # In-place variant.
    ts_inp = TrackState(; signal = GPSL1CA())
    ts_inp = add_satellite!(ts_inp; prn, code_phase = start_code_phase,
                   carrier_doppler = carrier_doppler - 20Hz)
    track!(measurement, ts_inp)
    @test abs(get_carrier_doppler(ts_inp, :default, prn) - carrier_doppler) < 50Hz
end

@testset "Bare-buffer track! rejects multi-band TrackState" begin
    ts = TrackState(; signals = (
        gps_l1 = (GPSL1CA(),),
        gps_l5 = (GPSL5I(),),
    ))
    # `_single_band` must error — no single band to auto-key the buffer to.
    @test_throws ArgumentError track!(zeros(ComplexF64, 4000), ts, 4e6Hz)
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
    track_state = add_satellite!(track_state;
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
    track_state = add_satellite!(track_state;
        prn = prn_l1, group = :legacy_gps_l1,
        code_phase = cp_l1, carrier_doppler = cd_l1 - 20Hz,
    )
    track_state = add_satellite!(track_state;
        prn = prn_l5, group = :gps_l5,
        code_phase = cp_l5, carrier_doppler = cd_l5 - 20.0Hz,
    )

    measurements = (
        l1 = BandMeasurement(signal_l1, fs_l1),
        l5 = BandMeasurement(signal_l5, fs_l5),
    )

    track!(measurements, track_state)

    # Each group's tracked Doppler should have moved toward its injected value.
    cd_l1_tracked = get_carrier_doppler(track_state, :legacy_gps_l1, prn_l1)
    cd_l5_tracked = get_carrier_doppler(track_state, :gps_l5, prn_l5)
    @test abs(cd_l1_tracked - cd_l1) < 50Hz
    @test abs(cd_l5_tracked - cd_l5) < 50.0Hz
end

@testset "BandMeasurement keys mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA())
    # Wrong key: TrackState has band L1 (key :l1) but we pass :l5.
    m = BandMeasurement(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((l5 = m,), ts)
end

@testset "Mismatched durations error (no tolerance)" begin
    ts = TrackState(; signals = (
        legacy_gps_l1 = (GPSL1CA(),),
        gps_l5        = (GPSL5I(),),
    ))
    ts = add_satellite!(ts; prn = 1, group = :legacy_gps_l1, carrier_doppler = 0Hz)
    ts = add_satellite!(ts; prn = 1, group = :gps_l5,        carrier_doppler = 0Hz)

    # L1 has 4000 samples @ 4 MHz = 1.000 ms; L5 has 25001 samples @ 25 MHz
    # ≈ 1.00004 ms. Off by one sample at L5's rate — must reject.
    m_l1 = BandMeasurement(zeros(ComplexF64, 4000),  4e6Hz)
    m_l5 = BandMeasurement(zeros(ComplexF64, 25001), 25e6Hz)
    @test_throws ArgumentError track!((l1 = m_l1, l5 = m_l5), ts)
end

@testset "Equal durations across mixed Float32/Float64 sampling rates" begin
    # Both bands span exactly 1 ms, but the L5 front-end reports its rate
    # in Float32. Identical real durations must pass the check even though
    # `4000 / 4e6Hz` (Float64) and `25000 / 25f6Hz` (Float32) are not
    # `==` after division.
    m_l1 = BandMeasurement(zeros(ComplexF64, 4000),  4e6Hz)
    m_l5 = BandMeasurement(zeros(ComplexF64, 25000), 25f6Hz)
    @test _validate_equal_durations((l1 = m_l1, l5 = m_l5)) === nothing

    # A genuine one-sample mismatch must still throw, and the message
    # should print durations in seconds, not Hz^-1.
    m_bad = BandMeasurement(zeros(ComplexF64, 25001), 25f6Hz)
    err = try
        _validate_equal_durations((l1 = m_l1, l5 = m_bad))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test !occursin("Hz^-1", err.msg)
end

@testset "Antenna shape mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA(), num_ants = NumAnts(2))
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 0Hz)
    # 2-antenna group expects a Matrix with 2 columns; pass a Vector.
    m = BandMeasurement(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((l1 = m,), ts)
end

@testset "Per-band antenna counts via SignalGroup entries" begin
    # L1 group with 2 antennas, L5 group with 1. Constructed by passing
    # `SignalGroup` instances as `signals` entries.
    ts = TrackState(; signals = (
        legacy_gps_l1 = SignalGroup((GPSL1CA(),); num_ants = NumAnts(2)),
        gps_l5        = SignalGroup((GPSL5I(),);  num_ants = NumAnts(1)),
    ))
    ts = add_satellite!(ts; prn = 1, group = :legacy_gps_l1, carrier_doppler = 0Hz)
    ts = add_satellite!(ts; prn = 1, group = :gps_l5,        carrier_doppler = 0Hz)

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
