module MultiBandTest

using Test: @test, @testset, @inferred, @test_throws
using GNSSSignals: CBOC
using Unitful: Hz
using GNSSSignals:
    GPSL1CA,
    GPSL5I,
    GalileoE1B,
    L1,
    L5,
    gen_code,
    get_band_id,
    get_chip_rate,
    get_chip_rate_to_center_frequency_ratio

using Tracking:
    TrackState,
    BandSamples,
    SignalGroup,
    track,
    track!,
    add_satellite!,
    get_band_ids,
    get_carrier_doppler,
    _validate_equal_durations
using TrackingLoops: NumAntennas, get_num_antennas

# Samples for one PRN at a Doppler and starting code phase.
function _make_signal(
    ::Type{S},
    prn,
    carrier_doppler,
    start_code_phase,
    num_samples,
    sampling_frequency,
) where {S}
    s = S()
    code_freq =
        carrier_doppler * get_chip_rate_to_center_frequency_ratio(s) + get_chip_rate(s)
    range = 0:(num_samples-1)
    cis.(2π .* carrier_doppler .* range ./ sampling_frequency) .*
    gen_code(num_samples, s, prn, sampling_frequency, code_freq, start_code_phase)
end

@testset "BandSamples takes the intermediate frequency as a keyword" begin
    buf = zeros(ComplexF64, 100)
    low_if = BandSamples(buf, 4e6Hz; intermediate_frequency = 1.5e6Hz)
    @test low_if.samples === buf
    @test low_if.sampling_frequency == 4e6Hz
    @test low_if.intermediate_frequency == 1.5e6Hz

    baseband = BandSamples(buf, 4e6Hz)
    @test baseband.intermediate_frequency == 0.0Hz
    # Positionally, an intermediate frequency would read as a second rate.
    @test_throws MethodError BandSamples(buf, 4e6Hz, 1.5e6Hz)
end

@testset "BandSamples promotes mixed frequency types" begin
    buf = zeros(ComplexF64, 100)
    # An integer IF with a Float64 sampling frequency promotes.
    band_samples = @inferred BandSamples(buf, 4e6Hz; intermediate_frequency = 100Hz)
    @test band_samples.sampling_frequency === 4e6Hz
    @test band_samples.intermediate_frequency === 100.0Hz
end

@testset "track with integer-typed intermediate_frequency" begin
    fs = 4e6Hz
    buf = _make_signal(GPSL1CA, 1, 200Hz, 100.0, 4000, fs)
    ts = TrackState(; signal = GPSL1CA())
    ts = add_satellite!(ts; prn = 1, code_phase = 100.0, carrier_doppler = 180Hz)
    new_ts = track(buf, ts, fs; intermediate_frequency = 0Hz)
    @test new_ts isa TrackState
end

@testset "BandSamples rejects non-dense sample buffers (issue #134)" begin
    # The kernels assume dense columns, so a strided view would silently
    # correlate the wrong data.
    buf = zeros(ComplexF64, 100, 4)
    vec_buf = zeros(ComplexF64, 100)

    @test BandSamples(vec_buf, 4e6Hz) isa BandSamples
    @test BandSamples(buf, 4e6Hz) isa BandSamples
    @test BandSamples(view(vec_buf, 1:50), 4e6Hz) isa BandSamples
    @test BandSamples(view(buf, :, 2), 4e6Hz) isa BandSamples
    @test BandSamples(view(buf, :, 2:3), 4e6Hz) isa BandSamples

    @test_throws ArgumentError BandSamples(view(vec_buf, 1:2:99), 4e6Hz)
    @test_throws ArgumentError BandSamples(view(buf, 1:50, 2:3), 4e6Hz)
    @test_throws ArgumentError BandSamples(transpose(buf), 4e6Hz)
end

@testset "get_band_id for L1 and L5" begin
    @test @inferred(get_band_id(L1())) === :L1
    @test @inferred(get_band_id(L5())) === :L5
end

@testset "Single-BandSamples track/track! shortcut on single-band TrackState" begin
    fs = 4e6Hz
    num_samples = 4000
    prn = 1
    carrier_doppler = 200Hz
    start_code_phase = 100.0
    buf = _make_signal(GPSL1CA, prn, carrier_doppler, start_code_phase, num_samples, fs)
    band_samples = BandSamples(buf, fs)

    ts_imm = TrackState(; signal = GPSL1CA())
    ts_imm = add_satellite!(
        ts_imm;
        prn,
        code_phase = start_code_phase,
        carrier_doppler = carrier_doppler - 20Hz,
    )
    new_ts = track(band_samples, ts_imm)
    @test abs(get_carrier_doppler(new_ts, :default, prn) - carrier_doppler) < 50Hz

    ts_inp = TrackState(; signal = GPSL1CA())
    ts_inp = add_satellite!(
        ts_inp;
        prn,
        code_phase = start_code_phase,
        carrier_doppler = carrier_doppler - 20Hz,
    )
    track!(band_samples, ts_inp)
    @test abs(get_carrier_doppler(ts_inp, :default, prn) - carrier_doppler) < 50Hz
end

@testset "Bare-buffer track! rejects multi-band TrackState" begin
    ts = TrackState(; signals = (gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
    # No single band to key a bare buffer to.
    @test_throws ArgumentError track!(zeros(ComplexF64, 4000), ts, 4e6Hz)
end

@testset "get_band_ids for single-band TrackState" begin
    ts = TrackState(; signal = GPSL1CA())
    @test get_band_ids(ts) == (:L1,)
end

@testset "get_band_ids for multi-band TrackState" begin
    ts = TrackState(; signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
    @test get_band_ids(ts) == (:L1, :L5)
end

@testset "Two groups sharing a band yield one band key" begin
    ts = TrackState(; signals = (legacy_gps = (GPSL1CA(),), another_l1 = (GPSL1CA(),)))
    @test get_band_ids(ts) == (:L1,)
end

@testset "Single-band track! still accepts bare buffer" begin
    sampling_frequency = 4e6Hz
    num_samples = 4000
    prn = 1
    # Integer-typed Hz, next to Float64-typed ones.
    carrier_doppler = 200Hz
    start_code_phase = 100.0

    track_state = TrackState(; signal = GPSL1CA())
    track_state = add_satellite!(
        track_state;
        prn,
        code_phase = start_code_phase,
        carrier_doppler = carrier_doppler - 20Hz,
    )

    samples = _make_signal(
        GPSL1CA,
        prn,
        carrier_doppler,
        start_code_phase,
        num_samples,
        sampling_frequency,
    )

    track!(samples, track_state, sampling_frequency)
    @test abs(get_carrier_doppler(track_state, :default, prn) - carrier_doppler) < 50Hz
end

@testset "Multi-band track! with two synthetic bands" begin
    fs_l1 = 4e6Hz
    fs_l5 = 25e6Hz
    n_l1 = 4000   # 4 MHz × 1 ms
    n_l5 = 25000  # 25 MHz × 1 ms

    prn_l1 = 1
    prn_l5 = 1
    # Integer- and Float64-typed Hz, both end to end.
    cd_l1 = 200Hz
    cd_l5 = -150.0Hz
    cp_l1 = 100.0
    cp_l5 = 250.0

    signal_l1 = _make_signal(GPSL1CA, prn_l1, cd_l1, cp_l1, n_l1, fs_l1)
    signal_l5 = _make_signal(GPSL5I, prn_l5, cd_l5, cp_l5, n_l5, fs_l5)

    track_state =
        TrackState(; signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
    track_state = add_satellite!(
        track_state;
        prn = prn_l1,
        group_key = :legacy_gps_l1,
        code_phase = cp_l1,
        carrier_doppler = cd_l1 - 20Hz,
    )
    track_state = add_satellite!(
        track_state;
        prn = prn_l5,
        group_key = :gps_l5,
        code_phase = cp_l5,
        carrier_doppler = cd_l5 - 20.0Hz,
    )

    multi_band_samples =
        (L1 = BandSamples(signal_l1, fs_l1), L5 = BandSamples(signal_l5, fs_l5))

    track!(multi_band_samples, track_state)

    cd_l1_tracked = get_carrier_doppler(track_state, :legacy_gps_l1, prn_l1)
    cd_l5_tracked = get_carrier_doppler(track_state, :gps_l5, prn_l5)
    @test abs(cd_l1_tracked - cd_l1) < 50Hz
    @test abs(cd_l5_tracked - cd_l5) < 50.0Hz
end

@testset "BandSamples keys mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA())
    m = BandSamples(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((L5 = m,), ts)
end

@testset "Mismatched durations error (no tolerance)" begin
    ts = TrackState(; signals = (legacy_gps_l1 = (GPSL1CA(),), gps_l5 = (GPSL5I(),)))
    ts = add_satellite!(ts; prn = 1, group_key = :legacy_gps_l1, carrier_doppler = 0Hz)
    ts = add_satellite!(ts; prn = 1, group_key = :gps_l5, carrier_doppler = 0Hz)

    # One L5 sample too many.
    m_l1 = BandSamples(zeros(ComplexF64, 4000), 4e6Hz)
    m_l5 = BandSamples(zeros(ComplexF64, 25001), 25e6Hz)
    @test_throws ArgumentError track!((L1 = m_l1, L5 = m_l5), ts)
end

@testset "Equal durations across mixed Float32/Float64 sampling rates" begin
    # Equal durations pass with a Float32 rate, though `4000 / 4e6Hz` and
    # `25000 / 25f6Hz` are not `==`.
    m_l1 = BandSamples(zeros(ComplexF64, 4000), 4e6Hz)
    m_l5 = BandSamples(zeros(ComplexF64, 25000), 25.0f6Hz)
    @test _validate_equal_durations((L1 = m_l1, L5 = m_l5)) === nothing

    # A real mismatch still throws, in seconds rather than Hz^-1.
    m_bad = BandSamples(zeros(ComplexF64, 25001), 25.0f6Hz)
    err = try
        _validate_equal_durations((L1 = m_l1, L5 = m_bad))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test !occursin("Hz^-1", err.msg)
end

@testset "Antenna shape mismatch errors" begin
    ts = TrackState(; signal = GPSL1CA(), num_antennas = NumAntennas(2))
    ts = add_satellite!(ts; prn = 1, carrier_doppler = 0Hz)
    m = BandSamples(zeros(ComplexF64, 4000), 4e6Hz)
    @test_throws ArgumentError track!((L1 = m,), ts)
end

@testset "Per-band antenna counts via SignalGroup entries" begin
    ts = TrackState(;
        signals = (
            legacy_gps_l1 = SignalGroup((GPSL1CA(),); num_antennas = NumAntennas(2)),
            gps_l5 = SignalGroup((GPSL5I(),); num_antennas = NumAntennas(1)),
        ),
    )
    ts = add_satellite!(ts; prn = 1, group_key = :legacy_gps_l1, carrier_doppler = 0Hz)
    ts = add_satellite!(ts; prn = 1, group_key = :gps_l5, carrier_doppler = 0Hz)

    @test get_num_antennas(ts, :legacy_gps_l1, 1) == 2
    @test get_num_antennas(ts, :gps_l5, 1) == 1
end

@testset "Same band, mismatched NumAntennas is rejected" begin
    # Two groups on L1 with different antenna counts.
    @test_throws ArgumentError TrackState(;
        signals = (
            gps = SignalGroup((GPSL1CA(),); num_antennas = NumAntennas(2)),
            galileo = SignalGroup((GalileoE1B(CBOC),); num_antennas = NumAntennas(1)),
        ),
    )
end

@testset "Two groups on the same band with matching NumAntennas is fine" begin
    ts = TrackState(;
        signals = (
            gps = SignalGroup((GPSL1CA(),); num_antennas = NumAntennas(2)),
            galileo = SignalGroup((GalileoE1B(CBOC),); num_antennas = NumAntennas(2)),
        ),
    )
    @test get_band_ids(ts) == (:L1,)
end

end # module
