module FlexibandIntegrationTest

# End-to-end test on the Fraunhofer IIS "Flexiband" reference capture
# `20171017_09-51-43_III-7a_short.zip` (Hanoi, 2017-10-17): GPS L1 C/A and
# Galileo E1B on the L1 band and GPS L5I on the L5 band, at two sample rates.
#
# Opt-in via `ENV["TRACKING_RUN_INTEGRATION_TEST"] = "true"` (a 1.6 GB download,
# cached with Scratch.jl); CI runs it on one matrix job.
#
# ── Capture format (Flexiband / ION SDR-metadata standard) ─────────────────
# The `.usb` file is a stream of fixed 1024-byte USB frames:
#
#     [0x55 0xAA][32-bit big-endian counter][1014-byte payload][0xDEADBEEF]
#      └ 2 bytes ┘└──── 4 bytes ───────────┘└── 169×6 bytes ──┘└─ 4 bytes ─┘
#
# The 1014-byte payload is 169 repeats of a 6-byte "chunk". Each chunk holds
# one base-clock cycle (freqbase = 20.25 MHz) of all three bands, in lump
# order, each band contributing `ratefactor` samples:
#
#     byte:   0        1     2     3     4        5
#     band:   S    | L1_0  L1_1  L1_2  L1_3 |  L5(2 samples)
#     rate:   1×   |        4×              |     2×
#
#   • S  (IRNSS, 20.25 MHz) — 1 sample/chunk, 4-bit I/Q. Not used here.
#   • L1 (81 MHz, GPS L1 C/A + Galileo E1) — 4 samples/chunk, 4-bit signed
#     I/Q packed as I=high nibble, Q=low nibble, two's complement.
#   • L5 (40.5 MHz, GPS L5) — 2 samples/chunk packed into one byte: 2-bit
#     signed I/Q, sample 0 in the high bits (wordshift=Left). The 2-bit
#     levels are sign-magnitude {00,01,10,11} → {+1,+3,-3,-1} (NOT two's
#     complement — that would carry a -0.5 DC bias and spawn false peaks).
#
# Band metadata is in the companion `.xml` (ION SDR-metadata standard):
#   L1: center 1580.0 MHz       → GPS L1 / Galileo E1 IF = 1575.42 - 1580.0  = -4.58   MHz
#   L5: center 1173.546875 MHz  → GPS L5            IF = 1176.45 - 1173.5469 = +2.903125 MHz
# (The XML's L5 translatedfreq, 2.903125e-1, is off by 10×.)

using Test
using Downloads: Downloads
using Scratch: get_scratch!
using ZipFile: ZipFile
using Unitful: Hz, ustrip
using GNSSSignals: GPSL1CA, GPSL5I, GalileoE1B
using Acquisition: acquire, is_detected
using Tracking:
    Tracking, TrackState, BandMeasurement, add_satellite!, track!, get_carrier_doppler
using TrackingLoops: estimate_cn0

# ── Download / cache ───────────────────────────────────────────────────────

const ZIP_URL =
    "https://www.iis.fraunhofer.de/content/dam/iis/de/doc/lv/los/lokalisierung/" *
    "SatNAV/Flexiband%20reference%20Data%2020171017_09-51-43_III-7a_short.zip"
const ZIP_NAME = "20171017_09-51-43_III-7a_short.zip"
const ZIP_SIZE = 1589435475          # exact byte count, used to verify integrity

# The server drops long connections, so fetch byte ranges and assemble them
# once each has its exact size. Returns the path to the cached zip.
function download_capture()
    dir = get_scratch!(Tracking, "flexiband_iii_7a")
    path = joinpath(dir, ZIP_NAME)
    isfile(path) && filesize(path) == ZIP_SIZE && return path

    chunk = 200_000_000
    parts = mktempdir(dir)
    try
        idx, start = 0, 0
        partfiles = String[]
        while start < ZIP_SIZE
            stop = min(start + chunk - 1, ZIP_SIZE - 1)
            want = stop - start + 1
            part = joinpath(parts, "part_$(lpad(idx, 2, '0'))")
            tries = 0
            while !(isfile(part) && filesize(part) == want)
                tries += 1
                tries > 30 && error("Flexiband download: chunk $idx failed after 30 tries")
                try
                    Downloads.download(
                        ZIP_URL,
                        part;
                        headers = ["Range" => "bytes=$start-$stop"],
                    )
                catch err
                    @warn "Flexiband download chunk $idx attempt $tries failed; retrying" err
                end
            end
            push!(partfiles, part)
            start = stop + 1
            idx += 1
        end
        tmp = path * ".part"
        open(tmp, "w") do out
            for p in partfiles
                write(out, read(p))
            end
        end
        filesize(tmp) == ZIP_SIZE ||
            error("Flexiband download: assembled size $(filesize(tmp)) != $ZIP_SIZE")
        mv(tmp, path; force = true)
    finally
        rm(parts; recursive = true, force = true)
    end
    return path
end

# ── Framing / demultiplexing ───────────────────────────────────────────────

const BLOCK = 1024          # bytes per USB frame
const HEADER = 6            # 0x55 0xAA + 4-byte counter
const CHUNK = 6             # bytes per base-clock cycle
const CYCLES = 169          # chunks per frame  (6 + 169*6 + 4 = 1024)
const PAYLOAD = CYCLES * CHUNK
const FREQBASE = 20.25e6    # Hz — base clock (XML <freqbase>)

# Number of frames spanning `seconds` of capture.
nframes(seconds) = ceil(Int, seconds * FREQBASE / CYCLES)

# The payload bytes of up to `nblocks` frames, streamed from the compressed
# `.usb` entry (the full file inflates to 1.9 GB).
#
# Only the first 1_638_984 frames (~13.68 s) are well-formed; the rest, to
# ~15.7 s, is corrupt and demuxes to noise, which was the lock loss of #157.
# Reading stops at the first malformed frame.
function read_payload(zippath, nblocks)
    reader = ZipFile.Reader(zippath)
    try
        usb = only(filter(f -> endswith(f.name, ".usb"), reader.files))
        nblocks = min(nblocks, Int(usb.uncompressedsize ÷ BLOCK))
        payload = Vector{UInt8}(undef, nblocks * PAYLOAD)
        frame = Vector{UInt8}(undef, BLOCK)
        valid = 0
        for b = 0:(nblocks-1)
            read!(usb, frame)
            (frame[1] == 0x55 && frame[2] == 0xAA) || break
            copyto!(payload, b * PAYLOAD + 1, frame, HEADER + 1, PAYLOAD)
            valid = b + 1
        end
        return resize!(payload, valid * PAYLOAD)
    finally
        close(reader)
    end
end

# Frames in the `.usb` entry, well-formed or not.
function num_usb_frames(zippath)
    reader = ZipFile.Reader(zippath)
    try
        usb = only(filter(f -> endswith(f.name, ".usb"), reader.files))
        return Int(usb.uncompressedsize ÷ BLOCK)
    finally
        close(reader)
    end
end

@inline _nibble(x) = (v = Int(x & 0x0f); v ≥ 8 ? v - 16 : v)          # 4-bit two's complement
const _L5_LEVELS = (1.0f0, 3.0f0, -3.0f0, -1.0f0)                      # 2-bit sign-magnitude
@inline _twobit(x) = @inbounds _L5_LEVELS[(x&0x03)+1]

# L1 (81 MHz): chunk bytes 2..5, 4-bit I (high nibble) + 4-bit Q (low nibble).
function demux_l1(payload)
    nchunks = length(payload) ÷ CHUNK
    out = Vector{ComplexF32}(undef, nchunks * 4)
    @inbounds for c = 0:(nchunks-1), k = 0:3
        byte = payload[c*CHUNK+2+k]
        out[c*4+k+1] = ComplexF32(_nibble(byte >> 4), _nibble(byte))
    end
    return out
end

# L5 (40.5 MHz): chunk byte 6 holds two 2-bit I/Q samples, sample 0 in MSBs.
function demux_l5(payload)
    nchunks = length(payload) ÷ CHUNK
    out = Vector{ComplexF32}(undef, nchunks * 2)
    @inbounds for c = 0:(nchunks-1)
        byte = payload[c*CHUNK+6]
        out[c*2+1] = ComplexF32(_twobit(byte >> 6), _twobit(byte >> 4))
        out[c*2+2] = ComplexF32(_twobit(byte >> 2), _twobit(byte))
    end
    return out
end

# ── Band parameters ────────────────────────────────────────────────────────

const FS_L1 = 81.0e6Hz
const FS_L5 = 40.5e6Hz
const IF_L1 = -4.58e6Hz       # GPS L1 / Galileo E1 within the 1580 MHz band
const IF_L5 = 2.903125e6Hz    # GPS L5 within the 1173.546875 MHz band

# ── Test ───────────────────────────────────────────────────────────────────

if get(ENV, "TRACKING_RUN_INTEGRATION_TEST", "false") != "true"
    @info "Skipping Flexiband integration test (downloads a 1.6 GB capture). " *
          "Set ENV[\"TRACKING_RUN_INTEGRATION_TEST\"] = \"true\" to run it."
else
    @testset "Flexiband III-7a multi-band acquire + track" begin
        zippath = download_capture()
        @test filesize(zippath) == ZIP_SIZE

        # Acquire on 40 ms, then track 100 ms chunks to see each loop settle,
        # not merely land near the coarse acquisition (#152).
        n_chunks = 3
        chunk_seconds = 0.100
        payload = read_payload(zippath, nframes(n_chunks * chunk_seconds + 0.005))
        l1 = demux_l1(payload)
        l5 = demux_l5(payload)

        # Exactly 2:1, so both bands span the same time, as `track!` requires.
        @test length(l1) == 2 * length(l5)
        @test length(l1) / ustrip(Hz, FS_L1) ≈ length(l5) / ustrip(Hz, FS_L5)

        n_acq_l1 = round(Int, ustrip(Hz, FS_L1) * 0.040)   # 40 ms
        n_acq_l5 = round(Int, ustrip(Hz, FS_L5) * 0.040)

        # The coherent length sets the Doppler resolution handed to tracking:
        # 4 ms on L1 C/A and E1B, a full 10 ms NH10 period on L5I.
        detected(acqs) = filter(a -> is_detected(a; pfa = 1e-8), acqs)
        acq_l1 = detected(
            acquire(
                GPSL1CA(),
                (@view l1[1:n_acq_l1]),
                FS_L1,
                1:32;
                interm_freq = IF_L1,
                num_coherently_integrated_code_periods = 4,
                num_noncoherent_accumulations = 2,
            ),
        )
        acq_e1 = detected(
            acquire(
                GalileoE1B(),
                (@view l1[1:n_acq_l1]),
                FS_L1,
                1:36;
                interm_freq = IF_L1,
                num_coherently_integrated_code_periods = 1,
                num_noncoherent_accumulations = 4,
            ),
        )
        acq_l5 = detected(
            acquire(
                GPSL5I(),
                (@view l5[1:n_acq_l5]),
                FS_L5,
                1:32;
                interm_freq = IF_L5,
                num_coherently_integrated_code_periods = 10,
                num_noncoherent_accumulations = 1,
            ),
        )

        prns_l1 = sort!([a.prn for a in acq_l1])
        prns_e1 = sort!([a.prn for a in acq_e1])
        prns_l5 = sort!([a.prn for a in acq_l5])
        @info "Flexiband acquisition" GPS_L1CA = prns_l1 Galileo_E1B = prns_e1 GPS_L5I =
            prns_l5

        # A fixed set, so a demux, band or acquisition regression shows. In
        # Oct 2017 only Block IIF carried L5: PRN 6, and PRN 9, which is below
        # the L1 threshold.
        @test prns_l1 == [2, 5, 6, 12, 17, 19]
        @test prns_e1 == [1, 7, 8, 26]
        @test prns_l5 == [6, 9]
        # PRN 6 on both bands: both demuxes are right and time-aligned.
        @test 6 in prns_l1 && 6 in prns_l5

        track_state = TrackState(;
            signals = (
                gps_l1 = (GPSL1CA(),),
                galileo = (GalileoE1B(),),
                gps_l5 = (GPSL5I(),),
            ),
        )
        for a in acq_l1
            track_state = add_satellite!(track_state, a; group = :gps_l1)
        end
        for a in acq_e1
            track_state = add_satellite!(track_state, a; group = :galileo)
        end
        for a in acq_l5
            track_state = add_satellite!(track_state, a; group = :gps_l5)
        end

        # Record each satellite's carrier Doppler after every chunk.
        groups_acqs = ((:gps_l1, acq_l1), (:galileo, acq_e1), (:gps_l5, acq_l5))
        nl1 = round(Int, ustrip(Hz, FS_L1) * chunk_seconds)
        nl5 = nl1 ÷ 2
        doppler_trail =
            Dict((group, a.prn) => Float64[] for (group, acqs) in groups_acqs for a in acqs)
        for c = 1:n_chunks
            l1c = @view l1[((c-1)*nl1+1):(c*nl1)]
            l5c = @view l5[((c-1)*nl5+1):(c*nl5)]
            track_state = track!(
                (
                    L1 = BandMeasurement(l1c, FS_L1, IF_L1),
                    L5 = BandMeasurement(l5c, FS_L5, IF_L5),
                ),
                track_state,
            )
            for (group, acqs) in groups_acqs, a in acqs
                push!(
                    doppler_trail[(group, a.prn)],
                    ustrip(Hz, get_carrier_doppler(track_state, group, a.prn)),
                )
            end
        end

        # Every satellite must stay within half a 250 Hz acquisition bin plus
        # margin of the coarse acquisition (#152), have settled (the last two
        # estimates agree well within a bin), and keep C/N0 above the floor.
        for (group, acqs) in groups_acqs
            for a in acqs
                trail = doppler_trail[(group, a.prn)]
                @test abs(trail[end] - ustrip(Hz, a.carrier_doppler)) < 150
                @test abs(trail[end] - trail[end-1]) < 75
                @test ustrip(estimate_cn0(track_state, group, a.prn)) > 34
            end
        end
    end

    # ── Full-span lock hold, all signals (#157) ─────────────────────────────
    # Every loop must hold lock across the whole valid capture in 0.2 s chunks.
    # Demuxing it at once would take ~13 GB, so each 0.2 s window is demuxed
    # from a view into the ~1.7 GB of packed bytes.
    @testset "All signals hold lock across the full valid capture (#157)" begin
        zippath = download_capture()

        payload = read_payload(zippath, num_usb_frames(zippath))
        total_cycles = length(payload) ÷ CHUNK        # base-clock cycles (6 B each)

        acq_cycles = round(Int, FREQBASE * 0.040)
        acq_bytes = @view payload[1:(acq_cycles*CHUNK)]
        l1_acq = demux_l1(acq_bytes)
        l5_acq = demux_l5(acq_bytes)
        detected(acqs) = filter(a -> is_detected(a; pfa = 1e-8), acqs)
        acq_l1 = detected(
            acquire(
                GPSL1CA(),
                l1_acq,
                FS_L1,
                1:32;
                interm_freq = IF_L1,
                num_coherently_integrated_code_periods = 4,
                num_noncoherent_accumulations = 2,
            ),
        )
        acq_e1 = detected(
            acquire(
                GalileoE1B(),
                l1_acq,
                FS_L1,
                1:36;
                interm_freq = IF_L1,
                num_coherently_integrated_code_periods = 1,
                num_noncoherent_accumulations = 4,
            ),
        )
        acq_l5 = detected(
            acquire(
                GPSL5I(),
                l5_acq,
                FS_L5,
                1:32;
                interm_freq = IF_L5,
                num_coherently_integrated_code_periods = 10,
                num_noncoherent_accumulations = 1,
            ),
        )
        @test sort!([a.prn for a in acq_l1]) == [2, 5, 6, 12, 17, 19]
        @test sort!([a.prn for a in acq_e1]) == [1, 7, 8, 26]
        @test sort!([a.prn for a in acq_l5]) == [6, 9]

        track_state = TrackState(;
            signals = (
                gps_l1 = (GPSL1CA(),),
                galileo = (GalileoE1B(),),
                gps_l5 = (GPSL5I(),),
            ),
        )
        for a in acq_l1
            track_state = add_satellite!(track_state, a; group = :gps_l1)
        end
        for a in acq_e1
            track_state = add_satellite!(track_state, a; group = :galileo)
        end
        for a in acq_l5
            track_state = add_satellite!(track_state, a; group = :gps_l5)
        end

        groups_acqs = ((:gps_l1, acq_l1), (:galileo, acq_e1), (:gps_l5, acq_l5))
        chunk_seconds = 0.2
        # Whole cycles keep both bands time-aligned.
        cycles_per_chunk = round(Int, FREQBASE * chunk_seconds)
        n_chunks = total_cycles ÷ cycles_per_chunk
        @test n_chunks >= 60

        # The corrupt tail took all three signals to ~20-30 dB-Hz at once.
        for c = 1:n_chunks
            window =
                @view payload[((c-1)*cycles_per_chunk*CHUNK+1):(c*cycles_per_chunk*CHUNK)]
            l1c = demux_l1(window)
            l5c = demux_l5(window)
            track_state = track!(
                (
                    L1 = BandMeasurement(l1c, FS_L1, IF_L1),
                    L5 = BandMeasurement(l5c, FS_L5, IF_L5),
                ),
                track_state,
            )
            for (group, acqs) in groups_acqs, a in acqs
                doppler = ustrip(Hz, get_carrier_doppler(track_state, group, a.prn))
                @test abs(doppler - ustrip(Hz, a.carrier_doppler)) < 150
                @test ustrip(estimate_cn0(track_state, group, a.prn)) > 34
            end
        end
    end
end

end # module
