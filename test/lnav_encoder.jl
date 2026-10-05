# A GPS L1 C/A LNAV encoder (IS-GPS-200, section 20.3) — the one TrackingLoops'
# tests use (its `test/lnav_encoder.jl`) — for the simulations: it
# encodes a decoded ephemeris back into the subframes a satellite broadcasts, so a
# simulated satellite can carry real navigation bits and the receiver under test
# decodes them itself. Subframes 1–3 carry the clock and ephemeris; subframes 4 and 5
# carry a page the decoder ignores. Every word has its parity, and words 2 and 10 the
# two solved-for bits that zero their last parity bits.
using GNSSDecoder: GPSL1CAData

const LNAV_PREAMBLE = 0b10001011

# The unsigned `num_bits`-bit field of `value / scale`, which must be an exact
# multiple of it (a decoded value always is), as two's complement when `signed`.
function lnav_field(value, scale, num_bits; signed = false)
    value = something(value)
    count = round(Int64, value / scale)
    @assert isapprox(count * scale, value; rtol = 1e-12, atol = abs(scale) * 1e-6) "$value is no multiple of $scale"
    if signed
        @assert -(Int64(1) << (num_bits - 1)) <= count < Int64(1) << (num_bits - 1)
        count = mod(count, Int64(1) << num_bits)
    else
        @assert 0 <= count < Int64(1) << num_bits
    end
    UInt32(count)
end

# The 24 data bits of a word from `(field, num_bits)` pairs, most significant first;
# the bits left over are zero.
function lnav_word(fields...)
    word = UInt32(0)
    used = 0
    for (field, num_bits) in fields
        word =
            (word << num_bits) |
            (UInt32(something(field)) & ((UInt32(1) << num_bits) - UInt32(1)))
        used += num_bits
    end
    @assert used <= 24
    word << (24 - used)
end

lnav_bit(data, n) = (data >> (24 - n)) & UInt32(1) == 1

# The six parity bits of the source data `data` (24 bits), after the previous word's
# last two transmitted bits (Table 20-XIV).
function lnav_parity(data, D29, D30)
    d(n) = lnav_bit(data, n)
    D25 =
        D29 ⊻ d(1) ⊻ d(2) ⊻ d(3) ⊻ d(5) ⊻ d(6) ⊻ d(10) ⊻ d(11) ⊻ d(12) ⊻ d(13) ⊻ d(14) ⊻
        d(17) ⊻ d(18) ⊻ d(20) ⊻ d(23)
    D26 =
        D30 ⊻ d(2) ⊻ d(3) ⊻ d(4) ⊻ d(6) ⊻ d(7) ⊻ d(11) ⊻ d(12) ⊻ d(13) ⊻ d(14) ⊻ d(15) ⊻
        d(18) ⊻ d(19) ⊻ d(21) ⊻ d(24)
    D27 =
        D29 ⊻ d(1) ⊻ d(3) ⊻ d(4) ⊻ d(5) ⊻ d(7) ⊻ d(8) ⊻ d(12) ⊻ d(13) ⊻ d(14) ⊻ d(15) ⊻
        d(16) ⊻ d(19) ⊻ d(20) ⊻ d(22)
    D28 =
        D30 ⊻ d(2) ⊻ d(4) ⊻ d(5) ⊻ d(6) ⊻ d(8) ⊻ d(9) ⊻ d(13) ⊻ d(14) ⊻ d(15) ⊻ d(16) ⊻
        d(17) ⊻ d(20) ⊻ d(21) ⊻ d(23)
    D29_ =
        D30 ⊻ d(1) ⊻ d(3) ⊻ d(5) ⊻ d(6) ⊻ d(7) ⊻ d(9) ⊻ d(10) ⊻ d(14) ⊻ d(15) ⊻ d(16) ⊻
        d(17) ⊻ d(18) ⊻ d(21) ⊻ d(22) ⊻ d(24)
    D30_ =
        D29 ⊻ d(3) ⊻ d(5) ⊻ d(6) ⊻ d(8) ⊻ d(9) ⊻ d(10) ⊻ d(11) ⊻ d(13) ⊻ d(15) ⊻ d(19) ⊻
        d(22) ⊻ d(23) ⊻ d(24)
    (D25, D26, D27, D28, D29_, D30_)
end

# The 300 transmitted bits of a subframe from its ten 24-bit data words; words 2 and 10
# get bits 23–24 solved so that their parity bits 29 and 30 are zero.
function lnav_subframe_bits(words)
    bits = Bool[]
    D29, D30 = false, false
    for (index, data) in enumerate(words)
        if index == 2 || index == 10
            data &= ~UInt32(0b11)
            for t = UInt32(0):UInt32(3)
                parity = lnav_parity(data | t, D29, D30)
                if !parity[5] && !parity[6]
                    data |= t
                    break
                end
            end
        end
        parity = lnav_parity(data, D29, D30)
        for n = 1:24
            push!(bits, lnav_bit(data, n) ⊻ D30)
        end
        append!(bits, parity)
        D29, D30 = parity[5], parity[6]
    end
    @assert length(bits) == 300
    bits
end

# The TLM and HOW words of the subframe that starts at time of week `tow` (a multiple
# of 6 s): the HOW carries the count of the *next* subframe's start.
function lnav_tlm_and_how(tow, subframe_id)
    tlm = lnav_word((LNAV_PREAMBLE, 8), (0, 14), (0, 1), (0, 1))
    count = mod(div(tow, 6) + 1, 100800)
    how = lnav_word((count, 17), (0, 1), (0, 1), (subframe_id, 3))
    tlm, how
end

const PI_GPS = 3.1415926535898

function lnav_subframe_words(data::GPSL1CAData, tow, subframe_id)
    tlm, how = lnav_tlm_and_how(tow, subframe_id)
    if subframe_id == 1
        iodc = something(data.IODC)
        words = (
            tlm,
            how,
            lnav_word(
                (mod(data.WN, 1024), 10),
                (something(data.code_on_L2, 1), 2),
                (data.URA_index, 4),
                (data.sv_health, 6),
                (iodc >> 8, 2),
            ),
            lnav_word((something(data.L2_P_data_flag, false), 1)),
            lnav_word(),
            lnav_word(),
            lnav_word((0, 16), (lnav_field(data.T_GD, 2.0^-31, 8; signed = true), 8)),
            lnav_word((iodc & 0xff, 8), (lnav_field(data.t_0c, 16, 16), 16)),
            lnav_word(
                (lnav_field(data.a_f2, 2.0^-55, 8; signed = true), 8),
                (lnav_field(data.a_f1, 2.0^-43, 16; signed = true), 16),
            ),
            lnav_word((lnav_field(data.a_f0, 2.0^-31, 22; signed = true), 22)),
        )
    elseif subframe_id == 2
        M_0 = lnav_field(data.M_0, PI_GPS * 2.0^-31, 32; signed = true)
        e = lnav_field(data.e, 2.0^-33, 32)
        sqrt_A = lnav_field(data.sqrt_A, 2.0^-19, 32)
        words = (
            tlm,
            how,
            lnav_word(
                (data.IODE_Sub_2, 8),
                (lnav_field(data.C_rs, 2.0^-5, 16; signed = true), 16),
            ),
            lnav_word(
                (lnav_field(data.Δn, PI_GPS * 2.0^-43, 16; signed = true), 16),
                (M_0 >> 24, 8),
            ),
            lnav_word((M_0 & 0xffffff, 24)),
            lnav_word(
                (lnav_field(data.C_uc, 2.0^-29, 16; signed = true), 16),
                (e >> 24, 8),
            ),
            lnav_word((e & 0xffffff, 24)),
            lnav_word(
                (lnav_field(data.C_us, 2.0^-29, 16; signed = true), 16),
                (sqrt_A >> 24, 8),
            ),
            lnav_word((sqrt_A & 0xffffff, 24)),
            lnav_word(
                (lnav_field(data.t_0e, 16, 16), 16),
                (something(data.fit_interval, false), 1),
                (div(something(data.AODO, 0), 900), 5),
            ),
        )
    elseif subframe_id == 3
        Ω_0 = lnav_field(data.Ω_0, PI_GPS * 2.0^-31, 32; signed = true)
        i_0 = lnav_field(data.i_0, PI_GPS * 2.0^-31, 32; signed = true)
        ω = lnav_field(data.ω, PI_GPS * 2.0^-31, 32; signed = true)
        words = (
            tlm,
            how,
            lnav_word(
                (lnav_field(data.C_ic, 2.0^-29, 16; signed = true), 16),
                (Ω_0 >> 24, 8),
            ),
            lnav_word((Ω_0 & 0xffffff, 24)),
            lnav_word(
                (lnav_field(data.C_is, 2.0^-29, 16; signed = true), 16),
                (i_0 >> 24, 8),
            ),
            lnav_word((i_0 & 0xffffff, 24)),
            lnav_word((lnav_field(data.C_rc, 2.0^-5, 16; signed = true), 16), (ω >> 24, 8)),
            lnav_word((ω & 0xffffff, 24)),
            lnav_word((lnav_field(data.Ω_dot, PI_GPS * 2.0^-43, 24; signed = true), 24)),
            lnav_word(
                (data.IODE_Sub_3, 8),
                (lnav_field(data.i_dot, PI_GPS * 2.0^-43, 14; signed = true), 14),
            ),
        )
    else
        # Subframes 4 and 5: data ID 01 and a page ID the decoder does not read.
        words = (tlm, how, lnav_word((1, 2), (0, 6)), ntuple(_ -> lnav_word(), 7)...)
    end
    words
end

"""
    LNAVStream(data::GPSL1CAData)

The LNAV bit stream of a satellite broadcasting the clock and ephemeris `data`:
`lnav_bit(stream, k)` is the bit transmitted over `[k, k + 1) / 50` seconds of GPS
week, `true` for a one. Subframes are encoded on demand and cached.
"""
struct LNAVStream
    data::GPSL1CAData
    subframes::Dict{Int,Vector{Bool}}
end

LNAVStream(data::GPSL1CAData) = LNAVStream(data, Dict{Int,Vector{Bool}}())

function lnav_bit(stream::LNAVStream, k::Integer)
    subframe = fld(k, 300)
    bits = get!(stream.subframes, subframe) do
        tow = 6 * subframe
        lnav_subframe_bits(lnav_subframe_words(stream.data, tow, mod(subframe, 5) + 1))
    end
    bits[k-300*subframe+1]
end

# The time of week a subframe 1 starts at, at or after `tow`.
next_subframe1_start(tow) = 30 * cld(tow, 30)
