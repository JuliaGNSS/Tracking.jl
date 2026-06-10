"""
$(SIGNATURES)

Secondary-code sync detector for GPS L1C-P.

L1C-P broadcasts a per-PRN 1800-chip overlay code (IS-GPS-800G §3.2.2.1.2)
on top of the 10-ms primary code, giving an 18-second cycle. To lock the
overlay we wait for the sliding `code_block_bits` window to fill to 1800
primary periods, then run a single 1800-phase shifted Hamming-distance
sweep against the PRN's known overlay pattern.

Returns `SyncResult(false, 0, 0)` until 1800 blocks have been buffered.
Once that horizon is reached the generic [`_secondary_code_search`](@ref)
rotation sweep picks the alignment whose Hamming distance to the overlay
(or its negation) is minimal; if that distance is within the 2.5 %
tolerance (≤ 45 errors) it reports `SyncResult(true, phase, ±1)` where
`phase` is the secondary-chip offset of the *upcoming* integration, which
downstream code uses to anchor the shared `sat.code_phase`.
"""
function detect_bit_or_secondary_code_sync(
    signal::GPSL1C_P,
    prn::Integer,
    code_block_bits::UInt1800,
    num_code_blocks::Integer,
)
    secondary_code_length = get_secondary_code_length(signal)   # 1800
    num_code_blocks < secondary_code_length && return SyncResult(false, 0, Int8(0))
    reference = _pack_overlay(signal, prn)
    # Tolerance is a percentage of the 1800-chip overlay window. The
    # default 2.5 % discretizes to floor(0.025 × 1800) = 45 bit-flips
    # allowed. Adjustable via
    # `get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P)`.
    max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * secondary_code_length)
    _secondary_code_search(code_block_bits, reference, secondary_code_length, max_errors)
end

# Pack the ±1 secondary-code column for `prn` into a UInt1800 in the same
# newest-first order the prompt buffer fills (bit 0 = the most recent
# block): bit `i` is `1` iff overlay chip `1799 - i` is `+1`. This makes
# the packed reference directly comparable to the prompt buffer when the
# most recent 1800 blocks end on the overlay's last chip — see
# [`_secondary_code_search`](@ref). Each rebuild walks 1800 chips and is
# hit once per sat at sync time, so the cost is negligible.
@inline function _pack_overlay(signal::GPSL1C_P, prn::Integer)
    codes = signal.overlay_codes        # 1800 × 63 Int8 matrix
    x = UInt1800(0)
    @inbounds for k in 1:1800
        if codes[k, prn] > 0
            x |= UInt1800(1) << (1800 - k)
        end
    end
    x
end

# GPS L1C-P is TMBOC(6,1,4/33): 29/33 of the spreading symbols are BOC(1,1)
# and 4/33 are BOC(6,1), so its autocorrelation has the same BOC(1,1)-style
# narrow main peak with side-lobes (slightly sharpened by the BOC(6,1) chips).
# As for L1C-D, the C/A-style 0.5-chip early-late spacing would land the taps
# on the side-lobes and bias the DLL discriminator; a narrow 0.1-chip spacing
# keeps them on the main peak.
function get_default_correlator(gpsl1c_p::GPSL1C_P, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants, preferred_early_late_to_prompt_code_shift = 0.1)
end

# 1800-chip overlay search needs an exact-width 1800-bit container. The
# `UInt1800` alias is defined in the top-level Tracking module via
# `BitIntegers.@define_integers 1800` and is what `BitBuffer{B}` carries
# for L1C-P throughout the tracker.
@inline get_code_block_buffer_type(::GPSL1C_P) = UInt1800
