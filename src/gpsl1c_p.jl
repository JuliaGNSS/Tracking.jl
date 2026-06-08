"""
$(SIGNATURES)

Secondary-code sync detector for GPS L1C-P.

L1C-P broadcasts a per-PRN 1800-chip overlay code (IS-GPS-800G §3.2.2.1.2)
on top of the 10-ms primary code, giving an 18-second cycle. To lock the
overlay we wait for the sliding `code_block_bits` window to fill to 1800
primary periods, then run a single 1800-phase shifted Hamming-distance
sweep against the PRN's known overlay pattern.

Returns `SyncResult(false, 0, 0)` until 1800 blocks have been buffered.
Once that horizon is reached, the sweep picks the rotation `k ∈ 0:1799`
whose Hamming distance to the overlay (or its negation) is minimal; if
that distance is within the 2 % tolerance (≤ 36 errors) the result is
`SyncResult(true, k, ±1)`. The `phase` field encodes the secondary-chip
offset of the *first* sample currently in the buffer; downstream code
uses it to seed the shared `sat.code_phase`.
"""
function is_upcoming_integration_new_bit(
    signal::GPSL1C_P,
    prn::Integer,
    code_block_bits::UInt1800,
    num_code_blocks::Integer,
)
    num_code_blocks < 1800 && return SyncResult(false, 0, Int8(0))
    overlay = _pack_overlay(signal, prn)
    # Tolerance is a percentage of the 1800-chip overlay window. The
    # default 2.5 % discretizes to floor(0.025 × 1800) = 45 bit-flips
    # allowed. Adjustable via
    # `get_bit_edge_or_secondary_code_tolerance(::GPSL1C_P)`.
    max_errors = floor(Int, get_bit_edge_or_secondary_code_tolerance(signal) * 1800)
    _full_phase_search(code_block_bits, overlay, max_errors)
end

# Pack the ±1 secondary-code column for `prn` into a UInt1800 where bit
# `k` is `1` iff the overlay chip at secondary-chip index `k` is `+1`.
# Each rebuild walks 1800 chips and is hit once per sat at sync time, so
# the cost is negligible.
@inline function _pack_overlay(signal::GPSL1C_P, prn::Integer)
    codes = signal.overlay_codes        # 1800 × 63 Int8 matrix
    x = UInt1800(0)
    @inbounds for k in 1:1800
        if codes[k, prn] > 0
            x |= UInt1800(1) << (k - 1)
        end
    end
    x
end

"""
$(SIGNATURES)

Full 1800-phase shifted Hamming-distance search over the L1C-P overlay.

For each candidate rotation `k ∈ 0:1799`, rotate `received` left by `k`
bits within its 1800-bit window and count differing bits against
`overlay`. Tracks the best positive- and negative-polarity matches in
the same pass. Returns `SyncResult(found, phase, polarity)` where
`found = true` iff the best Hamming distance does not exceed
`max_errors`.

The 1800-bit container is exact-width (no padding), so the rotate uses
the natural `<<` / `>>` semantics of `UInt1800` and no mask is needed
on the XOR. Benchmark: ~71 μs for the full sweep on Zen 4.
"""
@inline function _full_phase_search(received::UInt1800, overlay::UInt1800, max_errors::Int)
    best_k = 0
    best_dist = 1800
    best_pol = Int8(0)
    @inbounds for k in 0:1799
        # Rotate-left by k within the 1800-bit container. Special-case
        # k = 0 because `received >> 1800` is undefined for an
        # exact-width 1800-bit integer.
        shifted = k == 0 ? received : ((received << k) | (received >> (1800 - k)))
        d_pos = count_ones(shifted ⊻ overlay)
        d_neg = 1800 - d_pos
        if d_pos < best_dist
            best_dist = d_pos
            best_k = k
            best_pol = Int8(+1)
        end
        if d_neg < best_dist
            best_dist = d_neg
            best_k = k
            best_pol = Int8(-1)
        end
    end
    if best_dist <= max_errors
        return SyncResult(true, best_k, best_pol)
    else
        return SyncResult(false, 0, Int8(0))
    end
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
