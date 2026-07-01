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
    _detect_secondary_code_sync(signal, prn, code_block_bits, num_code_blocks)
end

# Pack the ±1 overlay-code column for `prn` into a UInt1800 in the
# newest-first order [`_packed_secondary_code`](@ref) requires: bit `i` is
# `1` iff overlay chip `1799 - i` is `+1`. This makes the packed reference
# directly comparable to the prompt buffer when the most recent 1800
# blocks end on the overlay's last chip — see [`_secondary_code_search`](@ref).
# The reference is rebuilt on every detector call once the 1800-block
# window has filled, until lock; each rebuild walks 1800 chips, which is
# small next to the 1800-phase Hamming sweep that follows it. Under
# weak-signal conditions (many detector calls before lock) a per-PRN cache
# would shave that rebuild, but the sweep still dominates, so it's not
# worth the state.
@inline function _packed_secondary_code(::Type{UInt1800}, signal::GPSL1C_P, prn::Integer)
    codes = signal.overlay_codes        # 1800 × 63 Int8 matrix
    x = UInt1800(0)
    @inbounds for k = 1:1800
        if codes[k, prn] > 0
            x |= UInt1800(1) << (1800 - k)
        end
    end
    x
end

# `get_default_correlator(::GPSL1C_P)` — the VeryEarlyPromptLate BOC default —
# is defined jointly with GPS L1C-D in `l1c_d.jl` (one `Union` method).

# 1800-chip overlay search needs an exact-width 1800-bit container. The
# `UInt1800` alias is defined in the top-level Tracking module via
# `BitIntegers.@define_integers 1800` and is what `BitBuffer{B}` carries
# for L1C-P throughout the tracker.
@inline get_code_block_buffer_type(::GPSL1C_P) = UInt1800
