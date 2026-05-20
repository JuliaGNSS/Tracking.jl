"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL1C_P.

GPS L1C-P is the *pilot* channel of L1C: it carries no navigation data
(`get_data_frequency` = 0 Hz) and instead has a 1800-bit prn-specific
overlay (Neuman-Hofman) secondary code with `get_secondary_code_length` =
1800. A real receiver would search for that overlay's phase via cross-
correlation against the 1800-bit pattern to recover the secondary-code
boundary, after which longer coherent integrations across the overlay
period become possible.

That overlay-code sync is not yet implemented in Tracking.jl — the
existing `BitBuffer.code_block_buffer::UInt128` only holds 128
primary-code-block signs, which is not enough for the 1800-bit search.
Returning `false` here keeps `bit_buffer.found = false`, which makes the
inner loop integrate one primary code period (10 ms) at a time forever.
That's enough to track the pilot's PLL/DLL through `signals[1]` of a
multi-signal sat without extended coherent integration.

Tracked status: pilot tracking works, navigation-bit recovery via the
secondary-code overlay is deferred to a follow-up.
"""
function is_upcoming_integration_new_bit(
    ::GPSL1C_P,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    # Step 4 of the sync-detection redesign replaces this with the
    # 1800-phase overlay search. Until then the pilot's bit buffer never
    # locks and the inner loop integrates one primary code period (10 ms)
    # at a time.
    SyncResult(false, 0, Int8(0))
end

function get_default_correlator(gpsl1c_p::GPSL1C_P, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# Placeholder until step 4 of the sync-detection redesign brings
# BitIntegers' `UInt1800`. Tracking still works with `bit_buffer.found =
# false` forever — the 1800-chip overlay search is not yet wired.
@inline get_code_block_buffer_type(::GPSL1C_P) = UInt64
