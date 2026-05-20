"""
$(SIGNATURES)

Checks if upcoming integration is a new bit for GPSL1C_D.

GPS L1C-D has a 10 ms primary code period (10230 chips at 1.023 Mcps) and a
50 Hz data rate — so each navigation bit spans **two** primary code blocks,
not one. Bit-sync therefore requires a real preamble search (the L1C CNAV-2
subframe carries a Barker-like preamble) which Tracking.jl does not yet
implement. Returning `false` here keeps `bit_buffer.found = false`, which
makes the inner loop integrate one primary code period (10 ms) at a time
forever. PLL/DLL tracking still works at 10-ms boundaries; bit recovery is
deferred to a follow-up.
"""
function is_upcoming_integration_new_bit(
    ::GPSL1C_D,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
)
    # Step 5 of the sync-detection redesign replaces this with
    # `SyncResult(true, 0, Int8(+1))` — L1C-D broadcasts one channel
    # symbol per primary code period (100 sps), so there's no sub-symbol
    # boundary to find. Returning `found = false` here keeps tracking
    # working at 10 ms primary-code boundaries until that change lands.
    SyncResult(false, 0, Int8(0))
end

function get_default_correlator(gpsl1c_d::GPSL1C_D, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# 1 channel symbol = 1 primary code period (100 sps, 10 ms primary
# period). No sub-symbol boundary to search for — the sync-search buffer
# is dead state at runtime; we still pick a concrete type to keep the
# `BitBuffer{B}` parameter chain stable. UInt8 is the smallest legal
# Unsigned.
@inline get_code_block_buffer_type(::GPSL1C_D) = UInt8
