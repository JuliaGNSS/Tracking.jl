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
    gpsl1c_d::GPSL1C_D,
    code_block_bits,
    num_code_blocks,
)
    false
end

function get_default_correlator(gpsl1c_d::GPSL1C_D, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end
