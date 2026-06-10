const GalileoE1BAny = Union{GalileoE1B,GalileoE1B_BOC11}

"""
$(SIGNATURES)

Symbol-sync detector for Galileo E1B (both [`GalileoE1B`](@ref) and the
BOC(1,1) approximation [`GalileoE1B_BOC11`](@ref)).

E1B broadcasts one I/NAV channel symbol per 4 ms primary code period
(250 sym/s; Galileo OS SIS ICD Table 11 — symbol period = primary code
period, Table 15 — 4092 chips at 1.023 Mcps with no secondary code), so
the buffer of primary-block signs is itself the symbol stream — there
is no sub-symbol boundary to find. The detector therefore reports
`found = true` from the very first integration, leaving downstream
consumers (GNSSDecoder.jl) to resolve the residual ±1 polarity
ambiguity via the I/NAV preamble.

Same shape as the GPS L1C-D "1-block-per-symbol" case.
"""
@inline function detect_bit_or_secondary_code_sync(
    ::GalileoE1BAny,
    ::Integer,           # PRN — ignored
    ::Unsigned,
    ::Integer,
)
    SyncResult(true, 0, Int8(+1))
end

# TODO: Very early very late correlator?
function get_default_correlator(::GalileoE1BAny, num_ants::NumAnts = NumAnts(1))
    VeryEarlyPromptLateCorrelator(; num_ants)
end

# 1 channel symbol = 1 primary code period (250 sym/s, 4 ms primary
# period). No sub-symbol boundary to search for — the sync-search buffer
# is dead state at runtime; we still pick a concrete type to keep the
# `BitBuffer{B}` parameter chain stable. UInt8 is the smallest legal
# Unsigned.
@inline get_code_block_buffer_type(::GalileoE1BAny) = UInt8
