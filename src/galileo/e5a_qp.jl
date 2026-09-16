# Galileo E5a-QP is the E5a *quasi-pilot*: the short-code acquisition aid added
# in OS SIS ICD Issue 2.2 (§2.3.1.4), "designed to enable low complexity
# acquisition capability of Galileo signals". It is a dataless BPSK(5) component
# at 5.115 Mcps on the 1176.45 MHz carrier (`L5`) whose primary code is only 330
# chips — a 64.5 µs period, repeated 31 times within 2 ms with no overlay, so
# there is no secondary code and no data symbol.
#
# Two consequences shape everything below.
#
#   1. There is nothing to synchronise. No data bit, no overlay chip, so no
#      boundary to search for — every primary block boundary is equivalent. The
#      detector says so immediately (see below), the same "the symbol grid *is*
#      the code block grid" answer Galileo E1B / E6-B, GPS L1C-D / L2CM and
#      BeiDou B2b-I give, except that here it is the coherent-integration grid
#      rather than a symbol grid.
#
#   2. One primary code block is far too short to run a loop on. This package
#      integrates in whole primary code blocks and sizes the carrier loop off
#      the primary period, so a 64.5 µs block would mean a 15.5 kHz loop update
#      rate against a 279 Hz reference bandwidth. E5a-QP therefore integrates a
#      whole 31-block code cycle — 10230 chips, exactly 2 ms — by default, which
#      puts the loop at 500 Hz and, through the estimator's automatic 1/N
#      bandwidth scaling, at a 9 Hz effective carrier bandwidth. That is the
#      "short-period integration policy" of issue #236, expressed as the two
#      integration-length traits rather than as a special case in `track`.
#
# What E5a-QP is *not* is a third component of the E5a composite. It runs at
# half the E5a-I/E5a-Q chip rate, so it cannot share a satellite's
# `SignalGroup` with them (that would need issue #151), and the ICD claims no
# carrier phase relation to the E5a-I reference, so `get_carrier_phase_offset`
# is the bare 0.0 default rather than a quadrature offset. It is acquired and
# tracked on its own.
#
# Its acquisition-aid role is a property of the 330-chip code, not of an API
# here: a 330-chip search is 31× cheaper than E5a-I's 10230-chip one, and the
# recovered carrier Doppler transfers directly, while the code phase pins E5a-I
# only modulo 64.5 µs — 31 residual hypotheses across the 2 ms in which the two
# grids realign, instead of 10230 chip bins. Tracking E5a-QP itself yields code
# phase, carrier phase and Doppler observables like any other pilot; being
# dataless it contributes no navigation message, on its own or for E5a.
#
# Note the ICD defines 40 primary codes (SVID `n` uses code `n`), so PRNs 1-40
# are supported, and that E5a-QP is transmitted by a subset of the constellation
# and is to be superseded eventually. Both are upstream properties of the code
# table, not something the tracker screens for.

"""
$(SIGNATURES)

"Sync" detector for Galileo E5a-QP — reports `found = true` from the very first
integration.

E5a-QP is dataless (`get_data_frequency == 0`) and carries no overlay
(`get_secondary_code` is `NoSecondaryCode`), so there is no bit edge and no
secondary-code phase to locate: every 330-chip primary block boundary is
equivalent, and the only thing "sync" gates here is the switch from the
single-block pre-sync integration to the whole-cycle one
([`max_num_code_blocks_to_integrate`](@ref)). Reporting the lock immediately is
therefore correct rather than optimistic — unlike the GPS L2CL case, which is
the same dataless/overlay-free shape but whose 1.5 s primary period is already a
whole coherent integration, so there is nothing for a lock to unlock.

It is therefore the same "fires immediately, nothing to find" body the
one-symbol-per-code-period signals use, [`_detect_symbol_is_code_block_sync`](@ref):
`phase` is 0 (no secondary code to index) and `polarity` is `+1` (with no data
and no overlay there is no sign convention to recover). Returns
[`SyncResult`](@ref).
"""
@inline detect_bit_or_secondary_code_sync(
    signal::GalileoE5aQP,
    prn::Integer,
    code_block_bits::Unsigned,
    num_code_blocks::Integer,
) = _detect_symbol_is_code_block_sync(signal, prn, code_block_bits, num_code_blocks)

# Plain BPSK(5) (`LOC`) on the L5 carrier, so the C/A-style EarlyPromptLate
# default applies — as for both E5a components.
function get_default_correlator(::GalileoE5aQP, num_ants::NumAnts = NumAnts(1))
    EarlyPromptLateCorrelator(; num_ants)
end

# No bit and no secondary code to search for, so the packed sign window is dead
# state; a concrete type is still needed to keep the `BitBuffer{B}` parameter
# chain stable, and UInt8 is the smallest legal Unsigned — as for GPS L2CM /
# L2CL and BeiDou B2b-I.
@inline get_code_block_buffer_type(::GalileoE5aQP) = UInt8

# One coherent integration is one whole 31-block code cycle: 31 × 330 = 10230
# chips at 5.115 Mcps, exactly the 2 ms the ICD describes the code as repeating
# 31 times within (Galileo OS SIS ICD v2.2 Table 15). The generic ceiling — the
# data-bit period, or the secondary-code period for a pilot — degenerates to a
# single block here because E5a-QP has neither, which is precisely the case the
# trait exists to override.
@inline max_num_code_blocks_to_integrate(::GalileoE5aQP) = 31
# …and it is also the default, because a single 64.5 µs block is not a usable
# integration for any loop. Every other signal keeps the single-block default.
@inline default_num_code_blocks_to_integrate(::GalileoE5aQP) = 31
