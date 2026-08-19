# Code bias — a signal's differential payload group delay against its
# satellite's estimator-driver signal, in seconds.
#
# Only multi-signal discriminator combining reads it. A combined DLL locks the
# satellite-shared `code_phase` to a weighted average of every signal's code
# phase; those code phases differ by the satellite's differential payload group
# delay, so without correcting for it the shared phase sits at a bias that
# *moves as the weights move* — and a downstream consumer
# (`PositionVelocityTime.jl`) then applies the group-delay correction of
# whichever ranging signal it was told, to a phase that is no longer that
# signal's. Referring every passenger's code discriminator to the driver keeps
# `code_phase` meaning "the driver signal's code phase", which is what that
# downstream correction already assumes.
#
# This package deliberately holds no table of per-signal values and no notion of
# where one comes from. A bias may be the difference of two broadcast
# inter-signal corrections, a ground calibration, or a hard zero justified by the
# ICD — deciding that needs constellation knowledge and, usually, a decoded
# navigation message, neither of which belongs in a tracking loop. Every signal
# therefore starts at `nothing` and the caller supplies what it knows, via
# [`set_code_bias!`](@ref).

"""
$(SIGNATURES)

Convert a code bias in seconds to the code-phase offset it produces, in chips, at
the code frequency actually in effect (chip rate plus the satellite's code
Doppler).

Sign convention: a positive bias means the signal sits at a **larger** code phase
than the driver, by `bias · f_code` chips — exactly the amount to subtract from
its code discriminator to refer that discriminator to the driver's code phase.

That is fixed by the loop's own stability rather than by any ICD: a positive
`dll_disc` raises the code frequency, which advances the replica phase, so
`dll_disc` carries the sign of `(true phase − replica phase)`. A passenger whose
code phase exceeds the driver's therefore reads `e_driver + δ`, and `δ` comes off.

Callers deriving a bias from broadcast inter-signal corrections should note that
the ICDs do not agree on *their* sign: IS-GPS-705/800 states a per-component
`ISC_x` that enters the transmit-time correction as `−T_GD + ISC_x`, so the
difference of two GPS ISCs is already in this convention, whereas the BeiDou
B1C/B2a ICDs state `−T_GD_pilot − ISC_data` and so need negating first.
"""
# `bias_seconds` is a bare `Float64` (that is what `set_code_bias!` takes), so it
# is given its unit here before meeting the `Hz`-typed code frequency;
# `uconvert(NoUnits, …)` then hands back a plain `Float64` in chips that the
# discriminators can be corrected with directly.
@inline _group_delay_to_chips(bias_seconds, code_frequency) =
    uconvert(NoUnits, bias_seconds * 1.0s * code_frequency)
