# Differential group delay — a signal's differential payload group delay against its
# satellite's estimator-driver signal, in seconds. Read by multi-signal code-discriminator
# combining and by nothing else; see `set_differential_group_delay!` for the contract and
# the manual's "Multi-signal discriminator combining" for why the value is supplied rather
# than derived here.

# Normalize whatever the caller supplied to the field's own type, `typeof(1.0s)`.
#
# Every rejection happens *here*, on a `Number`-wide entry point, rather than by narrowing
# the callers' argument types: a signature that admits only `Union{Real,Unitful.Time}` turns
# the likely mistakes into a `MethodError` whose text is a page of `TrackState` type
# parameters. So the callers take `Maybe{Number}` and a wrong dimension gets a sentence.
_as_differential_group_delay(::Nothing) = nothing
_as_differential_group_delay(delay::Unitful.Time) = float(uconvert(s, delay))
_as_differential_group_delay(delay::Real) = _throw_unitless_differential_group_delay(delay)
_as_differential_group_delay(delay::Number) =
    _throw_wrong_dimension_differential_group_delay(delay)

# A bare number is refused rather than read as seconds: every dimensioned quantity in this
# package carries its unit, and a realistic value here is sub-nanosecond — the scale at
# which a silently assumed unit turns into a metre.
@noinline _throw_unitless_differential_group_delay(delay) = throw(
    ArgumentError(
        "a differential group delay is a time and needs its unit: got the bare " *
        "number $delay. " *
        "Write `$(delay)s` for seconds, or e.g. `$(delay)u\"ns\"` " *
        "(`using Unitful`) — a broadcast inter-signal correction arrives in " *
        "seconds, so `isc_difference * 1.0s`.",
    ),
)

# A length is the *likely* wrong input rather than an exotic one: "code bias" in the
# SSR/PPP world (Galileo HAS, IGS, and GNSSDecoder.jl's own HAS decoder) is a per-signal
# pseudorange correction in metres, so a caller arriving from that side holds metres. Name
# the conversion, and the datum trap that comes with it, instead of refusing blankly.
@noinline function _throw_wrong_dimension_differential_group_delay(delay)
    hint =
        dimension(delay) == Unitful.𝐋 ?
        " A per-signal code bias in metres (the SSR/PPP sense) converts with the " *
        "speed of light — `$delay / Unitful.c0` — but check the datum first: this " *
        "field is a differential against the estimator-driver signal (`signals[1]`), " *
        "whereas an SSR code bias is referenced to the product's own clock datum and " *
        "is not a difference of two signals at all." : ""
    throw(
        ArgumentError(
            "a differential group delay is a time: got $delay, which has dimension " *
            "$(dimension(delay)). Supply it in seconds — `1.2e-9s`, `-0.3u\"ns\"` — or " *
            "`nothing` to mark it unknown." *
            hint,
        ),
    )
end

"""
$(SIGNATURES)

Convert a differential group delay (a time) to the code-phase offset it produces,
in chips, at
the code frequency actually in effect (chip rate plus the satellite's code
Doppler).

Sign convention: a positive delay means the signal sits at a **larger** code phase
than the driver, by `delay · f_code` chips — exactly the amount to subtract from
its code discriminator to refer that discriminator to the driver's code phase.

That is fixed by the loop's own stability rather than by any ICD: a positive
`dll_disc` raises the code frequency, which advances the replica phase, so
`dll_disc` carries the sign of `(true phase − replica phase)`. A passenger whose
code phase exceeds the driver's therefore reads `e_driver + δ`, and `δ` comes off.

See [`set_differential_group_delay!`](@ref) for how to derive `δ` from broadcast
inter-signal corrections, whose own sign conventions differ between the GPS and
BeiDou ICDs.
"""
# Both arguments carry their units, so this is a plain multiply; `uconvert(NoUnits,
# …)` strips the (already dimensionless) product back to a `Float64` in chips that
# the discriminators can be corrected with directly.
@inline _group_delay_to_chips(delay, code_frequency) =
    uconvert(NoUnits, delay * code_frequency)
