"""
$(SIGNATURES)

One band's incoming sample buffer plus the front-end metadata needed to
process it. Bundles `samples` with the `sampling_frequency` and
`intermediate_frequency` they were captured at — these are inseparable
in practice, and the bundle removes the chance of mismatched parallel
NamedTuples in a multi-band `track` call.

Fields:
- `samples::S`: complex sample buffer (`Vector` for one antenna, `Matrix`
  with rows = samples and columns = antennas for an antenna array)
- `sampling_frequency::F`: the buffer's sample rate (e.g. `4e6Hz`)
- `intermediate_frequency::F`: the band's IF (defaults to `0.0Hz`)

In a multi-band call, one `Measurement` is built per band; a NamedTuple
of `Measurement`s keyed by band feeds `track`. For the single-band case
a plain buffer + scalar sample-rate keeps working unchanged.

```julia
Measurement(buf, 4e6Hz)                              # IF defaults to 0.0Hz
Measurement(buf, 4e6Hz, 1.575e6Hz)                   # explicit IF
Measurement(buf; sampling_frequency = 4e6Hz)         # kwarg form
```
"""
struct Measurement{S<:AbstractVecOrMat,F}
    samples::S
    sampling_frequency::F
    intermediate_frequency::F
end

# Promotes the numeric type so `sampling_frequency` and
# `intermediate_frequency` end up the same concrete type (`F`), e.g. an
# integer-typed IF alongside a `Float64` sampling frequency.
function Measurement(
    samples::AbstractVecOrMat,
    sampling_frequency,
    intermediate_frequency,
)
    Measurement(samples, promote(sampling_frequency, intermediate_frequency)...)
end

# Positional constructor with IF defaulting to 0.0Hz.
function Measurement(samples::AbstractVecOrMat, sampling_frequency)
    intermediate_frequency = zero(sampling_frequency)
    Measurement(samples, sampling_frequency, intermediate_frequency)
end

# Kwarg constructor. `intermediate_frequency` defaults to zero of the
# same type as `sampling_frequency` so the struct's `F` stays unified.
function Measurement(
    samples::AbstractVecOrMat;
    sampling_frequency,
    intermediate_frequency = zero(sampling_frequency),
)
    Measurement(samples, sampling_frequency, intermediate_frequency)
end

@inline get_samples(m::Measurement) = m.samples
@inline get_sampling_frequency(m::Measurement) = m.sampling_frequency
@inline get_intermediate_frequency(m::Measurement) = m.intermediate_frequency
@inline get_num_samples(m::Measurement) = get_num_samples(m.samples)

"""
$(SIGNATURES)

Map a GNSSSignals `Band` instance to the `Symbol` used as the
NamedTuple key in a multi-band measurements collection. Singleton
dispatch — fold to a compile-time constant when `band` has a concrete
type, so the per-call NamedTuple lookup is free.

Concrete bands define one method each. New bands added downstream
must extend this for the multi-band `track` call to find them.
"""
function band_key end
@inline band_key(::GNSSSignals.L1) = :l1
@inline band_key(::GNSSSignals.L5) = :l5

"""
$(SIGNATURES)

Type alias for a NamedTuple of `Measurement`s — the multi-band input
shape of `track` / `track!`. Keys are band symbols (see [`band_key`](@ref)).
"""
const Measurements = NamedTuple{<:Any,<:Tuple{Vararg{Measurement}}}
