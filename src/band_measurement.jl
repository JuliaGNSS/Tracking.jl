"""
$(SIGNATURES)

One band's incoming sample buffer plus the front-end metadata needed to
process it. Bundles `samples` with the `sampling_frequency` and
`intermediate_frequency` they were captured at — these are inseparable
in practice, and the bundle removes the chance of mismatched parallel
NamedTuples in a multi-band `track` call.

Fields:

  - `samples::S`: complex sample buffer (`Vector` for one antenna, `Matrix`
    with rows = samples and columns = antennas for an antenna array).
    Must be densely laid out in memory (unit row stride, columns packed
    back-to-back) — the SIMD downconvert/correlate kernels read the buffer
    through raw pointers with dense column-stride math, so a non-contiguous
    strided view would silently correlate the wrong samples. The constructor
    validates this and rejects non-dense buffers with an `ArgumentError`;
    contiguous `view`s (e.g. `view(buf, 1:4000)`) remain fine.
  - `sampling_frequency::F`: the buffer's sample rate (e.g. `4e6Hz`)
  - `intermediate_frequency::F`: the band's IF (defaults to `0.0Hz`)

In a multi-band call, one `BandMeasurement` is built per band; a NamedTuple
of `BandMeasurement`s keyed by band feeds `track`. For the single-band case
a plain buffer + scalar sample-rate keeps working unchanged.

```julia
BandMeasurement(buf, 4e6Hz)                              # IF defaults to 0.0Hz
BandMeasurement(buf, 4e6Hz, 1.575e6Hz)                   # explicit IF
BandMeasurement(buf; sampling_frequency = 4e6Hz)         # kwarg form
```
"""
struct BandMeasurement{S<:AbstractVecOrMat,F}
    samples::S
    sampling_frequency::F
    intermediate_frequency::F

    function BandMeasurement(
        samples::S,
        sampling_frequency::F,
        intermediate_frequency::F,
    ) where {S<:AbstractVecOrMat,F}
        _assert_dense_layout(samples)
        new{S,F}(samples, sampling_frequency, intermediate_frequency)
    end
end

# The downconvert/correlate kernels read `samples` via `pointer` with
# dense column-stride math, so the buffer must be unit-strided with
# packed columns. Dense arrays trivially qualify; strided views qualify
# only when contiguous. Anything else (a `1:2:end` view, a transpose, a
# row-selected matrix view, …) would silently correlate the wrong data —
# reject it here at construction.
@inline _assert_dense_layout(::DenseVecOrMat) = nothing
function _assert_dense_layout(samples::AbstractVecOrMat)
    if samples isa StridedVecOrMat &&
       stride(samples, 1) == 1 &&
       (samples isa AbstractVector || stride(samples, 2) == size(samples, 1))
        return nothing
    end
    throw(
        ArgumentError(
            string(
                "BandMeasurement `samples` must be densely laid out in memory ",
                "(unit row stride; for a matrix, columns packed back-to-back): ",
                "the SIMD kernels read it through raw pointers. Got ",
                typeof(samples),
                ". Pass a dense `Vector`/`Matrix` or a contiguous `view` of one.",
            ),
        ),
    )
end

# Promotes the numeric type so `sampling_frequency` and
# `intermediate_frequency` end up the same concrete type (`F`), e.g. an
# integer-typed IF alongside a `Float64` sampling frequency.
function BandMeasurement(
    samples::AbstractVecOrMat,
    sampling_frequency,
    intermediate_frequency,
)
    BandMeasurement(samples, promote(sampling_frequency, intermediate_frequency)...)
end

# Positional constructor with IF defaulting to 0.0Hz.
function BandMeasurement(samples::AbstractVecOrMat, sampling_frequency)
    intermediate_frequency = zero(sampling_frequency)
    BandMeasurement(samples, sampling_frequency, intermediate_frequency)
end

# Kwarg constructor. `intermediate_frequency` defaults to zero of the
# same type as `sampling_frequency` so the struct's `F` stays unified.
function BandMeasurement(
    samples::AbstractVecOrMat;
    sampling_frequency,
    intermediate_frequency = zero(sampling_frequency),
)
    BandMeasurement(samples, sampling_frequency, intermediate_frequency)
end

@inline get_samples(m::BandMeasurement) = m.samples
@inline get_sampling_frequency(m::BandMeasurement) = m.sampling_frequency
@inline get_intermediate_frequency(m::BandMeasurement) = m.intermediate_frequency
@inline get_num_samples(m::BandMeasurement) = get_num_samples(m.samples)

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
@inline band_key(::GNSSSignals.L2) = :l2
@inline band_key(::GNSSSignals.L5) = :l5
# Fallback for bands without a `band_key` method: a raw `MethodError` deep
# inside the multi-band machinery would be opaque, so point the caller at
# the fix.
band_key(band::GNSSSignals.Band) = throw(
    ArgumentError(
        string(
            "No `band_key` method for band `",
            band,
            "` (",
            typeof(band),
            "). ",
            "Define `Tracking.band_key(::",
            typeof(band),
            ") = :some_symbol` so the ",
            "multi-band `track` call can key its measurements by this band.",
        ),
    ),
)

"""
$(SIGNATURES)

Type alias for a NamedTuple of `BandMeasurement`s — the multi-band input
shape of `track` / `track!`. Keys are band symbols (see [`band_key`](@ref)).
"""
const BandMeasurements = NamedTuple{<:Any,<:Tuple{Vararg{BandMeasurement}}}
