"""
$(SIGNATURES)

One band's sample buffer with the `sampling_frequency` and
`intermediate_frequency` it was captured at.

Fields:

  - `samples::S`: complex sample buffer (`Vector` for one antenna, `Matrix`
    with rows = samples and columns = antennas for an antenna array).
    Must be dense (unit row stride, columns packed back-to-back), as the SIMD
    kernels read it through raw pointers; anything else is an `ArgumentError`.
    Contiguous `view`s (e.g. `view(buf, 1:4000)`) are fine.
  - `sampling_frequency::F`: the buffer's sample rate (e.g. `4e6Hz`)
  - `intermediate_frequency::F`: the band's IF (defaults to `0.0Hz`)

A multi-band `track` takes a NamedTuple of them keyed by band; a single band
may also pass a plain buffer and sample rate.

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

# Dense arrays and contiguous strided views pass; anything else (a `1:2:end`
# view, a transpose, …) would silently correlate the wrong samples through the
# kernels' pointer arithmetic.
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

# Promote both frequencies to one type `F`.
function BandMeasurement(
    samples::AbstractVecOrMat,
    sampling_frequency,
    intermediate_frequency,
)
    BandMeasurement(samples, promote(sampling_frequency, intermediate_frequency)...)
end

function BandMeasurement(samples::AbstractVecOrMat, sampling_frequency)
    intermediate_frequency = zero(sampling_frequency)
    BandMeasurement(samples, sampling_frequency, intermediate_frequency)
end

function BandMeasurement(
    samples::AbstractVecOrMat;
    sampling_frequency,
    intermediate_frequency = zero(sampling_frequency),
)
    BandMeasurement(samples, sampling_frequency, intermediate_frequency)
end

@inline get_num_samples(m::BandMeasurement) = get_num_samples(m.samples)

"""
$(SIGNATURES)

A NamedTuple of `BandMeasurement`s, the multi-band input of `track` / `track!`,
keyed by `GNSSSignals.get_band_id` (e.g. `:L1`, `:L5`), a compile-time constant.
"""
const BandMeasurements = NamedTuple{<:Any,<:Tuple{Vararg{BandMeasurement}}}
