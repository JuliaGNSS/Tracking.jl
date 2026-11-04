"""
    BandSamples(samples, sampling_frequency; intermediate_frequency = 0Hz)

One band's digitized samples with the `sampling_frequency` and
`intermediate_frequency` they were captured at: complex baseband samples at the
default intermediate frequency of zero, low-IF samples otherwise.

Fields:

  - `samples::S`: complex sample buffer (`Vector` for one antenna, `Matrix`
    with rows = samples and columns = antennas for an antenna array).
    Must be dense (unit row stride, columns packed back-to-back), as the SIMD
    kernels read it through raw pointers; anything else is an `ArgumentError`.
    Contiguous `view`s (e.g. `view(buf, 1:4000)`) are fine.
  - `sampling_frequency::F`: the buffer's sample rate (e.g. `4e6Hz`)
  - `intermediate_frequency::F`: the band's intermediate frequency (`0.0Hz` by
    default)

A multi-band `track` takes a NamedTuple of them keyed by band
([`MultiBandSamples`](@ref Tracking.MultiBandSamples)); a single band may also
pass a plain buffer and sample rate.

```julia
BandSamples(buf, 4e6Hz)                                     # baseband, IF 0.0Hz
BandSamples(buf, 4e6Hz; intermediate_frequency = 1.25e6Hz)  # low IF
```
"""
struct BandSamples{S<:AbstractVecOrMat,F}
    samples::S
    sampling_frequency::F
    intermediate_frequency::F

    function BandSamples{S,F}(
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
                "BandSamples `samples` must be densely laid out in memory ",
                "(unit row stride; for a matrix, columns packed back-to-back): ",
                "the SIMD kernels read it through raw pointers. Got ",
                typeof(samples),
                ". Pass a dense `Vector`/`Matrix` or a contiguous `view` of one.",
            ),
        ),
    )
end

# Both frequencies promoted to one type `F`.
function BandSamples(
    samples::AbstractVecOrMat,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
)
    sampling_frequency, intermediate_frequency =
        promote(sampling_frequency, intermediate_frequency)
    BandSamples{typeof(samples),typeof(sampling_frequency)}(
        samples,
        sampling_frequency,
        intermediate_frequency,
    )
end

@inline get_num_samples(band_samples::BandSamples) = get_num_samples(band_samples.samples)

"""
$(SIGNATURES)

A NamedTuple of [`BandSamples`](@ref), one per band, the multi-band input of
`track` / `track!`,
keyed by `GNSSSignals.get_band_id` (e.g. `:L1`, `:L5`), a compile-time constant.
"""
const MultiBandSamples = NamedTuple{<:Any,<:Tuple{Vararg{BandSamples}}}
