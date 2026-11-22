"""
$(SIGNATURES)

Track the samples of one or more bands: downconvert, correlate and step the
tracking loops of every satellite. Returns the updated `TrackState`.

The first argument is one of:

| Argument                                       | Meaning                                                                                    |
|:---------------------------------------------- |:------------------------------------------------------------------------------------------ |
| `AbstractVecOrMat`                             | Bare sample buffer. Single-band TrackState only.                                           |
| [`BandSamples`](@ref)                          | One band's samples with their sample rate. Single-band TrackState.                         |
| [`SamplesByBand`](@ref Tracking.SamplesByBand) | Multi-band: a NamedTuple of one `BandSamples` per band id (see `GNSSSignals.get_band_id`). |

The bare buffer takes `sampling_frequency` as the next positional argument and
an `intermediate_frequency` keyword (default `0`). Every form runs the same
loop over its NCO update spans: downconvert and correlate all groups, then step the
tracking loops.

The returned `TrackState` is *structurally* detached from the input: each
group's key set and slot vector are copied, so adding, removing or tracking
satellites on one never affects the other. The copy is shallow, though, and
these are shared with the input (treat it as a stale handle, or `deepcopy` it
first):

  - the per-satellite scratch vectors (`filtered_prompts`, `correlator_bank_outputs`,
    the soft bits, the C/N₀ estimator's buffer), overwritten by the next call on
    either state; a bare `downconvert_and_correlate` likewise appends to the
    input's `correlator_bank_outputs`;
  - each signal's noise estimator, which `track` advances in place: branched
    states divide by one noise reference, and advancing two concurrently races,
    so build one `TrackState` per thread (see `downconvert_and_correlate` for why
    the window is not copied);
  - the `tracking_loop` itself, so results can be read off the instance you
    passed in. One with state of its own (TrackingLoops' `VectorTrackingLoop`) is
    advanced by `track` as by [`track!`](@ref), and two branches would mix their
    integrations in one navigation filter.

For real-time loops, construct the backend once outside the loop and pass it as
`downconvert_and_correlator` (see [`track!`](@ref)).

The coherent-integration length is set per signal with
[`set_preferred_num_code_periods_to_integrate!`](@ref), not passed to `track`.
"""
function track(
    samples_by_band::SamplesByBand,
    track_state::TS;
    kwargs...,
) where {TS<:TrackState}
    # Detach keys and values once (#123), then run in place (#133).
    detached = TrackState(
        track_state;
        groups = _detach_groups_slot_vectors(track_state.groups),
        sample_index_offsets = copy(track_state.sample_index_offsets),
    )
    track!(samples_by_band, detached; kwargs...)::TS
end

# Key a single band's samples by the state's only band (errors if multi-band).
@inline function _as_samples_by_band(band_samples::BandSamples, track_state::TrackState)
    key = get_band_id(_single_band(track_state))
    NamedTuple{(key,)}((band_samples,))
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track(
    samples::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    band_samples = BandSamples(samples, sampling_frequency; intermediate_frequency)
    track(_as_samples_by_band(band_samples, track_state), track_state; kwargs...)
end

# Single-band convenience.
function track(band_samples::BandSamples, track_state::TrackState; kwargs...)
    track(_as_samples_by_band(band_samples, track_state), track_state; kwargs...)
end

"""
$(SIGNATURES)

In-place version of [`track`](@ref): overwrites the satellite slots of
`track_state` and returns it.

After one warmup call the single-threaded path is allocation-free. The threaded
path (`CPUThreadedDownconvertAndCorrelator`) keeps a residual of about 96 B per
group per correlate pass (one per `nco_update_interval` span, so it scales
with the buffer length), when running on more than one thread with more than one
work item in the group (the noise despreads count, riding the first group); see
[`CPUThreadedDownconvertAndCorrelator`](@ref) for why.

For real-time loops, **construct the backend once outside the loop**:

```julia
track_state = TrackState(signal, initial_sats)
dc = CPUThreadedDownconvertAndCorrelator()      # hoist!

while got_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

The per-thread scratch buffers grow on first use; the default regrows them every
call.
"""
function track!(
    samples_by_band::SamplesByBand,
    track_state::TrackState;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    nco_update_interval = nothing,
)
    _validate_samples_by_band(track_state, samples_by_band)
    _prepare_for_new_buffer!(track_state)
    nco_update_interval = _resolve_nco_update_interval(nco_update_interval, track_state)
    _validate_nco_update_interval(nco_update_interval, samples_by_band)
    # Per NCO update span, one correlate pass at fixed NCO Doppler, then one estimate that
    # updates every satellite at a common epoch. The pass stops at each satellite's
    # last completed integration boundary (`stop_before_partial`), so the next pass
    # covers boundary to boundary at the updated Doppler: every integration sees
    # one Doppler, and each correction takes effect at its completing boundary.
    # The buffers are fixed for the call, so sample-derived caches (the bit
    # backends' band pack) are built on the first pass only.
    nco_update_index = 0
    while _nco_updates_left(nco_update_interval, nco_update_index, samples_by_band)
        downconvert_and_correlate!(
            downconvert_and_correlator,
            samples_by_band,
            track_state;
            nco_update_index,
            nco_update_interval,
            stop_before_partial = true,
            samples_unchanged = nco_update_index > 0,
        )
        step_tracking_loops!(track_state, samples_by_band)
        nco_update_index += 1
    end
    # Drain each satellite's trailing partial into its correlator bank, carried into
    # the next call; a boundary exactly at the buffer end completes here, hence
    # the loop update. Noise is measured here only if no span ran, or every sample
    # would enter the noise window twice.
    downconvert_and_correlate!(
        downconvert_and_correlator,
        samples_by_band,
        track_state;
        samples_unchanged = nco_update_index > 0,
        measure_noise = nco_update_index == 0,
    )
    step_tracking_loops!(track_state, samples_by_band)
    _advance_sample_index_offsets!(track_state.sample_index_offsets, samples_by_band)
    return track_state
end

# The next samples of every band start where these ended.
@inline _advance_sample_index_offsets!(
    sample_index_offsets::Dict{Symbol,Int},
    samples_by_band::SamplesByBand,
) = _advance_sample_index_offsets!(
    sample_index_offsets,
    keys(samples_by_band),
    values(samples_by_band),
)
@inline _advance_sample_index_offsets!(
    sample_index_offsets::Dict{Symbol,Int},
    ::Tuple{},
    ::Tuple{},
) = sample_index_offsets
@inline function _advance_sample_index_offsets!(
    sample_index_offsets::Dict{Symbol,Int},
    bands::Tuple,
    samples_by_band::Tuple,
)
    band = first(bands)
    sample_index_offsets[band] =
        get(sample_index_offsets, band, 0) + get_num_samples(first(samples_by_band))
    _advance_sample_index_offsets!(
        sample_index_offsets,
        Base.tail(bands),
        Base.tail(samples_by_band),
    )
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track!(
    samples::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    band_samples = BandSamples(samples, sampling_frequency; intermediate_frequency)
    track!(_as_samples_by_band(band_samples, track_state), track_state; kwargs...)
end

# Single-band convenience.
function track!(band_samples::BandSamples, track_state::TrackState; kwargs...)
    track!(_as_samples_by_band(band_samples, track_state), track_state; kwargs...)
end

# True until the previous NCO update span reached the buffer end on every band, independent
# of per-satellite progress. For `nco_update_index == 0`: is any buffer non-empty.
@inline function _nco_updates_left(
    nco_update_interval,
    nco_update_index::Int,
    samples_by_band,
)
    for band_samples in samples_by_band
        n = get_num_samples(band_samples)
        last(
            _nco_update_sample_range(
                nco_update_interval,
                nco_update_index - 1,
                band_samples.sampling_frequency,
                n,
            ),
        ) < n && return true
    end
    false
end

# `nothing` => the smallest primary-code period of any tracked signal.
@inline _resolve_nco_update_interval(nco_update_interval, ::TrackState) =
    nco_update_interval
@inline function _resolve_nco_update_interval(::Nothing, track_state::TrackState)
    _smallest_code_period(track_state)
end

# Primary-code period (a time) of one signal.
@inline _code_period(signal::AbstractGNSSSignal) =
    get_code_length(signal) / get_chip_rate(signal)

@inline _min_signal_code_period(::Tuple{}, m) = m
@inline _min_signal_code_period(signals::Tuple, m) =
    _min_signal_code_period(Base.tail(signals), min(m, _code_period(first(signals))))

@inline _min_group_code_period(::Tuple{}, m) = m
@inline _min_group_code_period(groups::Tuple, m) = _min_group_code_period(
    Base.tail(groups),
    _min_signal_code_period(first(groups).signals, m),
)

function _smallest_code_period(track_state::TrackState)
    groups = Tuple(track_state.groups)
    # Seeded from the first signal, which always exists, for type stability.
    init = _code_period(first(first(groups).signals))
    _min_group_code_period(groups, init)
end

# An NCO update span must cover a sample on every band, or `track!` would not terminate. A
# plain number gets a clear error rather than a Unitful one.
function _validate_nco_update_interval(nco_update_interval, samples_by_band::SamplesByBand)
    dimension(nco_update_interval) == dimension(1.0s) || throw(
        ArgumentError(
            "nco_update_interval must be a time quantity, e.g. `1u\"ms\"` or `1e-3u\"s\"` " *
            "(with `using Unitful`); got $nco_update_interval.",
        ),
    )
    for band_samples in samples_by_band
        samples_per_nco_update =
            uconvert(NoUnits, nco_update_interval * band_samples.sampling_frequency)
        samples_per_nco_update >= 1 || throw(
            ArgumentError(
                "nco_update_interval $nco_update_interval is shorter than one sample period " *
                "at sampling frequency $(band_samples.sampling_frequency); pick a longer interval.",
            ),
        )
    end
    nothing
end
