"""
$(SIGNATURES)

Track one or more `BandMeasurement`s: downconvert, correlate and estimate the
Dopplers of every satellite, and return the updated `TrackState` with new
phase/Doppler estimates and decoded bits.

The first argument is one of:

| Argument                                | Meaning                                                                        |
|:--------------------------------------- |:------------------------------------------------------------------------------ |
| `AbstractVecOrMat`                      | Bare sample buffer. Single-band TrackState only.                               |
| `BandMeasurement`                       | One band's bundled buffer + sample rate. Single-band TrackState.               |
| `NamedTuple{...}` of `BandMeasurement`s | Multi-band: one `BandMeasurement` per band id (see `GNSSSignals.get_band_id`). |

The bare buffer takes `sampling_frequency` as the next positional argument and
an `intermediate_frequency` keyword (default `0`). Every form runs the same
per-chunk loop: downconvert and correlate all groups, then estimate.

The returned `TrackState` is *structurally* detached from the input: each
group's key set and slot vector are copied, so adding, removing or tracking
satellites on one never affects the other. The copy is shallow, though, and
these are shared with the input (treat it as a stale handle, or `deepcopy` it
first):

  - the per-satellite scratch vectors (`filtered_prompts`, `correlator_outputs`,
    the soft bits, the C/N₀ estimator's buffer), overwritten by the next call on
    either state; a bare `downconvert_and_correlate` likewise appends to the
    input's `correlator_outputs`;
  - each signal's noise estimator, which `track` advances in place: branched
    states divide by one noise reference, and advancing two concurrently races,
    so build one `TrackState` per thread (see `downconvert_and_correlate` for why
    the window is not copied);
  - the `doppler_estimator` itself, so results can be read off the instance you
    passed in. One with state of its own (TrackingLoops' `VectorPLLAndDLL`) is
    advanced by `track` as by [`track!`](@ref), and two branches would mix their
    records in one navigation filter.

For real-time loops, **construct the correlator once outside the loop**:

```julia
dc = CPUThreadedDownconvertAndCorrelator()
while got_chunk(rx)
    chunk = read_chunk!(rx)
    track_state =
        track(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

The default builds fresh per-thread scratch buffers on every call. See
[`track!`](@ref) for the in-place variant.

The coherent-integration length is set **per signal** with
[`set_preferred_num_code_blocks_to_integrate!`](@ref), not passed to `track`:

```julia
set_preferred_num_code_blocks_to_integrate!(track_state, :gps_l5, 1, GPSL5I, 10)  # PRN 1 L5I: 10 ms
```

The conventional estimator caps the carrier loop bandwidth at `0.09 / T_int` and
the code loop's at `0.018 / T_int`, so longer integration stays stable (see
[`ConventionalPLLAndDLL`](@extref TrackingLoops.ConventionalPLLAndDLL)).
"""
function track(
    measurements::BandMeasurements,
    track_state::TS;
    kwargs...,
) where {TS<:TrackState}
    # Detach keys and values once (#123), then run in place (#133).
    detached = TrackState(
        track_state;
        groups = _detach_groups_slot_vectors(track_state.groups),
        sample_offsets = copy(track_state.sample_offsets),
    )
    track!(measurements, detached; kwargs...)::TS
end

# Key a single measurement by the state's only band (errors if multi-band).
@inline function _single_band_measurements(
    measurement::BandMeasurement,
    track_state::TrackState,
)
    key = get_band_id(_single_band(track_state))
    NamedTuple{(key,)}((measurement,))
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    m = BandMeasurement(signal, sampling_frequency, intermediate_frequency)
    track(_single_band_measurements(m, track_state), track_state; kwargs...)
end

# Single-BandMeasurement convenience: still a single-band path.
function track(measurement::BandMeasurement, track_state::TrackState; kwargs...)
    track(_single_band_measurements(measurement, track_state), track_state; kwargs...)
end

"""
$(SIGNATURES)

In-place version of [`track`](@ref): overwrites the satellite slots of
`track_state` and returns it.

After one warmup call the single-threaded path is allocation-free. The threaded
path (`CPUThreadedDownconvertAndCorrelator`) keeps a residual of about 96 B per
group per correlate pass (one per `doppler_update_interval` chunk, so it scales
with the buffer length), when running on more than one thread with more than one
work item in the group (the noise despreads count, riding the first group); see
[`CPUThreadedDownconvertAndCorrelator`](@ref) for why.

For real-time loops, **construct the correlator once outside the loop**:

```julia
track_state = TrackState(signal, initial_sats)
dc = CPUThreadedDownconvertAndCorrelator()      # hoist!

while got_chunk(rx)
    chunk = read_chunk!(rx)
    track!(chunk, track_state, sampling_frequency; downconvert_and_correlator = dc)
end
```

Its per-thread scratch buffers grow on first use; the default would regrow them
every call.
"""
function track!(
    measurements::BandMeasurements,
    track_state::TrackState;
    downconvert_and_correlator::AbstractDownconvertAndCorrelator = CPUThreadedDownconvertAndCorrelator(),
    doppler_update_interval = nothing,
)
    _validate_measurements(track_state, measurements)
    reset_start_sample_and_bit_buffer!(track_state)
    chunk_duration = _resolve_doppler_update_interval(doppler_update_interval, track_state)
    _validate_doppler_update_interval(chunk_duration, measurements)
    # Per chunk, one correlate pass at fixed NCO Doppler, then one estimate that
    # updates every satellite at a common epoch. The pass stops at each
    # satellite's last completed code-block boundary (`stop_before_partial`), so
    # the next pass covers boundary to boundary at the updated Doppler: every
    # integration sees one Doppler, and each correction takes effect at its
    # completing boundary.
    #
    # The buffers are fixed for the call, so sample-derived caches (the bit
    # backends' band pack) are built on the first pass only.
    chunk_index = 0
    while _chunks_left(chunk_duration, chunk_index, measurements)
        downconvert_and_correlate!(
            downconvert_and_correlator,
            measurements,
            track_state;
            chunk_index,
            chunk_duration,
            stop_before_partial = true,
            samples_unchanged = chunk_index > 0,
        )
        estimate_dopplers_and_filter_prompt!(track_state, measurements)
        chunk_index += 1
    end
    # Drain each satellite's trailing partial into its accumulator, carried into
    # the next call; a boundary exactly at the buffer end completes here, hence
    # the fold. Noise is measured here only if no chunk ran, or every sample
    # would enter the noise window twice.
    downconvert_and_correlate!(
        downconvert_and_correlator,
        measurements,
        track_state;
        samples_unchanged = chunk_index > 0,
        measure_noise = chunk_index == 0,
    )
    estimate_dopplers_and_filter_prompt!(track_state, measurements)
    _advance_sample_offsets!(track_state.sample_offsets, measurements)
    return track_state
end

# The next measurement of every band starts where this one ended.
@inline _advance_sample_offsets!(
    sample_offsets::Dict{Symbol,Int},
    measurements::BandMeasurements,
) = _advance_sample_offsets!(sample_offsets, keys(measurements), values(measurements))
@inline _advance_sample_offsets!(sample_offsets::Dict{Symbol,Int}, ::Tuple{}, ::Tuple{}) =
    sample_offsets
@inline function _advance_sample_offsets!(
    sample_offsets::Dict{Symbol,Int},
    bands::Tuple,
    measurements::Tuple,
)
    band = first(bands)
    sample_offsets[band] =
        get(sample_offsets, band, 0) + get_num_samples(first(measurements))
    _advance_sample_offsets!(sample_offsets, Base.tail(bands), Base.tail(measurements))
end

# Bare-buffer convenience wrapper. Single-band TrackStates only.
function track!(
    signal::AbstractVecOrMat,
    track_state::TrackState,
    sampling_frequency;
    intermediate_frequency = zero(sampling_frequency),
    kwargs...,
)
    m = BandMeasurement(signal, sampling_frequency, intermediate_frequency)
    track!(_single_band_measurements(m, track_state), track_state; kwargs...)
end

# Single-BandMeasurement convenience: still a single-band path.
function track!(measurement::BandMeasurement, track_state::TrackState; kwargs...)
    track!(_single_band_measurements(measurement, track_state), track_state; kwargs...)
end

# True until the previous chunk reached the buffer end on every band, independent
# of per-satellite progress. For `chunk_index == 0`: is any buffer non-empty.
@inline function _chunks_left(chunk_duration, chunk_index::Int, measurements)
    for m in measurements
        n = get_num_samples(m)
        _chunk_last_sample(chunk_duration, chunk_index - 1, m.sampling_frequency, n) < n &&
            return true
    end
    false
end

# `nothing` => the smallest primary-code period of any tracked signal.
@inline _resolve_doppler_update_interval(doppler_update_interval, ::TrackState) =
    doppler_update_interval
@inline function _resolve_doppler_update_interval(::Nothing, track_state::TrackState)
    _smallest_code_period(track_state)
end

# Primary-code period (a time) of one signal.
@inline _code_period(signal::AbstractGNSSSignal) =
    get_code_length(signal) / get_code_frequency(signal)

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

# A chunk must cover a sample on every band, or `track!` would not terminate. A
# plain number gets a clear error rather than a Unitful one.
function _validate_doppler_update_interval(chunk_duration, measurements::BandMeasurements)
    dimension(chunk_duration) == dimension(1.0s) || throw(
        ArgumentError(
            "doppler_update_interval must be a time quantity, e.g. `1u\"ms\"` or `1e-3u\"s\"` " *
            "(with `using Unitful`); got $chunk_duration.",
        ),
    )
    for m in measurements
        samples_per_chunk = uconvert(NoUnits, chunk_duration * m.sampling_frequency)
        samples_per_chunk >= 1 || throw(
            ArgumentError(
                "doppler_update_interval $chunk_duration is shorter than one sample period " *
                "at sampling frequency $(m.sampling_frequency); pick a longer interval.",
            ),
        )
    end
    nothing
end
