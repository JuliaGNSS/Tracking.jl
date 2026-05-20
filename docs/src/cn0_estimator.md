# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of signal quality in dB-Hz. Each [`TrackedSignal`](@ref) on a [`TrackedSat`](@ref) holds its own CN0 estimator, so a multi-signal satellite produces one CN0 value per signal.

## Default Estimator

The default CN0 estimator is the [Moment Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1) implemented as `MomentsCN0Estimator`. It buffers prompt correlation values and estimates CN0 based on the signal and noise power moments.

```@docs
MomentsCN0Estimator
```

You can adjust the buffer size for the `MomentsCN0Estimator` by constructing the `TrackedSat` yourself and handing it to `add_satellite!`'s escape-hatch overload:

```jldoctest cn0_buffer
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> sat = TrackedSat(GPSL1CA(), 1, 50.0, 1000.0Hz;
                        num_prompts_for_cn0_estimation = 200,
                        doppler_estimator = ConventionalAssistedPLLAndDLL());

julia> add_satellite!(track_state, :default, sat);

julia> get_prn(track_state, :default, 1)
1
```

## Estimating CN0

To get the CN0 estimate from a tracking state:

```jldoctest cn0_example
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> track_state = TrackState(; signal = GPSL1CA());

julia> add_satellite!(track_state; prn = 1, code_phase = 50.0,  carrier_doppler = 1000.0Hz);

julia> add_satellite!(track_state; prn = 5, code_phase = 120.0, carrier_doppler = -500.0Hz);

julia> estimate_cn0(track_state, 1)  # access by PRN
0.0 dB-Hz

julia> estimate_cn0(track_state, 5)
0.0 dB-Hz
```

Note: CN0 will be 0.0 dB-Hz until actual signal data has been processed through the tracking loop.

```@docs
estimate_cn0
```

## Custom CN0 Estimators

You can implement your own estimator by creating a subtype of `AbstractCN0Estimator`
and implementing the following functions:

- `update(cn0_estimator::MyCN0Estimator, prompt)`: Update the estimator with a new prompt value
- `estimate_cn0(cn0_estimator::MyCN0Estimator, integration_time)`: Return the CN0 estimate

To plug your custom estimator into a satellite, build a [`TrackedSignal`](@ref) with the custom `cn0_estimator` (via the kwarg-update constructor on an existing `TrackedSignal`), splice it into a fresh [`TrackedSat`](@ref), and add via the escape-hatch overload — see the [`TrackedSat`](@ref) docstring for the full constructor surface.
