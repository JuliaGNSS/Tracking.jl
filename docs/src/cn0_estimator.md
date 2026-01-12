# CN0 Estimator

The CN0 (Carrier-to-Noise density ratio) estimator provides a measure of signal quality in dB-Hz.

## Default Estimator

The default CN0 estimator is the [Moment Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1) implemented as `MomentsCN0Estimator`. It buffers prompt correlation values and estimates CN0 based on the signal and noise power moments.

```@docs
MomentsCN0Estimator
```
You can adjust the buffer size for the MomentsCN0Estimator when creating a `SatState`:
```julia
sat_state = SatState(
    system, prn, code_phase, carrier_doppler,
    num_prompts_for_cn0_estimation = 100  # Buffer size for MomentsCN0Estimator
)
```

## Estimating CN0

To get the CN0 estimate from a tracking state:

```jldoctest cn0_example
julia> using Tracking, GNSSSignals

julia> using Tracking: Hz

julia> system = GPSL1();

julia> sat_state1 = SatState(system, 1, 50.0, 1000.0Hz);

julia> sat_state2 = SatState(system, 5, 120.0, -500.0Hz);

julia> track_state = TrackState(system, [sat_state1, sat_state2]);

julia> estimate_cn0(track_state, 1)  # Access by PRN
0.0 dB-Hz

julia> estimate_cn0(track_state, 5)  # Access by PRN
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

Pass your custom estimator when creating a `SatState`:

```jldoctest cn0_example
julia> prn = 1;

julia> code_phase = 50.0;

julia> carrier_doppler = 1000.0Hz;

julia> num_prompts = 100;  # Buffer size for MomentsCN0Estimator

julia> sat_state = SatState(
           system, prn, code_phase, carrier_doppler,
           num_prompts_for_cn0_estimation = num_prompts
       );
```

Pass your custom estimator to a `SatState`:
```julia
sat_state = SatState(
    sat_state;
    cn0_estimator = MyCN0Estimator()
)
```