# CN0 Estimator

The default CN0 estimator is the [Moment Method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4621371&tag=1) called `MomentsCN0Estimator`.

You can easily add your own estimator. To do that you have to implement your own structure
`MyCN0Estimator <: AbstractCN0Estimator`
and the following functions: `update(cn0_estimator::MyCN0Estimator, prompt)` and
`estimate_cn0(cn0_estimator::MomentsCN0Estimator, integration_time)`.

Thereby, you can pass your CN0 estimator to the `TrackingState` by calling
```julia
TrackingState(GPSL1, carrier_doppler, code_phase, cn0_estimator = MyCN0Estimator())
```
