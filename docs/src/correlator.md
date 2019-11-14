# Correlator

The default correlator is given by
`get_default_correlator(::Type{AbstractGNSSSystem}, num_ants::NumAnts{N})`.
That is for example the `EarlyPromptLateCorrelator` for `GPSL1`.

You can easily add your own correlator. For that you have define your own
`MyCorrelator <: AbstractCorrelator` and implement various functions, that you find in
`src/correlator.jl`
