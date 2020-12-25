[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/dev)
[![Build Status](https://travis-ci.org/JuliaGNSS/Tracking.jl.svg?branch=master)](https://travis-ci.org/JuliaGNSS/Tracking.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/Tracking.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/Tracking.jl?branch=master)

# Tracking
This implements a basic tracking functionality for GNSS signals. The correlation is done in the interval of PRNs. Each call of the tracking function returns the current code phase, Doppler, the Carrier-to-Noise-Density-Ratio (CN0), data bits, number of data bits and the last correlator output.

## Features

* Supports GPS L1 / L5 and Galileo E1B
* CN0 estimation
* Secondary code detection
* Bit detection
* Phased array tracking

## Getting started

Install:
```julia
julia> ]
pkg> add Tracking
```

## Usage

```julia
using Tracking
using Tracking: Hz, GPSL1
carrier_doppler = 1000Hz
code_phase = 50
sampling_frequency = 2.5e6Hz
prn = 1
state = TrackingState(GPSL1, carrier_doppler, code_phase)
results = track(signal, state, prn, sampling_frequency)
next_results = track(next_signal, get_state(results), prn, sampling_frequency)
```

If you'd like to track several signals at once (e.g. in the case of phased antenna arrays), you will have to specify the optional parameter `num_ants::NumAnts{N}` and pass a beamforming function to the `track` function:

```julia
state = TrackingState(GPSL1, carrier_doppler, code_phase, num_ants = NumAnts(4)) # 4 antenna channels
results = track(signal, state, prn, sampling_frequency, post_corr_filter = x -> x[1]) # Post corr filter is optional
```
