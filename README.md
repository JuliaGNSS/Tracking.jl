[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/dev)
[![DOI](https://zenodo.org/badge/171484097.svg)](https://zenodo.org/badge/latestdoi/171484097)
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
* GPU acceleration (CUDA)

## Getting started

Install:
```julia
julia> ]
pkg> add Tracking
```

## Usage

```julia
using GNSSSignals
using Tracking
using Tracking: Hz
carrier_doppler = 1000Hz
code_phase = 50
sampling_frequency = 2.5e6Hz
prn = 1
gpsl1 = GPSL1()
state = TrackingState(prn, gpsl1, carrier_doppler, code_phase)
results = track(signal, state, sampling_frequency)
next_results = track(next_signal, get_state(results), sampling_frequency)
```

If you'd like to track several signals at once (e.g. in the case of phased antenna arrays), you will have to specify the optional parameter `num_ants::NumAnts{N}` and pass a beamforming function to the `track` function:

```julia
state = TrackingState(prn, gpsl1, carrier_doppler, code_phase, num_ants = NumAnts(4)) # 4 antenna channels
results = track(signal, state, sampling_frequency, post_corr_filter = x -> x[1]) # Post corr filter is optional
```

### Usage with `CUDA.jl`
This package supports accelerating the tracking loop by using the GPU. At the moment support is only provided for `CUDA.jl`. If you'd like to use this option, you'd have to opt-in by providing the following argument upon creating an `AbstractGNSS`:
``` julia
gpsl1_gpu = GPSL1(use_gpu = Val(true))
```
Beware that `num_samples` must be provided explicitly upon creating a `TrackingState`:
``` julia
state_gpu = TrackingState(prn, gpsl1_gpu, carrier_doppler, code_phase, num_samples = N)
```
Moreover, your signal must be a `StructArray{ComplexF32}` of `CuArray{Float32}` type:
``` julia
using StructArrays
signal_cu = CuArray{ComplexF32}(signal_cpu)
signal_gpu = StructArray(signal_cu)
```
Otherwise the usage is identical to the example provided above, including the case for multi-antenna tracking:
``` julia
results_gpu = track(signal_gpu, state_gpu, sampling_frequency)
next_results_gpu = track(next_signal_gpu, get_state(results_gpu), sampling_frequency)
```