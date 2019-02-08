[![pipeline status](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
# Tracking
This implements a basic tracking function of GNSS signals. The correlation is done in the interval of PRNs. Each call of the tracking function returns the current code phase, doppler, data bits, number of data bits and the last valid correlator output.

## Features

* Loop Filters of 1st, 2nd, and 3rd order, bilinear or boxcar
* GPSL1
* GPSL5

## Getting started

Install:
```julia
] add git@git.rwth-aachen.de:nav/Tracking.jl.git
```
If you have not added `GNSSSignals` before, you will also need to download `GNSSSignals`
```julia
] add git@git.rwth-aachen.de:nav/GNSSSignals.jl.git
```

## Usage

```julia
using Tracking, GNSSSignals
import Unitful: Hz
gpsl1 = GPSL1()
carrier_doppler = 100Hz
code_phase = 120
inits = Initials(gpsl1, carrier_doppler, code_phase)
sample_freq = 2.5e6Hz
interm_freq = 0Hz
prn = 1
track = init_tracking(gpsl1, inits, sample_freq, interm_freq, prn)
track, track_results = track(signal)
```

## Todo

* Galileo Signals

## License

MIT License
