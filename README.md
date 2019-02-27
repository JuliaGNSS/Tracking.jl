[![Build Status](https://travis-ci.org/JuliaGNSS/Tracking.jl.svg?branch=master)](https://travis-ci.org/JuliaGNSS/Tracking.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/Tracking.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/Tracking.jl?branch=master)

# Tracking
This implements a basic tracking functionality for GNSS signals. The correlation is done in the interval of PRNs. Each call of the tracking function returns the current code phase, doppler, the Carrier-to-Noise-Density-Ratio (CN0), data bits, number of data bits and the last valid correlator output.

## Features

* Supports Loop Filters of 1st, 2nd, and 3rd order, bilinear or boxcar
* Supports GPS L1 / L5
* CN0 estimation

## Getting started

Install:
```julia
] add https://github.com/JuliaGNSS/GNSSSignals.jl.git
] add https://github.com/JuliaGNSS/Tracking.jl.git
```

## Usage

```julia
using Tracking
import Tracking: MHz, Hz
gpsl1 = Tracking.GPSL1()
carrier_doppler = 100Hz
code_phase = 120
inits = TrackingInitials(gpsl1, carrier_doppler, code_phase)
sample_freq = 2.5MHz
interm_freq = 0Hz
prn = 1
track = init_tracking(gpsl1, inits, sample_freq, interm_freq, prn)
track, track_results = track(signal)
```

## Todo

* Support Galileo Signals

## License

MIT License
