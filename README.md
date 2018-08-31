[![pipeline status](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
# Tracking
Tracks GNSS signals.

## Features

* DLL/PLL Discriminators
* Loop Filters of 1st, 2nd, and 3rd order, bilinear or boxcar
* DLL/PLL loop functions
* Tracking loop
* Joined tracking of multiple GNSS systems

## Getting started

Install:
```julia
Pkg.clone("git@git.rwth-aachen.de:nav/GNSSSignals.jl.git")
Pkg.clone("git@git.rwth-aachen.de:nav/Tracking.jl.git")
```

## Usage

```julia
    using Tracking
    using GNSSSignals
    function beamform(x)
        [0.5 0.5 0.5 0.5] * x
    end
    inits = Initial(50, pi / 3, 0.0, 200)
    prn = 1
    track = init_tracking(GPSL1(), inits, 1000, 4e6, 18, 1, prn)
    track(signal)
```

## Todo

Provide very early and very late when necessary.

## License

MIT License
