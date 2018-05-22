[![pipeline status](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
# Tracking
Tracks GNSS signals. Currently it only provides dll and pll discriminators.

## Features

* DLL/PLL Discriminators

## Getting started

Install:
```julia
Pkg.clone("git@git.rwth-aachen.de:nav/Tracking.jl.git")
```

## Usage

```julia
using Tracking
chip_error = dll_disc(replica_code_phases)
phase_error = pll_disc(replica_code_phases)


## Todo

This is still missing:

* Loop Filter
* PLL/DLL
* Tracking Loop
* Add carrier aiding
* external velocity aiding

## Nice to have

* Multi Signal (satallite) support

## License

MIT License
