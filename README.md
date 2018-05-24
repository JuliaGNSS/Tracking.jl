[![pipeline status](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/Tracking.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/Tracking.jl/commits/master)
# Tracking
Tracks GNSS signals. Currently it only provides DLL and PLL discriminators.

## Features

* DLL/PLL Discriminators
* Loop Filters of 1st, 2nd, and 3rd order

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
loop_filter_function = init_3rd_order_loop_filter(bandwidth, Δt)
```

## Todo

This is still missing:
* PLL/DLL
* Tracking Loop
* Add carrier aiding
* external velocity aiding

## Nice to have

* Multi Signal (satellite) support

## License

MIT License
