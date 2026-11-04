[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7454547.svg)](https://doi.org/10.5281/zenodo.7454547)
[![codecov](https://codecov.io/gh/JuliaGNSS/Tracking.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNSS/Tracking.jl)
[![[Semantic Release]](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--release-e10079.svg)](https://github.com/semantic-release/semantic-release)

# Tracking.jl

GNSS satellite tracking, as part of a GNSS receiver: signal downconversion and
correlation (SIMD-optimized CPU kernels), and the code and carrier tracking loops
(by default an FLL-assisted PLL and a DLL). Both can be replaced
through multiple dispatch.

The per-integration loop arithmetic (correlators, discriminators, bit
synchronization, C/N₀ estimators and tracking loops) lives in [TrackingLoops.jl](https://github.com/JuliaGNSS/TrackingLoops.jl),
so a hardware correlator's loop can run the same code. Tracking.jl re-exports the
functions that read results (`estimate_cn0`, `get_prompt`, `get_soft_bits`, …);
load TrackingLoops next to it to configure a correlator, a tracking loop, a C/N₀
estimator or a beamformer (`using Tracking, TrackingLoops, GNSSSignals`).

## Features

* Every concrete signal GNSSSignals.jl defines (GPS L1 C/A, L1C, L2C, L5; Galileo E1, E5a incl. E5a-QP, E5b, E6; BeiDou B1I, B3I, B1C, B2a, B2b); see the [signal support matrix](https://JuliaGNSS.github.io/Tracking.jl/stable/signals.html)
* Multi-band tracking: one `track` call processes the samples of several RF bands (e.g. L1 + L5) in lockstep
* Multi-signal tracking on one satellite (e.g. L1 C/A + L1C-D + L1C-P) with one shared downconvert
* Coherent integration over `N` code periods once the bit or secondary code is found
* CN0 estimation
* Secondary code detection
* Bit detection
* Phased array tracking
* Multi-satellite and multi-constellation tracking

## Installation

```julia
julia> ]
pkg> add Tracking
```

## Documentation

For usage examples and detailed documentation, see the [stable docs](https://JuliaGNSS.github.io/Tracking.jl/stable) or [dev docs](https://JuliaGNSS.github.io/Tracking.jl/dev).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
