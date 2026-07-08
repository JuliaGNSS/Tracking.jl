[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7454547.svg)](https://doi.org/10.5281/zenodo.7454547)
[![codecov](https://codecov.io/gh/JuliaGNSS/Tracking.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNSS/Tracking.jl)
[![[Semantic Release]](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--release-e10079.svg)](https://github.com/semantic-release/semantic-release)

# Tracking.jl

This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch).
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation using hand-tuned SIMD intrinsics.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs with FLL-assisted carrier tracking as the default.
Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 C/A from 0 to `N`*1023). The factor `N` can be specified, but will be `1` as long as the bit start is unknown in order to find the bit start. Once that is done for every tracked satellite the result will be handed over to the code and carrier estimation block.
Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.
Multi-signal tracking is supported: a single satellite can be tracked on several signals at once (e.g. GPS L1 C/A together with L1C-D and L1C-P) sharing one carrier downconvert per outer iteration, with a per-signal correlator each.
Multi-band tracking is supported as well: one `track` call can process sample buffers from several RF bands at once (e.g. L1 and L5), each band bundled as a `Measurement` with its own sampling frequency and intermediate frequency.

## Features

* Supports GPS L1 C/A, GPS L1C-D, GPS L1C-P, GPS L5I, and Galileo E1B (including the BOC(1,1) approximation `GalileoE1B_BOC11`)
* Multi-band tracking: one `track` call processes measurements from multiple RF bands (e.g. L1 + L5) in lockstep
* Multi-signal tracking on a single satellite (e.g. L1 C/A + L1C-D + L1C-P together) with a single shared downconvert
* CN0 estimation
* Secondary code detection
* Bit detection
* Phased array tracking
* Multi-satellite and multi-system tracking

## Installation

```julia
julia> ]
pkg> add Tracking
```

## Documentation

For usage examples and detailed documentation, see the [stable docs](https://JuliaGNSS.github.io/Tracking.jl/stable) or [dev docs](https://JuliaGNSS.github.io/Tracking.jl/dev).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

<!-- ci probe (temporary): does base pass Julia 1.10 with today's deps? -->
