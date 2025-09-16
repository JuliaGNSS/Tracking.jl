[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaGNSS.github.io/Tracking.jl/dev)
[![DOI](https://zenodo.org/badge/171484097.svg)](https://zenodo.org/badge/latestdoi/171484097)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNSS/Tracking.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaGNSS/Tracking.jl?branch=master)
[![[Semantic Release]](https://img.shields.io/badge/%20%20%F0%9F%93%A6%F0%9F%9A%80-semantic--release-e10079.svg)](https://github.com/semantic-release/semantic-release)

# Tracking
This package implements the tracking functionality of GNSS satellites that's part of the larger GNSS receiver.
Tracking.jl primarily consists of two main blocks:

1. Signal down-conversion and correlation
2. Code and carrier estimation to generate replicas and close the loop

Tracking.jl provides defaults for both blocks, but it provides mechanisms to hook in your own implementation (by using multiple dispatch). 
For signal down-conversion and correlation Tracking.jl provides a highly optimized CPU implementation. There is also a GPU implementation that is mainly there for reference, but is not yet as fast as the CPU implementation.
With respect to the second block Tracking.jl provides conventional DLLs and PLLs.
Down-conversion and correlation is done in full code blocks meaning from code start to code end or multiples of that (e.g. in GPS L1 from 0 to `N`*1023). The factor `N` can be specified, but will be `1` as long as the bit start is unknown in order to find the bit start. Once that is done for every tracked satellite the result will be handed over to the code and carrier estimation block.
Moreover, Tracking.jl allows tracking of signals from phased antenna arrays meaning that they are down-converted and correlated by the very same replica to conserve phase relationships.

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