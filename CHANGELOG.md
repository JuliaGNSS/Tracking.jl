# Changelog

# [0.16.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.15.8...v0.16.0) (2025-11-19)


### Bug Fixes

* add Pkg.develop in CI to properly resolve Tracking dependencies ([6a7c299](https://github.com/JuliaGNSS/Tracking.jl/commit/6a7c299b17ae792bed86ed749d2beeadcc678bf4))
* export types and functions needed by CUDA extension ([dde6f40](https://github.com/JuliaGNSS/Tracking.jl/commit/dde6f4070a6016b19bd44c044b9b0253d1c634d8))
* import functions from Tracking to properly extend in GPU extension ([fbd09aa](https://github.com/JuliaGNSS/Tracking.jl/commit/fbd09aabeec55aa5839be9880183abe7fb6aecde))
* import get_num_samples and update from Tracking in GPU extension ([90865cc](https://github.com/JuliaGNSS/Tracking.jl/commit/90865cccef2fa15ffdf9ebdb518f23e81558d038))
* properly export GPU types from extension to parent module ([a50656d](https://github.com/JuliaGNSS/Tracking.jl/commit/a50656d4a6e92d6591452574f3fa348098e5cd9a))
* remove __init__() pattern and update tests to import from TrackingCUDAExt ([49d0c44](https://github.com/JuliaGNSS/Tracking.jl/commit/49d0c446cdc49d6c8790a69a8dbf2fd6ad117d3e))
* use Base.get_extension to access GPU types in extension tests ([3e6ccd7](https://github.com/JuliaGNSS/Tracking.jl/commit/3e6ccd74734c3b06d585291a66e047c924caecbd))


### Features

* add TrackingCUDAExt extension entry point ([37f3113](https://github.com/JuliaGNSS/Tracking.jl/commit/37f3113764c24a586309ed7acdf94f088f192714))

## [0.15.8](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.15.7...v0.15.8) (2025-10-04)


### Bug Fixes

* use correct earliest and latest sample shift in gen_code_replica function ([a0be2d7](https://github.com/JuliaGNSS/Tracking.jl/commit/a0be2d7c761e124d9b49e599070f7e636a24f43b))
* use correct latest sample shift in correlate function ([8b8b49b](https://github.com/JuliaGNSS/Tracking.jl/commit/8b8b49b9c402ae5409ab43bbd4231515d886d457))
