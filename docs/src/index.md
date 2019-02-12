# Tracking.jl

Modular tracking algorithm for various Global Navigation Satellite Systems (GNSS) (e.g. GPS L1, GPS L5, ...)

## Introduction

This package provides a basic tracking function for GNSS. It aims to be modular and performant. Various loop filters are implemented as well as an estimator for the Carrier-to-Noise-Density-Ratio (CN0).

## Installation

```julia-repl
pkg> add Tracking
```

## Usage

```julia
using Tracking, GNSSSignals
import Unitful: MHz, Hz
gpsl1 = GPSL1()
carrier_doppler = 100Hz
code_phase = 120
inits = Initials(gpsl1, carrier_doppler, code_phase)
sample_freq = 2.5MHz
interm_freq = 0Hz
prn = 1
track = init_tracking(gpsl1, inits, sample_freq, interm_freq, prn)
track, track_results = track(signal)
```
