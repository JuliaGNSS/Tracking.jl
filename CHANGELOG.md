# Changelog

## [3.2.3](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.2.2...v3.2.3) (2026-07-22)


### Performance Improvements

* **dc:** pack the bit-backend band sign planes once per track! call ([c57c926](https://github.com/JuliaGNSS/Tracking.jl/commit/c57c926e4d0652da6337487c7b5f2cbdab0a1fa0))

## [3.2.2](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.2.1...v3.2.2) (2026-07-15)


### Bug Fixes

* combined pilot+data decode for quadrature signals (GPS L5, Galileo E5a) ([32f7178](https://github.com/JuliaGNSS/Tracking.jl/commit/32f71785482c40d343f88603c044f4a03c2b7561))

## [3.2.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.2.0...v3.2.1) (2026-07-08)


### Performance Improvements

* **bit-buffer:** fix pre-sync per-code-block allocation in _buffer_find_bit ([a892769](https://github.com/JuliaGNSS/Tracking.jl/commit/a8927694793220cc8f71bf4ba796e002987e7e0e))

# [3.2.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.1.1...v3.2.0) (2026-07-08)


### Bug Fixes

* **twobit:** keep the threaded band pack off Polyester's closure fallback ([f453bbf](https://github.com/JuliaGNSS/Tracking.jl/commit/f453bbfdef7d01f2be629f7b02529e1826b73404))


### Features

* **tracking:** two-bit downconvert + correlate on SinCosLUT's native 2-bit NCO ([732d61c](https://github.com/JuliaGNSS/Tracking.jl/commit/732d61cd758a5d5a9b2db7d42635447b711eabdb)), closes [#162](https://github.com/JuliaGNSS/Tracking.jl/issues/162) [#160](https://github.com/JuliaGNSS/Tracking.jl/issues/160)


### Performance Improvements

* **twobit:** NEON addp-tree lowering for the sign and magnitude masks ([213ece9](https://github.com/JuliaGNSS/Tracking.jl/commit/213ece9729dd10fa6b906fdb88f3abcea8dd8650))

## [3.1.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.1.0...v3.1.1) (2026-07-08)


### Performance Improvements

* **int16:** NEON SMLAL widening accumulate on aarch64 ([9394948](https://github.com/JuliaGNSS/Tracking.jl/commit/9394948975603e70bdb729aba89ea5823ccdf223))
* **onebit:** faster NEON sign-mask via addp tree on aarch64 ([dd3d599](https://github.com/JuliaGNSS/Tracking.jl/commit/dd3d5995f8137c1ea104eefdf7e03a3e08bafa9e))

# [3.1.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.0.2...v3.1.0) (2026-07-07)


### Bug Fixes

* **bench:** rebuild track! state per sample so the loop can't drift to NaN ([ca4f0d7](https://github.com/JuliaGNSS/Tracking.jl/commit/ca4f0d74c673af917f64bddaccd3e4949e1c484b))
* **onebit:** adapt to master's get_band_id + required Int16 max_meas ([50b86ff](https://github.com/JuliaGNSS/Tracking.jl/commit/50b86ff8cd5fda3b7ae3baf217d775633920e245))
* **onebit:** correct tap funnel shift for sample-shift offset ≥ 64 ([2f7998c](https://github.com/JuliaGNSS/Tracking.jl/commit/2f7998cf08d7ce503ef16c2470a05a7e2c40cbb1))
* **onebit:** gate on modulation, not code element type, to reject CBOC ([1e92533](https://github.com/JuliaGNSS/Tracking.jl/commit/1e92533a20dd1b1f63cfe88f16528075f2b5ec82))
* **onebit:** gate shared-measurement packing to >1 sat (no single-sat regression) ([0280dd7](https://github.com/JuliaGNSS/Tracking.jl/commit/0280dd7fa41eafcdfefb8b615b20e4190f9081ad)), closes [#1](https://github.com/JuliaGNSS/Tracking.jl/issues/1) [3/#4](https://github.com/JuliaGNSS/Tracking.jl/issues/4)


### Features

* one-bit (bit-wise) downconvert + correlate backend ([e5f3eaa](https://github.com/JuliaGNSS/Tracking.jl/commit/e5f3eaaafa695b6e10c83921716ebf42e87ca8ed)), closes [#69](https://github.com/JuliaGNSS/Tracking.jl/issues/69) [90/#69](https://github.com/JuliaGNSS/Tracking.jl/issues/69)
* **onebit:** dynamic (runtime tap-count) correlator support ([e84e3eb](https://github.com/JuliaGNSS/Tracking.jl/commit/e84e3eb5195f9a0731bbc74b6ffb94ffb15e2e95)), closes [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126)


### Performance Improvements

* **onebit:** multi-signal-per-sat carrier+measurement tile-share ([f7a329e](https://github.com/JuliaGNSS/Tracking.jl/commit/f7a329ef64f336a75be1c89a48fc968b4a627e44))
* **onebit:** pack measurement sign planes once per band, share across sats ([3b1f53b](https://github.com/JuliaGNSS/Tracking.jl/commit/3b1f53ba59d3803b904360fa6bbe671afcbf3235))
* **onebit:** real multi-signal-per-sat carrier+measurement tile-share ([0ec6722](https://github.com/JuliaGNSS/Tracking.jl/commit/0ec67225453c5d9614ab9a783102b502e050ee3a))
* **onebit:** SIMD VPOPCNTQ correlation + pack-once code plane ([21ad38d](https://github.com/JuliaGNSS/Tracking.jl/commit/21ad38d683b241391b818d2918d96b5fe511174d))

## [3.0.2](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.0.1...v3.0.2) (2026-07-07)


### Bug Fixes

* **int16:** require max_meas to size the carrier amplitude ([333ccd0](https://github.com/JuliaGNSS/Tracking.jl/commit/333ccd0b9b062254a5d42fe212528fcc3cbed092))

## [3.0.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v3.0.0...v3.0.1) (2026-07-07)


### Bug Fixes

* **int16:** normalize the Int8 carrier amplitude out of the correlator ([cfe553b](https://github.com/JuliaGNSS/Tracking.jl/commit/cfe553b0d68a034ba101193add4739c97e1a0f2d)), closes [#165](https://github.com/JuliaGNSS/Tracking.jl/issues/165)
* **normalize:** divide out code amplitude for multi-level CBOC codes ([96c008e](https://github.com/JuliaGNSS/Tracking.jl/commit/96c008e9357f9b5e92b8bb8918853fef94b43391))
* **normalize:** source code amplitude from GNSSSignals.get_code_amplitude ([0fba2ee](https://github.com/JuliaGNSS/Tracking.jl/commit/0fba2ee86009d5e1cc498d53fc03fc90f1d30294))

# [3.0.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v2.3.0...v3.0.0) (2026-07-06)


* feat(bands)!: key multi-band measurements by GNSSSignals.get_band_id ([3d9a984](https://github.com/JuliaGNSS/Tracking.jl/commit/3d9a9844bb7572953dfdc74c4e157d9050e8b5ff))


### BREAKING CHANGES

* multi-band `track`/`track!` measurement NamedTuples are
now keyed by `GNSSSignals.get_band_id` (`:L1`, `:L2`, `:L5`) instead of
Tracking's lowercase `band_key` (`:l1`, `:l2`, `:l5`). `Tracking.band_key`
is removed — use `GNSSSignals.get_band_id`. `band_keys` now returns the
capitalized ids.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>

# [2.3.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v2.2.0...v2.3.0) (2026-07-06)


### Features

* **deps:** require GNSSSignals 3.3 and adopt its get_signal_id accessor ([823ba8b](https://github.com/JuliaGNSS/Tracking.jl/commit/823ba8b60ca71e9b7bf95bf407577daafe216a42))

# [2.2.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v2.1.0...v2.2.0) (2026-07-05)


### Features

* **signals:** support GPS L2CM/L2CL/L5Q and Galileo E1C/E5a ([8e03224](https://github.com/JuliaGNSS/Tracking.jl/commit/8e03224c1f7d1b8317530dc558faac15b3e7ab72))

# [2.1.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v2.0.1...v2.1.0) (2026-07-04)


### Bug Fixes

* **int16:** cap strip-mine block so single Int32 lane can't overflow ([#166](https://github.com/JuliaGNSS/Tracking.jl/issues/166)) ([76701a1](https://github.com/JuliaGNSS/Tracking.jl/commit/76701a1e24a59418e85444c3858314e39052e8a6))
* **int16:** reject `steps` that changes the carrier engine width ([#168](https://github.com/JuliaGNSS/Tracking.jl/issues/168)) ([79ceddb](https://github.com/JuliaGNSS/Tracking.jl/commit/79ceddb50f2766f7ddbdf47ee2405d913b69bc0b))
* **int16:** shrink strip-mine block on large-max_meas fallback to bound Int32 accumulator ([#167](https://github.com/JuliaGNSS/Tracking.jl/issues/167)) ([6b8f8cb](https://github.com/JuliaGNSS/Tracking.jl/commit/6b8f8cbd068a94fa71dccfcbd589b6b91d7b1145))
* **int16:** Val{W} barrier for dynamic fallback so W const-folds (no per-block alloc) ([62e858a](https://github.com/JuliaGNSS/Tracking.jl/commit/62e858aeb9613f091a626cf4b7c8caae9d6d6304))
* **int16:** validate blk kwarg at construction ([#169](https://github.com/JuliaGNSS/Tracking.jl/issues/169)) ([adb1bb4](https://github.com/JuliaGNSS/Tracking.jl/commit/adb1bb481158d37e81fe6815f967c9608e81e049))
* **int16:** widen block-flush accumulator to Int64 before horizontal sum ([#165](https://github.com/JuliaGNSS/Tracking.jl/issues/165)) ([3cf1c95](https://github.com/JuliaGNSS/Tracking.jl/commit/3cf1c95b95cff60148787220391fb45e6e6fde56))
* **test:** repin Acquisition to the reachable 2.6.0 release commit ([b4fcb22](https://github.com/JuliaGNSS/Tracking.jl/commit/b4fcb224f481685232d40ea9c767c99c44b62d36))


### Features

* **int16:** dynamic (runtime) correlator tap count (Phase 3) ([fd0cca2](https://github.com/JuliaGNSS/Tracking.jl/commit/fd0cca2657c2903b5cb1b8abf62c7d827b4a46e4)), closes [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126)
* **int16:** integer hybrid-blocked downconvert+correlate backend (Phases 0–2) ([028ef11](https://github.com/JuliaGNSS/Tracking.jl/commit/028ef11d741db74c2d0546e71609f77773f58e23)), closes [#90](https://github.com/JuliaGNSS/Tracking.jl/issues/90)
* **int16:** multi-signal-per-sat tile-share kernel ([b87e8ac](https://github.com/JuliaGNSS/Tracking.jl/commit/b87e8acf4f76da0af2d71a682ac070b795a17837))
* **int16:** multiple antenna channels (Phase 4) ([5cfb2d7](https://github.com/JuliaGNSS/Tracking.jl/commit/5cfb2d79c22c6befe31237253d3241fa16641ca2))
* **int16:** user-configurable measurement amplitude; fix CBOC Int32 overflow ([638bb61](https://github.com/JuliaGNSS/Tracking.jl/commit/638bb61ff4594dbf0ba8f9938faebea613405acc))


### Performance Improvements

* **int16:** make dynamic-shifts fallback scratch-reusing with per-block flush ([002a415](https://github.com/JuliaGNSS/Tracking.jl/commit/002a415a56a2846bd98acb57ceaded48bf410f99))
* **int16:** one-shot per-block code gen (alloc-free) — fixes low-OSR tie ([e1ee142](https://github.com/JuliaGNSS/Tracking.jl/commit/e1ee142351188f4bb0f623bfaab77196de215ab5))

## [2.0.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v2.0.0...v2.0.1) (2026-07-04)


### Bug Fixes

* **remove_satellite:** keep the satellite Dictionary hole-free so track! is safe ([#182](https://github.com/JuliaGNSS/Tracking.jl/issues/182)) ([f6c4c53](https://github.com/JuliaGNSS/Tracking.jl/commit/f6c4c53a55a690543a5c78481c9b1499be2f24d6))

# [2.0.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.5.0...v2.0.0) (2026-06-17)


* build(deps)!: bump to GNSSSignals v2 and update for renamed API ([de2acc9](https://github.com/JuliaGNSS/Tracking.jl/commit/de2acc9806dd649b35df91d547494b131f83956a))
* feat!: add allocation-free track! and TrackedSat wrapper ([8b9c634](https://github.com/JuliaGNSS/Tracking.jl/commit/8b9c634a3ebfeaec2682f11452b769bbc9e69d5d))
* feat!: add remove_satellite!/remove_satellite, drop filter_out_sats ([42e1f93](https://github.com/JuliaGNSS/Tracking.jl/commit/42e1f9395c45750712a0f4d07c2647d9902fef21))
* feat!: signals=(...) constructor; add_satellite!/add_satellite API ([910804b](https://github.com/JuliaGNSS/Tracking.jl/commit/910804baa0022f47357aebb9039dd7f6448ff9a7))
* feat(tracking)!: per-signal coherent integration length with auto-scaled loop bandwidth ([66ed6e3](https://github.com/JuliaGNSS/Tracking.jl/commit/66ed6e37ac0606fc6bb32937aaec010fc9255476))
* fix(tracking)!: infer or require group= in the by-hand add/remove API ([3bffef4](https://github.com/JuliaGNSS/Tracking.jl/commit/3bffef45b9725a17cf0a6567dcdd626488492c62))
* refactor!: cached per-thread scratch buffers, drop Bumper dep ([ec5b778](https://github.com/JuliaGNSS/Tracking.jl/commit/ec5b77849ee187941cfdfd3b7692666209764274))
* refactor!: drop prewarm! and unexport wrap_sats ([af1426a](https://github.com/JuliaGNSS/Tracking.jl/commit/af1426a1dd4b44ea3c7dc7c68fe76ac3775e723a))
* refactor!: drop Val(sampling_frequency) plumbing for GNSSSignals 1.0.3 ([018dd67](https://github.com/JuliaGNSS/Tracking.jl/commit/018dd67682798a7c4c1caeba885ccce592fb1e27)), closes [#52](https://github.com/JuliaGNSS/Tracking.jl/issues/52)
* refactor!: error on mismatched doppler_estimator type; drop silent reseed ([40b2240](https://github.com/JuliaGNSS/Tracking.jl/commit/40b22403b29a5aad3adb095d7f5b990aba59c7eb))
* refactor!: generalize TrackedSat.signals to Tuple{Vararg{TrackedSignal}} ([6894d56](https://github.com/JuliaGNSS/Tracking.jl/commit/6894d56ff4b464f709ac02772f1e41928d3a9c2f))
* refactor!: introduce TrackedSignal; flatten TrackedSat; drop TrackedSystem ([2943a81](https://github.com/JuliaGNSS/Tracking.jl/commit/2943a81abc04716cf19e2613a300556260b6d7a7))
* refactor!: remove CUDA extension ([c37b4cd](https://github.com/JuliaGNSS/Tracking.jl/commit/c37b4cdf2d25d56b27a9368e2f6d6677935be83b))
* refactor!: rename SystemSatsState → TrackedSystem (and consequents) ([de04101](https://github.com/JuliaGNSS/Tracking.jl/commit/de04101f3214d5fe55678f70069695b9ff3fdbed))


### Bug Fixes

* address accessor gaps, error messages, and review polish from [#113](https://github.com/JuliaGNSS/Tracking.jl/issues/113) ([#134](https://github.com/JuliaGNSS/Tracking.jl/issues/134)) ([c8fa8f8](https://github.com/JuliaGNSS/Tracking.jl/commit/c8fa8f8c4f9eae60efadeed2f87e813aa8e4b18e))
* **api:** settle public-API naming before v2 ([#132](https://github.com/JuliaGNSS/Tracking.jl/issues/132)) ([57b1c75](https://github.com/JuliaGNSS/Tracking.jl/commit/57b1c75956fff27f8124889cd2a082ba4ab5d98f)), closes [#113](https://github.com/JuliaGNSS/Tracking.jl/issues/113)
* **bench:** make benchmarks.jl portable across master/new constructor signatures ([f105c06](https://github.com/JuliaGNSS/Tracking.jl/commit/f105c0649c02ad2f9e57b2c9277f1f2a8661b177))
* **benchmark:** build per-signal-typed found BitBuffer in steady-state setup ([785cbdc](https://github.com/JuliaGNSS/Tracking.jl/commit/785cbdc3e605c2c311b708f6f4a83a77cd9e3d1c))
* **bit_buffer:** apply lock polarity to pre-sync recovered bits ([35c6c4f](https://github.com/JuliaGNSS/Tracking.jl/commit/35c6c4fda4b59b4f122434afac76144b23b52c02)), closes [#127](https://github.com/JuliaGNSS/Tracking.jl/issues/127)
* **bit_buffer:** clamp within-bin residual to keep CFAR variance real ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([0d62da8](https://github.com/JuliaGNSS/Tracking.jl/commit/0d62da8da6754c328b00ce5ae727fa207cfcdbfb)), closes [hi#SNR](https://github.com/hi/issues/SNR)
* **bit_buffer:** decode post-sync data bits at 1-block integration ([#125](https://github.com/JuliaGNSS/Tracking.jl/issues/125)) ([1949aab](https://github.com/JuliaGNSS/Tracking.jl/commit/1949aab80c60f06d6c86823a7aa2987c3aed3454))
* **bit_buffer:** hoist Acklam constants; drop version-fragile alloc test ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([f839449](https://github.com/JuliaGNSS/Tracking.jl/commit/f83944935fe6ba29e748f0173342950bc322f8b8))
* **bit_buffer:** reject/clamp non-divisor coherent-integration length ([a2519ac](https://github.com/JuliaGNSS/Tracking.jl/commit/a2519acb8b746c42372e08a020a082a85a15e776)), closes [#128](https://github.com/JuliaGNSS/Tracking.jl/issues/128)
* **bit_buffer:** Welford variance, confidence clamp, signal-agnostic soft bit-edge ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([7381cff](https://github.com/JuliaGNSS/Tracking.jl/commit/7381cff6b1b1a748089ba1ff6a6b74e2c8c58bf0))
* **ci:** docs, exports, and Julia 1.10 alloc-test gating ([3c70d4d](https://github.com/JuliaGNSS/Tracking.jl/commit/3c70d4dbbd688d8b0be3da195ca55b1eee5e4484))
* **code_phase:** widen wrap to one full symbol period post-sync ([14a4f63](https://github.com/JuliaGNSS/Tracking.jl/commit/14a4f637e20e361203e2ab62b2f2b2f03a8f1283))
* **correlator:** use narrow 0.1-chip DLL spacing for GPS L1C-D and L1C-P ([4a1b747](https://github.com/JuliaGNSS/Tracking.jl/commit/4a1b7473ef4e91fdaa90bebb0652f7e368f1fcc1))
* **docs:** add [@docs](https://github.com/docs) blocks for default_carrier/code_loop_filter_bandwidth ([398eaec](https://github.com/JuliaGNSS/Tracking.jl/commit/398eaec5c7187130d6161d237510d9d41e690892))
* **docs:** drop [@ref](https://github.com/ref) to internal ScratchBuffers from public docstrings ([9dc7985](https://github.com/JuliaGNSS/Tracking.jl/commit/9dc7985a96cc64eb0b4b6f1daefd70ab9537e133))
* **ext:** rename GPU kernel-launch helper, update GPU test for TrackedSat ([a531337](https://github.com/JuliaGNSS/Tracking.jl/commit/a5313378cf85ebc503b3032fb2698b28e198d268))
* **ext:** unwrap TrackedSat in GPUDownconvertAndCorrelator constructor ([72531f7](https://github.com/JuliaGNSS/Tracking.jl/commit/72531f72009d462673442e6c2d4f8ca3a2e71b8a))
* **galileo_e1b:** drop bit-edge search — symbol is aligned with primary code ([f99775f](https://github.com/JuliaGNSS/Tracking.jl/commit/f99775f5689966de6c7298bca9762b38cd2a243f))
* **gpsl1ca:** edge-lock the L1CA bit-sync detector ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([f660a28](https://github.com/JuliaGNSS/Tracking.jl/commit/f660a2837e4117cd5de2e7524c7403e48f59190c))
* **gpsl1ca:** replace L1CA bit-sync with soft-decision CFAR detector ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([fa7f783](https://github.com/JuliaGNSS/Tracking.jl/commit/fa7f78318c38004dd8779c8e80e37b4c6880bf1c)), closes [#138](https://github.com/JuliaGNSS/Tracking.jl/issues/138)
* **measurement:** promote frequency types and fix mixed-precision duration check ([b585f42](https://github.com/JuliaGNSS/Tracking.jl/commit/b585f42c6d30b7e571589a052a6e0737dd7230a2)), closes [#131](https://github.com/JuliaGNSS/Tracking.jl/issues/131)
* **sat_state:** detach Dictionary Indices in immutable slot-vector copy ([c9794e9](https://github.com/JuliaGNSS/Tracking.jl/commit/c9794e930ce4f7be95a0709554c39e1e2a4fc237)), closes [#123](https://github.com/JuliaGNSS/Tracking.jl/issues/123)
* **sat_state:** validate SignalGroup band/chip-rate and use lcm for code wrap ([f852045](https://github.com/JuliaGNSS/Tracking.jl/commit/f85204544d21c8bd353c656a628279f162363cff)), closes [#129](https://github.com/JuliaGNSS/Tracking.jl/issues/129)
* **test:** measure track! allocations from typed scope ([7d7c2d3](https://github.com/JuliaGNSS/Tracking.jl/commit/7d7c2d337d304ffd47c0bc425606f5ada13c4710))
* **tracking:** actionable error for bands without a band_key method ([9a2b5e2](https://github.com/JuliaGNSS/Tracking.jl/commit/9a2b5e206fff8c6ccfdf0fe5e7e33151b8d5e111))
* **tracking:** honor rebuilt estimator from add_satellite! handoff ([328a32b](https://github.com/JuliaGNSS/Tracking.jl/commit/328a32b0c64fd13f3152781bbf5f7ab8f4879b5e)), closes [#130](https://github.com/JuliaGNSS/Tracking.jl/issues/130)
* **tracking:** index code replica by absolute start_sample in tile-share kernel ([543cf10](https://github.com/JuliaGNSS/Tracking.jl/commit/543cf108852efc79e872678b01a55f8a36a8e939))
* **tracking:** index dynamic-shifts fused kernel by absolute start_sample; unify CPU per-signal path ([74206e1](https://github.com/JuliaGNSS/Tracking.jl/commit/74206e186970831107e82793b31f3116c902de8c)), closes [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126) [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126) [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126)
* **tracking:** make fused correlator width-independent ([#152](https://github.com/JuliaGNSS/Tracking.jl/issues/152)) ([300c397](https://github.com/JuliaGNSS/Tracking.jl/commit/300c397a24b368ea2b414840f9f82e51caf59036))
* **tracking:** preserve secondary-code phase in multi-block coherent integration ([cde8e53](https://github.com/JuliaGNSS/Tracking.jl/commit/cde8e539742a25eca1fcbc6c156d84b5c8698509))
* **tracking:** prevent deadlock when chunks are one code period ([#117](https://github.com/JuliaGNSS/Tracking.jl/issues/117)) ([83681d4](https://github.com/JuliaGNSS/Tracking.jl/commit/83681d493426b4bb0ea9a002161b329c31b18fec))
* **tracking:** reject ambiguous duplicate signal-type selector ([647eab3](https://github.com/JuliaGNSS/Tracking.jl/commit/647eab3ecd12379f3419cf9275ce5e9156d53916))
* **tracking:** throw a consistent KeyError from remove_satellite! ([d00d5f0](https://github.com/JuliaGNSS/Tracking.jl/commit/d00d5f04fbe0eb8e90ac2e627ca0da9c041a83ed))
* **tracking:** unify secondary-code sync into a phase-recovering rotation search ([d63bec7](https://github.com/JuliaGNSS/Tracking.jl/commit/d63bec704f1eb43d8155af508d98cbb12496af46))


### Features

* Acquisition-aware TrackState/add_satellite! constructors ([3ae9110](https://github.com/JuliaGNSS/Tracking.jl/commit/3ae9110ec8cf1875d68a07a360605d45e08483d3))
* **acquisition:** auto-route add_satellite! by signal type when group= omitted ([525f86d](https://github.com/JuliaGNSS/Tracking.jl/commit/525f86d7f2682ba70d2958111719ecdb8b502746))
* add update_estimator_on_handoff hook for shared estimator state ([8daba12](https://github.com/JuliaGNSS/Tracking.jl/commit/8daba1245cd81f8a9ef7b0d0f0c8642cf9d28d0e))
* **bit_buffer:** output soft bits alongside hard bits ([08c6602](https://github.com/JuliaGNSS/Tracking.jl/commit/08c6602526a272ffc49c7ee8b6b4f3a6e2e84ff9))
* **bit_buffer:** per-signal sync-search buffer width via trait ([f840c49](https://github.com/JuliaGNSS/Tracking.jl/commit/f840c494482ab780566235bf52cee45a22275c36))
* **bit_buffer:** user-adjustable sync tolerance via dispatch ([16baf56](https://github.com/JuliaGNSS/Tracking.jl/commit/16baf5675cb86a919862a6835ea5787e67c27e85))
* **galileo_e1b:** support GalileoE1B_BOC11 alongside GalileoE1B ([5f4feac](https://github.com/JuliaGNSS/Tracking.jl/commit/5f4feac4bf96ed426f79b85c26c596d5914b41ee)), closes [GNSSSignals.jl#59](https://github.com/GNSSSignals.jl/issues/59)
* GPSL1C_D/P signal support, per-signal loop-bandwidth defaults, doc fixes ([81c0e29](https://github.com/JuliaGNSS/Tracking.jl/commit/81c0e29644d9aca915f02fde9196e999ed26eaf0))
* multi-antenna support + antenna-outer correlate for tuple kernel ([a9e26fc](https://github.com/JuliaGNSS/Tracking.jl/commit/a9e26fcb2f6c56d9a7bda2456ef5b5625ae132b7))
* per-group NumAnts via SignalGroup entries + same-band agreement check ([efc0ed5](https://github.com/JuliaGNSS/Tracking.jl/commit/efc0ed5448c37b484d6b35e51a1f109612d3bf1f))
* per-signal accessors with integer or signal-type selector ([9fbd928](https://github.com/JuliaGNSS/Tracking.jl/commit/9fbd92857598ec9eebdf9f8e1a1ef37e7315a008))
* singular `signal` kwarg on TrackState for one-signal shortcut ([78afaa8](https://github.com/JuliaGNSS/Tracking.jl/commit/78afaa8a95af96d63a95210531761574958acc7d))
* **sync:** Hamming tolerance on bit-edge / NH10 detectors ([abe8053](https://github.com/JuliaGNSS/Tracking.jl/commit/abe8053a2e56b9883fb03d7741a1e53d39832e66))
* **sync:** L1C-D one-symbol-per-block sync (trivial) ([9e8add5](https://github.com/JuliaGNSS/Tracking.jl/commit/9e8add59801c21372e82e21e38cb7011c330c9db))
* **sync:** L1C-P 1800-chip overlay-code phase search ([33439d3](https://github.com/JuliaGNSS/Tracking.jl/commit/33439d36a38be686b5671f141d15ed340239fafa))
* **tracking:** add reset_loop_filters! for clean integration-length handoff ([0bb7965](https://github.com/JuliaGNSS/Tracking.jl/commit/0bb796513f9170768e3f66c2368137166243f517))
* **tracking:** auto-size loop bandwidth per estimator-driver signal ([5556104](https://github.com/JuliaGNSS/Tracking.jl/commit/5556104f3db64bb75573853d30eca3a707ac866b))
* **tracking:** immutable batch add_satellite for acquisition results ([d2adb2e](https://github.com/JuliaGNSS/Tracking.jl/commit/d2adb2e839b30c0b74594021a544fb5f88f132a4))


### Performance Improvements

* avoid Set allocations in _validate_measurements ([62627b8](https://github.com/JuliaGNSS/Tracking.jl/commit/62627b8a52f0ebaf0796343c2431fb58f0fa9751))
* **bit_buffer:** incremental O(L) soft bit-edge detector + drift-robust noise ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([0706462](https://github.com/JuliaGNSS/Tracking.jl/commit/070646211337d6775aa79238c713ea69cec9625e))
* **bit_buffer:** make the soft bit-edge detector allocation-free ([#124](https://github.com/JuliaGNSS/Tracking.jl/issues/124)) ([54ad9a4](https://github.com/JuliaGNSS/Tracking.jl/commit/54ad9a45a64eedb59238860217a28aee1c667dd8))
* **fused:** wire ScratchBuffers tile_re/tile_im through fused fallback ([c43d7c6](https://github.com/JuliaGNSS/Tracking.jl/commit/c43d7c6846b41d57fd9a752ffe091eb21067116d))
* **multi-system:** walk system tuple recursively, avoid boxing ([57cdcfb](https://github.com/JuliaGNSS/Tracking.jl/commit/57cdcfbc4f588b2680991b30d4b44a7268568ff8))
* **sat_state:** detach Dictionary Indices once at the track boundary, not per loop iteration ([c50e730](https://github.com/JuliaGNSS/Tracking.jl/commit/c50e730389fe70bc87c45b5f9f8dc59ad865e8a0)), closes [#123](https://github.com/JuliaGNSS/Tracking.jl/issues/123) [#123](https://github.com/JuliaGNSS/Tracking.jl/issues/123)
* tile-share fused kernel for multi-signal sats ([e0538ec](https://github.com/JuliaGNSS/Tracking.jl/commit/e0538ec8372d4addc627a5ff3818368a7bdfc2df))
* **track!:** eliminate per-call Core.Box allocation in bit_buffer ([1a32840](https://github.com/JuliaGNSS/Tracking.jl/commit/1a32840bc61bde24671d8e2975ea96184f57def2))
* **tracking:** keep dynamic fused kernel buffer-taking only to avoid benchmark regression ([6d4e903](https://github.com/JuliaGNSS/Tracking.jl/commit/6d4e9031de9ae97c7723c99178d7972b626b356f)), closes [#126](https://github.com/JuliaGNSS/Tracking.jl/issues/126)
* tuple-recursive antenna-shape validation to keep track! alloc-free ([0f438a1](https://github.com/JuliaGNSS/Tracking.jl/commit/0f438a17ed9f13f366e1cb36e627b43e8966019d))


### BREAKING CHANGES

* the keyword add/remove API no longer assumes a `:default`
group; multi-group TrackStates must pass `group=` explicitly (previously a
raw KeyError) and single named groups are now inferred.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
* `track` / `track!` no longer accept a
`preferred_num_code_blocks_to_integrate` keyword, and the internal
`downconvert_and_correlate[!]` / `estimate_dopplers_and_filter_prompt[!]`
signatures drop the trailing parameter — custom `AbstractDopplerEstimator`
implementations of `estimate_dopplers_and_filter_prompt[!]` must drop it too
and read each signal's field instead. Configure integration length via
`set_preferred_num_code_blocks_to_integrate!` (or the `TrackedSignal`
constructor kwarg) before tracking.

The conventional PLL/DLL estimator auto-scales each signal's effective loop
bandwidth by `1/N` for its integration length `N`, holding the loop's `BL·Δt`
stability product at its one-primary-period value. The configured bandwidth is
a per-primary-period reference; longer coherent integration (e.g. L5I 10 ms)
stays stable with no manual bandwidth change. For N=1 it divides by 1 and is
bit-identical.

Also scopes the secondary-/overlay-aware code-replica + integration-window wrap
to multi-block integration only (`num_blocks > 1`); a single block keeps the
primary-only replica so the per-block secondary sign stays handled by the bit
buffer — restoring bit-identical N=1 behavior and fixing the GPS L5I
secondary-code phase-recovery test that the unscoped wrap perturbed.

Full test suite (incl. Aqua) passes. Verified on a real GPS L5I capture:
per-signal mixed lengths (PRN 1 at 10 ms, others at 1 ms) and a staged
1 ms->10 ms promotion both reach full coherent gain with no manual tuning.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
* TrackedSat instances must be constructed with the
same doppler_estimator that the TrackState (or merge_sats target) uses.
Code that built sats with the default estimator and then passed them to
a TrackState with a custom estimator must now thread the estimator
through both calls explicitly. The keyword `add_satellite!(; prn, …)`
form is unaffected.
* filter_out_sats is removed. Use remove_satellite! (or
remove_satellite) instead. The new API takes one PRN at a time; for
batched removal, loop over the PRN list.
* `TrackState` now has three type parameters
(`{S, SG, DE}`) instead of two; external callers writing
`TrackState{S, DE}` type annotations must update them. The
`get_tracked_system` accessor was already removed in an earlier commit.
* `_update_tracked_sat_correlator`,
`_update_tracked_sat_doppler`, and `update` no longer accept a positional
`system` argument; the per-signal system is recovered from each
TrackedSignal's `.signal` field.
* SatState, TrackedSystem, get_tracked_system, and the
old wrapper-shape TrackedSat are gone. Callers that built SatState
directly must use `TrackedSat(signal, prn, code_phase, carrier_doppler;
doppler_estimator = ...)` instead. Code that read `sat.sat_state.x` or
`sat.estimator_state.y` now reads `sat.x` and `sat.doppler_estimator_state.y`
respectively.
* Tracking.jl no longer supports GNSSSignals v1 or
earlier. Downstream users must update their own signal-type names
and AbstractGNSS uses to match the new API.
* GPU tracking via TrackingCUDAExt is no longer
available. Use the CPU path for now; KernelAbstractions-based GPU
support will land in a follow-up.
* `SystemSatsState` is renamed to `TrackedSystem`.
Code that imports `Tracking.SystemSatsState` or constructs it
directly must update the symbol name.
* `MultipleSystemSatsState` is renamed to
`TrackedSystems`. Custom estimators with method signatures
constraining `MultipleSystemSatsState` need to update to
`TrackedSystems`.
* `get_system_sats_state` is renamed to
`get_tracked_system`. External callers must update.
* `TrackState`'s field
`multiple_system_sats_state::S` is renamed to `tracked_systems::S`.
The keyword argument on `TrackState(track_state; ...)` is also
renamed accordingly. Code that reaches into
`track_state.multiple_system_sats_state` directly must update
to `track_state.tracked_systems`.
* `MultipleSystemType` is no longer exported (and
the alias is gone from the source). The GPU extension inlines the
single use site.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
* `CPUDownconvertAndCorrelator` lost its `{B}` type
parameter and `buffer` field. The constructor still works as
`CPUDownconvertAndCorrelator()` (zero-arg), but external code that
matched on `CPUDownconvertAndCorrelator{B}` or read
`.buffer::SlabBuffer` will break.
* `CPUThreadedDownconvertAndCorrelator` lost its
`buffers::Vector{SlabBuffer}` field; the new field is
`buffers::Matrix{Vector{UInt8}}`. External code that read the
`buffers` field directly will need to adapt.
* `Bumper` is no longer a dependency. Code that relied
on Tracking re-exporting Bumper symbols (none did, but listed for
completeness) needs to import Bumper directly.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
* `prewarm!` is removed from the public API. Callers
should drop the `prewarm!(track_state, ...)` line from their setup;
the single-iteration warmup now happens automatically on the first
`track!` call.
* `wrap_sats` is no longer exported. Users who need
it (rare — the standard path is `SystemSatsState(estimator, system,
sats)`) must qualify as `Tracking.wrap_sats`.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
* `CPUDownconvertAndCorrelator` and
`CPUThreadedDownconvertAndCorrelator` no longer have an `MESF` type
parameter and their constructors no longer accept
`Val(sampling_frequency)`. Replace
`CPUDownconvertAndCorrelator(Val(sampling_frequency))` with
`CPUDownconvertAndCorrelator()` (and likewise for the threaded
variant).
* The 13-arg single-sat `downconvert_and_correlate!`
and `gen_code_replica!` no longer accept the trailing
`maximum_expected_sampling_frequency::Val` argument. Callers must
drop it.
* `GNSSSignals` compat is tightened to `"1.0.3"`.
Downstream packages that pin `GNSSSignals 0.17.x` cannot use this
release of Tracking.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
* `SystemSatsState`'s value type changed from
`SatState` to `TrackedSat`. Code that built a `SystemSatsState` from
a raw `Dictionary{I, SatState}` directly must now wrap each entry in
a `TrackedSat`, e.g. via `SystemSatsState(estimator, system, sats)`
or `wrap_sats(estimator, sats)`. External callers that go through
`merge_sats`, `filter_out_sats`, `get_sat_state`, etc. are
unaffected.
* `get_sat_states(track_state)` and
`get_sat_states(track_state, system_idx)` now return
`Dictionary{I, TrackedSat}` instead of `Dictionary{I, SatState}`.
Callers iterating these dictionaries must read `.sat_state` to reach
the underlying `SatState` (or use the new `get_sat_state(track_state,
prn)` accessor, which unwraps).
* `AbstractDopplerEstimator` lost its `{N, I}` type
parameters. Custom estimators declared as `MyEstimator <:
AbstractDopplerEstimator{N, I}` must drop the parameters.
* `ConventionalPLLAndDLL` and
`ConventionalAssistedPLLAndDLL` constructors no longer take a
`multiple_system_sats_state` positional argument. The estimator is
now config-only:

    # before
    ConventionalPLLAndDLL((system_sats_state,);
        carrier_loop_filter_bandwidth = 18.0Hz)

    # after
    ConventionalPLLAndDLL(;
        carrier_loop_filter_bandwidth = 18.0Hz)

Per-satellite estimator state is built lazily via
`init_estimator_state(estimator, sat_state)` when a sat enters the
track set, and lives inside its `TrackedSat` wrapper.
* `track_state.doppler_estimator.states` field is
removed. Per-satellite estimator state now lives on
`tracked_sat.estimator_state` for each `TrackedSat` in
`SystemSatsState.states`. Use the new `get_estimator_state` accessor
(`get_estimator_state(track_state, prn)`) to read it.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>

# [1.5.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.4.0...v1.5.0) (2026-04-30)

### Features

  - **deps:** support Acquisition v2 ([f58ee8e](https://github.com/JuliaGNSS/Tracking.jl/commit/f58ee8eae23f420dda224a88a96eea1eacdef4c3))

# [1.4.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.3.0...v1.4.0) (2026-04-27)

### Features

  - **pll_and_dll:** support per-satellite loop filter bandwidths ([5e4a9d3](https://github.com/JuliaGNSS/Tracking.jl/commit/5e4a9d32d1723be5783450832b5dc610a54f9508))

# [1.3.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.2.1...v1.3.0) (2026-04-27)

### Features

  - collect filtered prompt per correlation in SatState ([8a086b7](https://github.com/JuliaGNSS/Tracking.jl/commit/8a086b77c2f96eaf4727eb93a6047a2bb816d528))

## [1.2.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.2.0...v1.2.1) (2026-03-25)

### Bug Fixes

  - use Threads.maxthreadid() instead of Threads.nthreads() ([1f7d104](https://github.com/JuliaGNSS/Tracking.jl/commit/1f7d104f706439c17757ab7b8e14f1b9059249ff))

# [1.2.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.1.1...v1.2.0) (2026-03-25)

### Bug Fixes

  - remove objectid-based system lookup in CPUThreadedDownconvertAndCorrelator ([e4fa722](https://github.com/JuliaGNSS/Tracking.jl/commit/e4fa722487cd7e750c65dde1cc4ec1207e7975ee))

### Features

  - add CPUThreadedDownconvertAndCorrelator struct and constructor ([0b6ca67](https://github.com/JuliaGNSS/Tracking.jl/commit/0b6ca676ec0e74aaaa831255d1eb20a79db21387))
  - default to CPUThreadedDownconvertAndCorrelator in track() ([bd66b92](https://github.com/JuliaGNSS/Tracking.jl/commit/bd66b9242f3acd90591fd574a4634a573934e0cc))

### Performance Improvements

  - use Polyester [@batch](https://github.com/batch) instead of Threads.[@threads](https://github.com/threads) and update benchmarks ([01e190f](https://github.com/JuliaGNSS/Tracking.jl/commit/01e190f307b3bf371402f33d79ef952ab538ba59))

## [1.1.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.1.0...v1.1.1) (2026-03-25)

### Bug Fixes

  - drop conventional_pll_and_dll.jl change (no alloc reduction) ([6c84e3b](https://github.com/JuliaGNSS/Tracking.jl/commit/6c84e3bf91d761f1f97e7df7557506340b2dd339))

### Performance Improvements

  - reduce allocations in estimate_dopplers_and_filter_prompt ([894f3db](https://github.com/JuliaGNSS/Tracking.jl/commit/894f3db7357da705ec9db1e3f8675f4bf5c27560))
  - use map_unzip to avoid intermediate Dictionary in estimate_dopplers ([7137bd7](https://github.com/JuliaGNSS/Tracking.jl/commit/7137bd7c9414a8a84d811d39b1dfaf7a91c72278))

# [1.1.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.0.1...v1.1.0) (2026-03-21)

### Bug Fixes

  - compute SIMD width at [@generated](https://github.com/generated) time, not runtime ([2caba06](https://github.com/JuliaGNSS/Tracking.jl/commit/2caba0611e89cd39f406ceaa9535d53e38a1d26a))
  - use extension test path for v1.10 Buildkite step ([827e608](https://github.com/JuliaGNSS/Tracking.jl/commit/827e6085e0c6b80d4dd3266517e75ed52f5b6181))

### Features

  - add fused downconvert + correlate SIMD kernel ([960e472](https://github.com/JuliaGNSS/Tracking.jl/commit/960e4727fe5f6a90daab29ea46a3318aebe6a5b8))
  - fused carrier gen + downconvert with [@generated](https://github.com/generated) NumAnts unrolling ([9866a33](https://github.com/JuliaGNSS/Tracking.jl/commit/9866a33b50872d5a4236475ee61f6c2ea272f9f4))
  - restore VectorizationBase for optimal SIMD width detection ([23ff237](https://github.com/JuliaGNSS/Tracking.jl/commit/23ff2377136c1265f02e829ef7e6e65df161efd8))

### Performance Improvements

  - fix AVX-512 heap allocations in SIMD helpers ([2bf3ceb](https://github.com/JuliaGNSS/Tracking.jl/commit/2bf3ceb0e87821027409955e1883c12814663036))
  - optimize fused downconvert+correlate kernel ([a6a6213](https://github.com/JuliaGNSS/Tracking.jl/commit/a6a6213555ac0cfcd06b3b7bc0c8528127180bf9))

### Reverts

  - remove .codecov.yml, coverage fixed in test instead ([481047c](https://github.com/JuliaGNSS/Tracking.jl/commit/481047c939ec43a6520be60fbb9d5c71bc4585b1))

## [1.0.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v1.0.0...v1.0.1) (2026-02-25)

### Bug Fixes

  - correct DLL discriminator normalization and update docstrings ([bb9fa6f](https://github.com/JuliaGNSS/Tracking.jl/commit/bb9fa6ff7e1dd716548e5e9ea49e0df9dd8a5b16))

# [0.18.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.17.1...v0.18.0) (2026-01-11)

### Features

  - add type-stable loop filter selection for ConventionalPLLAndDLL ([18f3bc2](https://github.com/JuliaGNSS/Tracking.jl/commit/18f3bc2a5cd714fdbff1623e05bf4a84b8974277))

## [0.17.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.17.0...v0.17.1) (2025-12-20)

### Bug Fixes

  - initial fll discriminator ([3eb5a13](https://github.com/JuliaGNSS/Tracking.jl/commit/3eb5a13decd952e72685bac4e0b799fbd0a7977b))

# [0.17.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.16.1...v0.17.0) (2025-12-19)

### Bug Fixes

  - added missing imports for fll discriminator test ([8242dab](https://github.com/JuliaGNSS/Tracking.jl/commit/8242dab6840f3f95f19683a5c7ce0f34b4eaca7a))
  - correct units in fll discriminator test ([1b19ed1](https://github.com/JuliaGNSS/Tracking.jl/commit/1b19ed182faca872c81f84b1250586cda3f09ffd))
  - simplify fll discriminator ([51ba424](https://github.com/JuliaGNSS/Tracking.jl/commit/51ba424ac1295ee6ce0dd6febc6fab9cbe16d5da))
  - simplify units of fll discriminator ([28dfa8f](https://github.com/JuliaGNSS/Tracking.jl/commit/28dfa8f04cdc607098e69d5e54d7507dc4570eae))

### Features

  - fll discriminator ([b9b7124](https://github.com/JuliaGNSS/Tracking.jl/commit/b9b7124af03f1a5626330eed882f934a0d46d156))

## [0.16.1](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.16.0...v0.16.1) (2025-12-09)

### Bug Fixes

  - update acquisition to 0.2 ([1c756c0](https://github.com/JuliaGNSS/Tracking.jl/commit/1c756c0c9e70ce6e71e38a491f41e93643746eda))

# [0.16.0](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.15.8...v0.16.0) (2025-11-19)

### Bug Fixes

  - add Pkg.develop in CI to properly resolve Tracking dependencies ([6a7c299](https://github.com/JuliaGNSS/Tracking.jl/commit/6a7c299b17ae792bed86ed749d2beeadcc678bf4))
  - export types and functions needed by CUDA extension ([dde6f40](https://github.com/JuliaGNSS/Tracking.jl/commit/dde6f4070a6016b19bd44c044b9b0253d1c634d8))
  - import functions from Tracking to properly extend in GPU extension ([fbd09aa](https://github.com/JuliaGNSS/Tracking.jl/commit/fbd09aabeec55aa5839be9880183abe7fb6aecde))
  - import get_num_samples and update from Tracking in GPU extension ([90865cc](https://github.com/JuliaGNSS/Tracking.jl/commit/90865cccef2fa15ffdf9ebdb518f23e81558d038))
  - properly export GPU types from extension to parent module ([a50656d](https://github.com/JuliaGNSS/Tracking.jl/commit/a50656d4a6e92d6591452574f3fa348098e5cd9a))
  - remove __init__() pattern and update tests to import from TrackingCUDAExt ([49d0c44](https://github.com/JuliaGNSS/Tracking.jl/commit/49d0c446cdc49d6c8790a69a8dbf2fd6ad117d3e))
  - use Base.get_extension to access GPU types in extension tests ([3e6ccd7](https://github.com/JuliaGNSS/Tracking.jl/commit/3e6ccd74734c3b06d585291a66e047c924caecbd))

### Features

  - add TrackingCUDAExt extension entry point ([37f3113](https://github.com/JuliaGNSS/Tracking.jl/commit/37f3113764c24a586309ed7acdf94f088f192714))

## [0.15.8](https://github.com/JuliaGNSS/Tracking.jl/compare/v0.15.7...v0.15.8) (2025-10-04)

### Bug Fixes

  - use correct earliest and latest sample shift in gen_code_replica function ([a0be2d7](https://github.com/JuliaGNSS/Tracking.jl/commit/a0be2d7c761e124d9b49e599070f7e636a24f43b))
  - use correct latest sample shift in correlate function ([8b8b49b](https://github.com/JuliaGNSS/Tracking.jl/commit/8b8b49b9c402ae5409ab43bbd4231515d886d457))
