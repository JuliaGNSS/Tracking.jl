# One-bit (hard-limited) bit-wise downconvert + correlate backend.
#
# The bit-wise counterpart of the Float32 fused / integer paths. It takes the same
# `Complex{Int16}` (12-bit ADC) sample buffers, but keeps only the SIGN of every
# operand (measurement, carrier, code), so the carrier wipe-off collapses to XOR
# and the tap accumulate to popcount:
#
#   value +1 → bit 0, −1 → bit 1 (signbit); a ±1 product is an XOR of sign bits,
#   and Σ over N samples of a ±1 sequence is  N − 2·popcount.
#
#   DI = mᵣ·cos + mᵢ·sin,  DQ = mᵢ·cos − mᵣ·sin      (all ±1 in 1-bit form) ⇒
#   Iₓ = Σ codeₓ·DI = (N − 2·pc(codeₓ⊻mᵣ⊻cos)) + (N − 2·pc(codeₓ⊻mᵢ⊻sin))
#   Qₓ = Σ codeₓ·DQ = (N − 2·pc(codeₓ⊻mᵢ⊻cos)) − (N − 2·pc(codeₓ⊻mᵣ⊻sin))
#
# Per block: pack the measurement sign planes (per antenna), generate the carrier
# sign planes straight off the NCO's top bit (sign(sin)=MSB(acc), sign(cos)=MSB(acc+¼)),
# pack the code sign plane per tap (movemask of the GNSSSignals Int8 ±1 replica), then
# XOR + `count_ones` accumulate. Strip-mined into `blk`-sample blocks with L1-resident
# scratch reused across blocks, mirroring the integer hybrid-blocked backend.
#
# 1-bit quantisation costs ≈2 dB of SNR (plus ≈1 dB for the 1-bit carrier), the classic
# trade for bit-wise speed and 1-bit memory bandwidth. Accumulators are converted to
# `ComplexF64` (M=1) / `SVector{M,ComplexF64}` (M>1) at finalize; every downstream
# consumer (discriminators, C/N0, bit buffer) is ratio/normalised, so the absolute 1-bit
# scale is immaterial. Scope mirrors the integer backend: one signal per correlate call
# (a multi-signal sat correlates its signals in turn), static tap counts, any antenna M.

# Default strip-mine block length (samples); a multiple of 64 so every block packs a
# whole number of UInt64 words.
const _ONEBIT_BLK = 8192

# ── sign-mask helpers: pack the sign bit of 64 lanes into a UInt64 (bit j ⇔ lane j < 0) ──
@inline _ob_mm(v::SIMD.Vec{64,Int8}) = Base.llvmcall(("""
    define i64 @entry(<64 x i8> %v) #0 { %c = icmp slt <64 x i8> %v, zeroinitializer
      %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""", "entry"),
    UInt64, Tuple{NTuple{64,Base.VecElement{Int8}}}, v.data)
@inline _ob_mm(v::SIMD.Vec{64,Int16}) = Base.llvmcall(("""
    define i64 @entry(<64 x i16> %v) #0 { %c = icmp slt <64 x i16> %v, zeroinitializer
      %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""", "entry"),
    UInt64, Tuple{NTuple{64,Base.VecElement{Int16}}}, v.data)
@inline _ob_mm32(v::SIMD.Vec{64,UInt32}) = Base.llvmcall(("""
    define i64 @entry(<64 x i32> %v) #0 { %c = icmp slt <64 x i32> %v, zeroinitializer
      %m = bitcast <64 x i1> %c to i64 ret i64 %m } attributes #0={alwaysinline}""", "entry"),
    UInt64, Tuple{NTuple{64,Base.VecElement{Int32}}}, reinterpret(SIMD.Vec{64,Int32}, v).data)

@inline @generated _ob_iota() =
    :(SIMD.Vec{64,UInt32}($(Expr(:tuple, (UInt32(j) for j = 0:63)...))))
@inline _ob_masklast(w::UInt64, r::Int) = r == 0 ? w : w & ((UInt64(1) << r) - UInt64(1))

# cycles → UInt32 NCO word (wraps via two's complement for negative Doppler).
@inline _ob_word(cycles::Real) = unsafe_trunc(UInt32, round(Int64, cycles * 4.294967296e9))

# Fill `nwb` carrier sign words. Lane j of word w has NCO accumulator base+w·(64·fw)+j·fw;
# sign(sin)=MSB(acc), sign(cos)=MSB(acc+¼cycle). O(1)/word at any frequency.
@inline function _ob_fill_carrier!(
    sinw::Vector{UInt64},
    cosw::Vector{UInt64},
    base0::UInt32,
    fw::UInt32,
    nwb::Int,
    r::Int,
)
    ramp = _ob_iota() * SIMD.Vec{64,UInt32}(fw)
    step = fw * UInt32(64)
    q = SIMD.Vec{64,UInt32}(0x40000000)
    base = base0
    @inbounds for w = 1:nwb
        ph = SIMD.Vec{64,UInt32}(base) + ramp
        s = _ob_mm32(ph)
        c = _ob_mm32(ph + q)
        if w == nwb
            s = _ob_masklast(s, r);
            c = _ob_masklast(c, r)
        end
        sinw[w] = s;
        cosw[w] = c
        base += step
    end
    nothing
end

# Pack tap code sign plane: bit for output sample n = sign(extb[byteoff+n]). `extb` is
# padded ≥ byteoff+nwb·64 so the last (masked) word can be read whole.
@inline function _ob_pack_code!(
    codeb::Vector{UInt64},
    koff::Int,
    extb::Vector{Int8},
    byteoff::Int,
    nwb::Int,
    r::Int,
)
    @inbounds for w = 0:(nwb-1)
        wd = _ob_mm(SIMD.vload(SIMD.Vec{64,Int8}, extb, byteoff + (w << 6) + 1))
        codeb[koff+w+1] = w == nwb - 1 ? _ob_masklast(wd, r) : wd
    end
    nothing
end

# Pack one antenna's measurement sign planes (sign(real)→mrb, sign(imag)→mib). Full
# 64-sample words via a deinterleaving vector load; a scalar tail for the last < 64.
@inline function _ob_pack_meas!(
    mrb::Vector{UInt64},
    mib::Vector{UInt64},
    joff::Int,
    signal,
    j::Int,
    p_sig::Ptr{Int16},
    colbase::Int,
    base::Int,
    len::Int,
    nwb::Int,
    r::Int,
)
    full = fld(len, 64)
    @inbounds for w = 0:(full-1)
        n = (w << 6) + 1
        byte_off = (colbase + base + n - 2) * 2 * sizeof(Int16)
        re, im = _deinterleave_load(SIMD.Vec{64,Int16}, p_sig, byte_off)
        mrb[joff+w+1] = _ob_mm(re);
        mib[joff+w+1] = _ob_mm(im)
    end
    if r != 0
        wr = zero(UInt64);
        wi = zero(UInt64)
        @inbounds for i = 0:(r-1)
            sig = signal[base+(full<<6)+i, j]
            (real(sig) < 0) && (wr |= UInt64(1) << i)
            (imag(sig) < 0) && (wi |= UInt64(1) << i)
        end
        mrb[joff+full+1] = wr;
        mib[joff+full+1] = wi
    end
    nothing
end

# Scratch: block code buffer `extb` (Int8), carrier sign planes, per-tap code sign planes
# (`codeb`, NC·wpb), per-antenna measurement sign planes (`mrb`/`mib`, M·wpb). Grown lazily
# and reused, so a hoisted backend is allocation-free in steady state.
struct OneBitScratchBuffers
    extb::Vector{Int8}
    sinw::Vector{UInt64}
    cosw::Vector{UInt64}
    codeb::Vector{UInt64}
    mrb::Vector{UInt64}
    mib::Vector{UInt64}
end
OneBitScratchBuffers() =
    OneBitScratchBuffers(Int8[], UInt64[], UInt64[], UInt64[], UInt64[], UInt64[])

"""
$(SIGNATURES)

One-bit (hard-limited) bit-wise CPU downconvert + correlate backend
(single-threaded). Opt-in alternative to [`CPUDownconvertAndCorrelator`] for
`Complex{Int16}` sample buffers: it 1-bit-quantises the measurement, carrier and
code and correlates with XOR + popcount. Construct **once outside** the `track!`
loop and pass it via the `downconvert_and_correlator` keyword for an
allocation-free steady state. 1-bit quantisation trades ≈2–3 dB of SNR for
bit-wise speed; downstream consumers are ratio-normalised, so the coarse
amplitude is immaterial.
"""
struct OneBitDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::OneBitScratchBuffers
    blk::Int
end

"""
$(SIGNATURES)

Multi-threaded one-bit bit-wise backend. One `OneBitScratchBuffers` per thread
(indexed by `Threads.threadid()` inside `@batch`). See
[`OneBitDownconvertAndCorrelator`](@ref).
"""
struct OneBitThreadedDownconvertAndCorrelator <: AbstractDownconvertAndCorrelator
    buffers::Vector{OneBitScratchBuffers}
    blk::Int
end

OneBitDownconvertAndCorrelator(; blk::Integer = _ONEBIT_BLK) =
    OneBitDownconvertAndCorrelator(OneBitScratchBuffers(), Int(blk))

OneBitThreadedDownconvertAndCorrelator(; blk::Integer = _ONEBIT_BLK) =
    OneBitThreadedDownconvertAndCorrelator(
        [OneBitScratchBuffers() for _ = 1:Threads.maxthreadid()],
        Int(blk),
    )

const _OneBitDC =
    Union{OneBitDownconvertAndCorrelator,OneBitThreadedDownconvertAndCorrelator}

@inline _ob_scratch(dc::OneBitDownconvertAndCorrelator) = dc.buffers
@inline _ob_scratch(dc::OneBitThreadedDownconvertAndCorrelator) =
    dc.buffers[Threads.threadid()]

@inline _ob_num_ants_val(::AbstractCorrelator{M}) where {M} = NumAnts{M}()

# ── The one-bit hybrid-blocked kernel ─────────────────────────────────────────
# Returns this integration's correlation contribution: `Vector` of `NC` complex
# sums (one per tap) — `ComplexF64` (M=1) or `SVector{M,ComplexF64}` (M>1) — to be
# added to the correlator's running accumulators. `@generated` over (NC, M): the
# per-(antenna, tap) popcount accumulators live in named locals and unroll.
@generated function _onebit_hybrid_blocked!(
    dc::_OneBitDC,
    signal::AbstractVecOrMat{Complex{Int16}},
    ::NumAnts{M},
    signal_type,
    prn::Integer,
    sample_shifts::SVector{NC},
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample::Integer,
    num_samples::Integer,
) where {M,NC}
    # Accumulator declarations: A/B (→ I) and C/E (→ Q) popcounts per (antenna j, tap k).
    init = Expr(:block)
    for k = 1:NC
        push!(init.args, :($(Symbol("off_$k")) = Int(sample_shifts[$k]) - min_shift))
    end
    for j = 1:M, k = 1:NC
        for s in ("A", "B", "C", "E")
            push!(init.args, :($(Symbol("$(s)_$(j)_$k")) = zero(Int64)))
        end
    end

    # Per-word unrolled accumulate: shared sample⊻carrier products per antenna, then
    # one XOR + count_ones per (tap, product).
    body = Expr(:block)
    for j = 1:M
        push!(body.args, :(mr = mrb[$(j - 1) * wpb + w]))
        push!(body.args, :(mi = mib[$(j - 1) * wpb + w]))
        push!(body.args, :(prc = mr ⊻ cosv))
        push!(body.args, :(pis = mi ⊻ sinv))
        push!(body.args, :(pic = mi ⊻ cosv))
        push!(body.args, :(prs = mr ⊻ sinv))
        for k = 1:NC
            # tap offset is baked into the packed plane (see _ob_pack_code!); index by (k-1)·wpb+w
            push!(body.args, :(cw = codeb[$(k - 1) * wpb + w]))
            push!(body.args, :($(Symbol("A_$(j)_$k")) += count_ones(cw ⊻ prc)))
            push!(body.args, :($(Symbol("B_$(j)_$k")) += count_ones(cw ⊻ pis)))
            push!(body.args, :($(Symbol("C_$(j)_$k")) += count_ones(cw ⊻ pic)))
            push!(body.args, :($(Symbol("E_$(j)_$k")) += count_ones(cw ⊻ prs)))
        end
    end

    # Finalize: Iₓ = 2N − 2(A+B); Qₓ = 2(E − C).
    function tapval(k)
        if M == 1
            :(complex(
                Float64(2 * N - 2 * ($(Symbol("A_1_$k")) + $(Symbol("B_1_$k")))),
                Float64(2 * ($(Symbol("E_1_$k")) - $(Symbol("C_1_$k")))),
            ))
        else
            ant = [
                :(complex(
                    Float64(2 * N - 2 * ($(Symbol("A_$(j)_$k")) + $(Symbol("B_$(j)_$k")))),
                    Float64(2 * ($(Symbol("E_$(j)_$k")) - $(Symbol("C_$(j)_$k")))),
                )) for j = 1:M
            ]
            :(SVector{$M,ComplexF64}(tuple($(ant...))))
        end
    end
    ret =
        M == 1 ? :(ComplexF64[$([tapval(k) for k = 1:NC]...)]) :
        :(SVector{M,ComplexF64}[$([tapval(k) for k = 1:NC]...)])

    quote
        get_code_type(signal_type) <: Integer || throw(
            ArgumentError(
                string(
                    "OneBitDownconvertAndCorrelator supports BPSK (±1 integer code) signals ",
                    "only; got ", typeof(signal_type), " with code type ",
                    get_code_type(signal_type),
                    ". Bit-wise correlation is awkward for non-binary (CBOC/BOC) modulations.",
                ),
            ),
        )
        N = Int(num_samples)
        min_shift = minimum(sample_shifts)
        span = maximum(sample_shifts) - min_shift
        $init

        sampling_freq = Float64(upreferred(sampling_frequency / Hz))
        carrier_freq = Float64(upreferred(carrier_frequency / Hz))
        cps_car = carrier_freq / sampling_freq
        fw = _ob_word(cps_car)
        phase0 = _ob_word(Float64(carrier_phase))
        cps_code = Float64(upreferred(code_frequency / Hz)) / sampling_freq
        code_phase0 = Float64(code_phase)
        num_rows = size(signal, 1)
        p_sig = Ptr{Int16}(pointer(signal))

        bufs = _ob_scratch(dc)
        blk = dc.blk
        wpb = cld(blk, 64)                      # UInt64 words per block plane
        length(bufs.extb) < blk + span + 64 && resize!(bufs.extb, blk + span + 64)
        length(bufs.sinw) < wpb && resize!(bufs.sinw, wpb)
        length(bufs.cosw) < wpb && resize!(bufs.cosw, wpb)
        length(bufs.codeb) < $NC * wpb && resize!(bufs.codeb, $NC * wpb)
        length(bufs.mrb) < $M * wpb && resize!(bufs.mrb, $M * wpb)
        length(bufs.mib) < $M * wpb && resize!(bufs.mib, $M * wpb)
        extb = bufs.extb;
        sinw = bufs.sinw;
        cosw = bufs.cosw
        codeb = bufs.codeb;
        mrb = bufs.mrb;
        mib = bufs.mib

        blk_off = 0
        @inbounds while blk_off < N
            len = min(blk, N - blk_off)
            nwb = cld(len, 64)
            r = len & 63
            # code (Int8 ±1) for samples [min_shift, len+span), then per-tap sign planes
            gen_code!(
                view(extb, 1:(len+span)),
                signal_type,
                prn,
                sampling_frequency,
                code_frequency,
                code_phase0 + cps_code * blk_off,
                min_shift,
            )
            $(Expr(:block, (:(_ob_pack_code!(codeb, ($k - 1) * wpb, extb, $(Symbol("off_$k")), nwb, r)) for k = 1:NC)...))
            # carrier sign planes (shared across antennas & taps)
            base0 = phase0 + fw * UInt32(blk_off)
            _ob_fill_carrier!(sinw, cosw, base0, fw, nwb, r)
            # measurement sign planes, per antenna
            base = signal_start_sample + blk_off
            $(Expr(:block, (:(_ob_pack_meas!(mrb, mib, $(j - 1) * wpb, signal, $j, p_sig, $(j - 1) * num_rows, base, len, nwb, r)) for j = 1:M)...))
            # XOR + popcount accumulate
            for w = 1:nwb
                cosv = cosw[w];
                sinv = sinw[w]
                $body
            end
            blk_off += len
        end
        $ret
    end
end

# ── Correlate / plumbing (mirrors the CPU/Int16 backends) ─────────────────────
# Single signal per sat: one kernel call.
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal},
    per_signal_completed::Tuple{Bool},
    dc::_OneBitDC,
    signal,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    head = signals[1]
    p = _signal_replica_params(head, code_doppler, code_phase, sampling_frequency, num_samples_signal)
    new_acc = _onebit_hybrid_blocked!(
        dc, signal, _ob_num_ants_val(head.correlator), head.signal, prn,
        p.sample_shifts, p.signal_code_phase, carrier_phase, p.code_frequency,
        carrier_frequency, sampling_frequency, signal_start_sample, samples_to_integrate,
    )
    new_corr = update_accumulator(head.correlator, get_accumulators(head.correlator) .+ new_acc)
    ((new_corr, per_signal_completed[1]),)
end

# Multiple signals per sat: correlate each in turn (no carrier-share tile; correct,
# simpler than the integer tile-share — a future optimisation).
@inline function _correlate_signals(
    signals::Tuple{TrackedSignal,TrackedSignal,Vararg{TrackedSignal}},
    per_signal_completed::Tuple,
    dc::_OneBitDC,
    signal,
    code_doppler,
    code_phase,
    carrier_frequency,
    carrier_phase,
    sampling_frequency,
    signal_start_sample,
    samples_to_integrate,
    prn,
    num_samples_signal,
)
    new_corrs = map(signals) do head
        p = _signal_replica_params(head, code_doppler, code_phase, sampling_frequency, num_samples_signal)
        new_acc = _onebit_hybrid_blocked!(
            dc, signal, _ob_num_ants_val(head.correlator), head.signal, prn,
            p.sample_shifts, p.signal_code_phase, carrier_phase, p.code_frequency,
            carrier_frequency, sampling_frequency, signal_start_sample, samples_to_integrate,
        )
        update_accumulator(head.correlator, get_accumulators(head.correlator) .+ new_acc)
    end
    map(tuple, new_corrs, per_signal_completed)
end

function _update_tracked_sat_correlator(
    sat::TrackedSat,
    dc::_OneBitDC,
    signal,
    num_samples_signal,
    sampling_frequency,
    intermediate_frequency,
)
    samples_to_integrate, per_signal_completed = _calc_min_samples_and_completed(
        sat.signals, sat.signal_start_sample, sampling_frequency,
        sat.code_doppler, sat.code_phase, num_samples_signal,
    )
    samples_to_integrate == 0 && return sat
    carrier_frequency = sat.carrier_doppler + intermediate_frequency
    new_signals_data = _correlate_signals(
        sat.signals, per_signal_completed, dc, signal, sat.code_doppler, sat.code_phase,
        carrier_frequency, sat.carrier_phase, sampling_frequency, sat.signal_start_sample,
        samples_to_integrate, sat.prn, num_samples_signal,
    )
    update(sat, samples_to_integrate, intermediate_frequency, sampling_frequency, new_signals_data)
end

@inline function _dc_one_group!(g::SignalGroup, dc::_OneBitDC, measurements::BandMeasurements)
    vals = g.satellites.values
    isempty(vals) && return nothing
    m = measurements[band_key(g.band)]
    eltype(m.samples) === Complex{Int16} || throw(
        ArgumentError(
            string(
                "OneBitDownconvertAndCorrelator requires `Complex{Int16}` measurement ",
                "samples (12-bit ADC); got element type ", eltype(m.samples),
                ". Use a CPU(Threaded)DownconvertAndCorrelator for floating-point samples.",
            ),
        ),
    )
    _dc_group_loop!(
        dc, vals, m.samples, get_num_samples(m), m.sampling_frequency, m.intermediate_frequency,
    )
end

@inline function _dc_group_loop!(dc::OneBitDownconvertAndCorrelator, vals, args::Vararg{Any,4})
    @inbounds for i in eachindex(vals)
        vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

@inline function _dc_group_loop!(
    dc::OneBitThreadedDownconvertAndCorrelator,
    vals,
    args::Vararg{Any,4},
)
    @batch for i = 1:length(vals)
        @inbounds vals[i] = _update_tracked_sat_correlator(vals[i], dc, args...)
    end
    return nothing
end

"""
$(SIGNATURES)

Downconvert and correlate all satellites with the one-bit bit-wise backend.
"""
function downconvert_and_correlate(
    dc::_OneBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    new_track_state =
        TrackState(track_state; groups = _copy_groups_slot_vectors(track_state.groups))
    downconvert_and_correlate!(dc, measurements, new_track_state)
end

"""
$(SIGNATURES)

In-place one-bit downconvert and correlate. Returns the same `track_state`.
"""
function downconvert_and_correlate!(
    dc::_OneBitDC,
    measurements::BandMeasurements,
    track_state::TrackState,
)
    _foreach_group!(_dc_one_group!, track_state.groups, dc, measurements)
    return track_state
end
