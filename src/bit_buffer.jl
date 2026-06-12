"""
    SyncResult

Outcome of a per-signal bit-sync / secondary-code-sync detector call.

Fields:

- `found::Bool` — whether the detector locked on this update.
- `phase::Int` — when `found = true`, the secondary-code chip the
  *upcoming* integration aligns to, in `0:secondary_code_length-1`
  (recovered by the rotation search in [`_secondary_code_search`](@ref)).
  Zero for signals without a secondary code (the L1 C/A bit-edge case
  fires at the data-bit boundary, where the upcoming integration starts a
  new bit, not at a secondary-chip offset).
- `polarity::Int8` — `+1` or `-1`; which match orientation the detector
  locked. Carries through to the post-sync prompt accumulator so that a
  negative-polarity lock doesn't trip the downstream bit decoder.
"""
struct SyncResult
    found::Bool
    phase::Int
    polarity::Int8
end

"""
$(SIGNATURES)

Fixed-alignment Hamming-tolerance template matcher for the GPS L1 C/A
bit-edge detector (`detect_bit_or_secondary_code_sync(::GPSL1CA, …)`).
Signals with a periodic secondary code (GPS L5I, GPS L1C-P) instead use
[`_secondary_code_search`](@ref), which sweeps all rotations.

Compares `code_block_bits & mask` against `template` for the "positive
polarity" hit and against `template ⊻ mask` for the "negative polarity"
hit. The first orientation whose Hamming distance does not exceed
`max_errors` wins; if neither fits within tolerance the result reports
`found = false`. `phase` is always 0: the bit-edge template only matches
at the data-bit boundary, so the upcoming integration starts a new bit
(there is no secondary-chip offset to recover).

Inlined so the per-signal template / mask / tolerance constants fold at
the call site.
"""
@inline function _try_match(
    code_block_bits::B,
    template::B,
    mask::B,
    max_errors::Int,
) where {B<:Unsigned}
    masked = code_block_bits & mask
    dist_pos = count_ones(masked ⊻ template)
    dist_pos <= max_errors && return SyncResult(true, 0, Int8(+1))
    dist_neg = count_ones(masked ⊻ (template ⊻ mask))
    dist_neg <= max_errors && return SyncResult(true, 0, Int8(-1))
    return SyncResult(false, 0, Int8(0))
end

"""
$(SIGNATURES)

Generic secondary-code sync detector shared by every signal that locks
onto a periodic secondary / overlay code — GPS L5I's NH10 and GPS L1C-P's
1800-chip overlay both route through here. Unlike [`_try_match`](@ref)
(which matches one fixed template alignment and is used for the GPS L1 C/A
*bit-edge*), this runs a full rotation search, so it locks after a single
secondary-code period in the worst case instead of two and recovers the
true secondary-code phase.

`received` is the prompt-sign sliding window with the newest block in
bit 0. `reference` is the secondary code packed in that **same**
newest-first order — bit `i` holds secondary chip `(N - 1 - i)` — so that
when the most recent `N` blocks span exactly one period ending on its
last chip, `received & mask == reference`. `N` is the secondary-code
length.

The search rotates the low `N` bits of `received` left by `d ∈ 0:N-1`,
tracking the best positive- and negated-polarity Hamming match in one
pass. The winning rotation `d` is how far the buffer leads the reference,
which maps to the secondary-chip offset of the **upcoming** integration
as `phase = mod(N - d, N)` — exactly the value the post-sync `code_phase`
snap ([`_snap_code_phase_from_synced_signal`](@ref)) anchors on. Returns
`SyncResult(false, 0, 0)` when the best distance exceeds `max_errors`.

Inlined so the per-signal `reference` / `N` constants fold at the call
site (the L5I window is 10, the L1C-P window 1800).
"""
@inline function _secondary_code_search(
    received::B,
    reference::B,
    secondary_code_length::Int,
    max_errors::Int,
) where {B<:Unsigned}
    N = secondary_code_length
    # `reference` lives in the low `N` bits; mask the search to that window.
    # `one(B) << N` is undefined when `N` equals the full bit width
    # (e.g. UInt1800 for L1C-P), so special-case the exact-width buffer.
    mask = N == 8 * sizeof(B) ? ~zero(B) : (one(B) << N) - one(B)
    masked = received & mask
    best_d = 0
    best_dist = N + 1
    best_pol = Int8(0)
    @inbounds for d in 0:(N-1)
        # Rotate-left by `d` within the N-bit window. `d == 0` is special-
        # cased because `masked >> N` is undefined for an exact-width buffer.
        shifted = d == 0 ? masked : ((masked << d) | (masked >> (N - d))) & mask
        dist_pos = count_ones(shifted ⊻ reference)
        # Both operands occupy only the low `N` bits, so the negated-polarity
        # distance is the complement within that window.
        dist_neg = N - dist_pos
        if dist_pos < best_dist
            best_dist = dist_pos
            best_d = d
            best_pol = Int8(+1)
        end
        if dist_neg < best_dist
            best_dist = dist_neg
            best_d = d
            best_pol = Int8(-1)
        end
    end
    best_dist > max_errors && return SyncResult(false, 0, Int8(0))
    SyncResult(true, mod(N - best_d, N), best_pol)
end

"""
$(SIGNATURES)

BitBuffer to buffer bits.

The `code_block_buffer` field is the sync-search sliding window — its
width `B` is chosen per signal by [`get_code_block_buffer_type`](@ref) so
that a single integer can hold the entire pre-sync search horizon (one
NH10 period for GPS L5I, 40 primary blocks for GPS L1 C/A, 1800 chips
for the GPS L1C-P overlay, etc.). After sync the field is dead state and
the decoded navigation bits accumulate in the fixed-width
`buffer::UInt128` instead.
"""
struct BitBuffer{B<:Unsigned}
    code_block_buffer::B
    code_block_buffer_lengh::Int
    found::Bool
    secondary_phase::Int      # 0 until found; secondary-chip offset post-sync
    polarity::Int8            # +1 or -1 once found; 0 before sync
    buffer::UInt128
    length::Int
    prompt_accumulator::ComplexF64
    prompt_accumulator_integrated_code_blocks::Int
    soft_bits::Vector{Float32}
end

# Default constructor preserves the pre-refactor `UInt128`-backed search
# buffer. Once `get_code_block_buffer_type` lands (Step 2) the per-signal
# `TrackedSignal` constructor picks the right width instead.
function BitBuffer()
    BitBuffer{UInt128}(
        zero(UInt128), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0, Float32[],
    )
end

# Typed empty constructor used by the per-signal `TrackedSignal` path.
function BitBuffer{B}() where {B<:Unsigned}
    BitBuffer{B}(
        zero(B), 0, false, 0, Int8(0), zero(UInt128), 0, complex(0.0, 0.0), 0, Float32[],
    )
end

# Convenience outer constructor matching the legacy 7-arg form (no phase /
# polarity arguments — assumed zero; soft bits default to empty). Used by test
# code that builds a `BitBuffer` from raw integer / Complex{Int} literals.
function BitBuffer(
    code_block_buffer::B,
    code_block_buffer_lengh::Integer,
    found::Bool,
    buffer::Integer,
    length::Integer,
    prompt_accumulator::Complex,
    prompt_accumulator_integrated_code_blocks::Integer,
) where {B<:Unsigned}
    BitBuffer{B}(
        code_block_buffer,
        Int(code_block_buffer_lengh),
        found,
        0,
        Int8(0),
        UInt128(buffer),
        Int(length),
        ComplexF64(prompt_accumulator),
        Int(prompt_accumulator_integrated_code_blocks),
        Float32[],
    )
end

@inline get_bits(bit_buffer::BitBuffer) = bit_buffer.buffer
@inline length(bit_buffer::BitBuffer) = bit_buffer.length
@inline has_bit_or_secondary_code_been_found(bit_buffer::BitBuffer) = bit_buffer.found

# Get the soft bits, i.e. the accumulated (summed) filtered prompt of each
# completed bit. The sign of each soft bit corresponds to the respective hard
# bit returned by `get_bits`. The buffer is reset to length 0 at the start of
# each `track` call, mirroring the hard bit buffer. Kept as a plain comment (not
# a docstring) to match the sibling accessors `get_bits` / `get_num_bits`, which
# `checkdocs = :exports` would otherwise require to appear in the manual.
@inline get_soft_bits(bit_buffer::BitBuffer) = bit_buffer.soft_bits

"""
$(SIGNATURES)

Width of the sync-search sliding window for `signal`, returned as a
concrete `Unsigned` subtype.

Each signal's bit/secondary-code sync detector matches the buffered
primary-code-block signs against a per-signal template. This trait
picks the smallest integer width that holds one full search horizon for
that signal:

| Signal            | Returns   | Window |
|-------------------|-----------|--------|
| GPS L1 C/A        | `UInt64`  | 40 blocks (2 × 20 blocks/symbol) |
| Galileo E1B       | `UInt8`   | 8 blocks (2 × 4 blocks/symbol)   |
| GPS L5I           | `UInt32`  | 20 blocks (2 × NH10 length)      |
| GPS L1C-D         | `UInt8`   | n/a — symbol = primary period, buffer unused |
| GPS L1C-P         | `UInt1800`| 1800-chip overlay (added in step 4) |

The default for any signal not specialized below is `UInt64`. The width
flows through `BitBuffer{B}` and `TrackedSignal{Sig, B, C, PCF}` so the
parameter chain stays type-stable at construction.
"""
@inline get_code_block_buffer_type(::AbstractGNSSSignal) = UInt64

"""
$(SIGNATURES)

Per-signal Hamming tolerance used by the bit-edge / secondary-code
sync-search detectors, expressed as a fraction of the search window.

Returns the largest **fraction** of bit-flips the per-signal
`detect_bit_or_secondary_code_sync` accepts before reporting
`found = true`. Each detector converts this to an integer error budget
at its call site: `max_errors = floor(Int, tolerance × window_size)`.

Default is `0.025` (2.5 %), which discretizes per-signal as:

| Signal      | Window (blocks) | Effective `max_errors` |
|-------------|-----------------|------------------------|
| GPS L1 C/A  | 40              | 1                      |
| GPS L5I     | 10              | 0 (exact match)        |
| GPS L1C-P   | 1800            | 45                     |
| Galileo E1B | n/a — trivial   | unused                 |
| GPS L1C-D   | n/a — trivial   | unused                 |

Galileo E1B and GPS L1C-D broadcast one channel symbol per primary
code period, so their detectors return `SyncResult(true, 0, +1)`
unconditionally and never invoke [`_try_match`](@ref) — the trait
default applies but the value is ignored.

# Overriding

To loosen the tolerance for low-C/N₀ work, dispatch the trait on the
signal type in your own module:

```julia
Tracking.get_bit_edge_or_secondary_code_tolerance(::GPSL1CA) = 0.05
```

The override takes effect at the next call to
`detect_bit_or_secondary_code_sync` — there is no need to rebuild any
TrackState. The trait is `@inline`'d so the override folds at the
detector's call site.
"""
@inline get_bit_edge_or_secondary_code_tolerance(::AbstractGNSSSignal) = 0.025

# Number of primary-code blocks that form one navigation bit.
# Returns 0 for pilot signals (`data_frequency = 0`), where the concept is
# undefined; callers must guard for that case before using the result.
@inline function _calc_num_code_blocks_that_form_a_bit(signal::AbstractGNSSSignal)
    data_freq = get_data_frequency(signal)
    iszero(data_freq) && return 0
    Int(get_code_frequency(signal) / (get_code_length(signal) * data_freq))
end

"""
$(SIGNATURES)

Buffer data bits based on the prompt accumulation and the current prompt value.
"""
function buffer(signal::AbstractGNSSSignal, prn::Integer, bit_buffer::BitBuffer{B}, integrated_code_blocks, prompt) where {B<:Unsigned}
    # The divide is deferred to the helper — pilot signals
    # (`get_data_frequency = 0`) would otherwise blow up here with `Int(Inf)`.
    num_code_blocks_that_form_a_bit = _calc_num_code_blocks_that_form_a_bit(signal)

    if (bit_buffer.found == false)
        return _buffer_find_bit(
            signal, prn, bit_buffer, num_code_blocks_that_form_a_bit,
            integrated_code_blocks, prompt,
        )
    end

    # Pilot signals (e.g. GPS L1C-P) carry no navigation data; once their
    # secondary code is synced there is nothing to decode, so the buffer is
    # left untouched (the `found` / `secondary_phase` / `polarity` state stays
    # for code-phase anchoring and longer coherent integration).
    num_code_blocks_that_form_a_bit == 0 && return bit_buffer

    prompt_accumulator = bit_buffer.prompt_accumulator + prompt
    prompt_accumulator_integrated_code_blocks =
        bit_buffer.prompt_accumulator_integrated_code_blocks + integrated_code_blocks

    if prompt_accumulator_integrated_code_blocks == num_code_blocks_that_form_a_bit
        # Flip the decoded bit if the detector locked at negative polarity:
        # the prompt accumulator's real-part sign is then inverted relative
        # to the data symbol's "0/1" convention.
        bit_acc = bit_buffer.polarity < 0 ?
                  -real(prompt_accumulator) : real(prompt_accumulator)
        bit = bit_acc > 0
        # Store the polarity-corrected accumulation so the soft bit's sign
        # matches the decoded hard bit.
        push!(bit_buffer.soft_bits, Float32(bit_acc))
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            get_bits(bit_buffer) << 1 + UInt64(bit),
            length(bit_buffer) + 1,
            zero(prompt_accumulator),
            0,
            bit_buffer.soft_bits,
        )
    else
        return BitBuffer{B}(
            bit_buffer.code_block_buffer,
            bit_buffer.code_block_buffer_lengh,
            true,
            bit_buffer.secondary_phase,
            bit_buffer.polarity,
            bit_buffer.buffer,
            bit_buffer.length,
            prompt_accumulator,
            prompt_accumulator_integrated_code_blocks,
            bit_buffer.soft_bits,
        )
    end
end

function _buffer_find_bit(signal, prn::Integer, bit_buffer::BitBuffer{B}, num_code_blocks_that_form_a_bit,
                          integrated_code_blocks, prompt) where {B<:Unsigned}
    if (integrated_code_blocks != 1)
        error(
            "The number code blocks must be equal to 1 if bit or secondary code hasn't been found yet.",
        )
    end
    code_block_buffer = bit_buffer.code_block_buffer << 1 + B(real(prompt) > 0)
    code_block_buffer_lengh = bit_buffer.code_block_buffer_lengh + 1
    sync = detect_bit_or_secondary_code_sync(
        signal,
        prn,
        code_block_buffer,
        code_block_buffer_lengh,
    )
    if !sync.found
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_lengh,
            false,
            0,
            Int8(0),
            zero(UInt128),
            0,
            complex(0.0, 0.0),
            0,
            bit_buffer.soft_bits,
        )
    end
    if get_secondary_code_length(signal) > 1
        # Secondary-code signals (GPS L5I, GPS L1C-P): the buffered pre-sync
        # prompt signs are modulated by the secondary code, not the navigation
        # data, so there are no data bits to recover from them. Data-bit
        # decoding starts fresh post-sync. The rotation search locks at *any*
        # secondary chip, so the upcoming integration starts `sync.phase` blocks
        # into the current data bit — not at its boundary. Seed the accumulator
        # block count with `sync.phase` so the first emitted bit completes
        # exactly `num_code_blocks_that_form_a_bit − sync.phase` blocks later, on
        # the data-bit boundary, instead of `num_code_blocks_that_form_a_bit`
        # blocks after the lock instant (issue #125). For pilots
        # (`num_code_blocks_that_form_a_bit == 0`) the post-sync `buffer` path
        # returns early and never consults this seed, so it is harmless there.
        return BitBuffer{B}(
            code_block_buffer,
            code_block_buffer_lengh,
            true,
            sync.phase,
            sync.polarity,
            zero(UInt128),
            0,
            complex(0.0, 0.0),
            sync.phase,
            bit_buffer.soft_bits,
        )
    end
    num_bits = min(
        div(code_block_buffer_lengh, num_code_blocks_that_form_a_bit),
        div(sizeof(code_block_buffer) * 8, num_code_blocks_that_form_a_bit),
    )
    bits = reduce(num_bits:-1:1; init = UInt128(0)) do bits, bit_index
        # Don't shadow the outer `integrated_code_blocks` argument here:
        # because this closure reassigns its own local of the same name,
        # Julia conservatively boxes the outer parameter even though this
        # closure path does not always run, leading to a per-call
        # Core.Box allocation.
        bit_sum = sum(0:num_code_blocks_that_form_a_bit-1) do code_block_index
            buffer_code_block_index =
                (bit_index - 1) * num_code_blocks_that_form_a_bit + code_block_index
            ((code_block_buffer & (one(B) << buffer_code_block_index)) > 0) * 2 - 1
        end
        push!(bit_buffer.soft_bits, Float32(bit_sum))
        bits << 1 + (bit_sum > 0)
    end
    return BitBuffer{B}(
        code_block_buffer,
        code_block_buffer_lengh,
        true,
        sync.phase,
        sync.polarity,
        bits,
        num_bits,
        complex(0, 0),
        0,
        bit_buffer.soft_bits,
    )
end

function reset(bit_buffer::BitBuffer{B}) where {B<:Unsigned}
    empty!(bit_buffer.soft_bits)
    BitBuffer{B}(
        bit_buffer.code_block_buffer,
        bit_buffer.code_block_buffer_lengh,
        bit_buffer.found,
        bit_buffer.secondary_phase,
        bit_buffer.polarity,
        zero(UInt128),
        0,
        bit_buffer.prompt_accumulator,
        bit_buffer.prompt_accumulator_integrated_code_blocks,
        bit_buffer.soft_bits,
    )
end
