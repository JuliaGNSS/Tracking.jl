"""
$(SIGNATURES)

TrackingResults that hold the correlation results.
"""
struct TrackingResults{
        TS <: TrackingState,
        C <: AbstractCorrelator,
        CS,
        ELI
    }
    state::TS
    correlator::C
    correlator_sample_shifts::CS
    early_late_index_shift::ELI
    correlator_carrier_frequency::typeof(1.0Hz)
    correlator_carrier_phase::Float64
    got_correlator::Bool
    bit_buffer::BitBuffer
    cn0::typeof(1.0dBHz)
end

"""
$(SIGNATURES)

Get state of the tracking result.
"""
@inline get_state(results::TrackingResults) = results.state

"""
$(SIGNATURES)

Get the carrier doppler of the tracking result.
"""
@inline get_carrier_doppler(results::TrackingResults) = get_carrier_doppler(results.state)

"""
$(SIGNATURES)

Get the carrier phase of the tracking result.
"""
@inline get_carrier_phase(results::TrackingResults) = get_carrier_phase(results.state) * 2π

"""
$(SIGNATURES)

Get the code doppler of the tracking result.
"""
@inline get_code_doppler(results::TrackingResults) = get_code_doppler(results.state)

"""
$(SIGNATURES)

Get the code phase of the tracking result.
"""
@inline get_code_phase(results::TrackingResults) = get_code_phase(results.state)

"""
$(SIGNATURES)

Get the correlator of the tracking result.
"""
@inline get_correlator(results::TrackingResults) = results.correlator

"""
$(SIGNATURES)

Get the correlator sample shifts of the tracking result.
"""
@inline get_correlator_sample_shifts(results::TrackingResults) = results.correlator_sample_shifts

"""
$(SIGNATURES)

Get the correlator sample shifts of the tracking result.
"""
@inline get_early_late_index_shift(results::TrackingResults) = results.early_late_index_shift

"""
$(SIGNATURES)

Get carrier phase at the same timestamp when the correlator is created of the tracking
result.
"""
@inline function get_correlator_carrier_phase(results::TrackingResults)
    results.correlator_carrier_phase * 2π
end

"""
$(SIGNATURES)

Get carrier frequency at the same timestamp when the correlator is created of the tracking
result.
"""
@inline function get_correlator_carrier_frequency(results::TrackingResults)
    results.correlator_carrier_frequency
end

"""
$(SIGNATURES)

Get the early of the tracking result.
"""
@inline get_early(results::TrackingResults) = get_early(
    get_correlator(results),
    get_correlator_sample_shifts(results),
    get_early_late_index_shift(results)
)

"""
$(SIGNATURES)

Get the prompt of the tracking result.
"""
@inline get_prompt(results::TrackingResults) = get_prompt(get_correlator(results), get_correlator_sample_shifts(results))

"""
$(SIGNATURES)

Get the late of the tracking result.
"""
@inline get_late(results::TrackingResults) = get_late(
    get_correlator(results),
    get_correlator_sample_shifts(results),
    get_early_late_index_shift(results)
)

"""
$(SIGNATURES)

Get the bits of the tracking result.
"""
@inline get_bits(results::TrackingResults) = get_bits(results.bit_buffer)

"""
$(SIGNATURES)

Get the number of bits of the tracking result.
"""
@inline get_num_bits(results::TrackingResults) = length(results.bit_buffer)

"""
$(SIGNATURES)

Get the carrier to noise density ratio of the tracking result.
"""
@inline get_cn0(results::TrackingResults) = results.cn0

"""
$(SIGNATURES)

Check if the secondary code or bit has been found.
"""
@inline get_secondary_code_or_bit_found(results::TrackingResults) =
    found(get_sc_bit_detector(results.state))
