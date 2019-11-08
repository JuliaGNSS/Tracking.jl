"""
$(SIGNATURES)

SecondaryCodeOrBitDetector to detect the secondary code or bits shift.
"""
struct SecondaryCodeOrBitDetector
    buffer::UInt64
    length::Int
    found::Bool
end

function SecondaryCodeOrBitDetector()
    SecondaryCodeOrBitDetector(0, 0, false)
end

@inline get_buffer(sc_bit_detector::SecondaryCodeOrBitDetector) = sc_bit_detector.buffer
@inline found(sc_bit_detector::SecondaryCodeOrBitDetector) = sc_bit_detector.found
@inline length(sc_bit_detector::SecondaryCodeOrBitDetector) = sc_bit_detector.length

"""
$(SIGNATURES)

Finds the secondary code or the bit shift.
"""
function find(::Type{S}, detector, prompt_correlator) where S <: AbstractGNSSSystem
    found(detector) && return detector
    prns = get_buffer(detector) << 1 + UInt64(real(prompt_correlator) > 0)
    num_prns = length(detector) + 1
    bit_found = is_upcoming_integration_new_bit(S, prns, num_prns)
    SecondaryCodeOrBitDetector(prns, num_prns, bit_found)
end
