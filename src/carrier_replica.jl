"""
$(SIGNATURES)

Updates the carrier phase.
"""
function update_carrier_phase(
    num_samples,
    carrier_frequency,
    sampling_frequency,
    start_phase,
)
    phase = num_samples * carrier_frequency / sampling_frequency + start_phase
    mod(phase + 0.5, 1) - 0.5
end
