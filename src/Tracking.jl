module Tracking
    using DocStringExtensions, GNSSSignals, LinearAlgebra, Statistics, StaticArrays, Unitful
    import Unitful: Hz, s, ms

    export prompt,
        init_tracking,
        Initials

    struct TrackingResults{P}
        carrier_doppler::typeof(1.0Hz)
        carrier_phase::Float64
        code_doppler::typeof(1.0Hz)
        code_phase::Float64
        prompt::P
        data_bits::UInt
        num_bits::Int
        num_processed_prns::Int
        cn0::Float64
    end

    struct CodeShift{N}
        samples::Int
        actual_shift::Float64
    end

    struct Phases
        carrier::Float64
        code::Float64
    end

    struct Dopplers
        carrier::typeof(1.0Hz)
        code::typeof(1.0Hz)
    end

    struct Initials
        carrier_doppler::typeof(1.0Hz)
        carrier_phase::Float64
        code_doppler::typeof(1.0Hz)
        code_phase::Float64
    end

    struct DataBits{T<:AbstractGNSSSystem}
        synchronisation_buffer::UInt
        num_bits_in_synchronisation_buffer::Int
        first_found_after_num_prns::Int
        prompt_accumulator::Float64
        buffer::UInt
        num_bits_in_buffer::Int
    end

    function Initials(res::TrackingResults)
        Initials(res.carrier_doppler, res.carrier_phase, res.code_doppler, res.code_phase)
    end

    function Initials(carrier_doppler, code_phase)
        Initials(carrier_doppler, 0.0, 0.0Hz, code_phase)
    end

    function Initials(system::AbstractGNSSSystem, carrier_doppler, code_phase)
        Initials(carrier_doppler, 0.0, carrier_doppler * (system.code_freq / system.center_freq), code_phase)
    end

    function CodeShift{N}(system::AbstractGNSSSystem, sample_freq, preferred_code_shift) where N
        sample_shift = round(preferred_code_shift * sample_freq / system.code_freq)
        actual_shift = sample_shift * system.code_freq / sample_freq
        CodeShift{N}(sample_shift, actual_shift)
    end

    function Phases(inits)
        Phases(inits.carrier_phase, inits.code_phase)
    end

    function Dopplers(inits)
        Dopplers(inits.carrier_doppler, inits.code_doppler)
    end

    function TrackingResults(dopplers::Dopplers, phases::Phases, correlator_outputs, data_bits::DataBits, num_integrated_prns, cn0)
        TrackingResults(dopplers.carrier, phases.carrier, dopplers.code, phases.code, correlator_outputs, data_bits.buffer, data_bits.num_bits_in_buffer, num_integrated_prns, cn0)
    end

    include("discriminators.jl")
    include("loop_filters.jl")
    include("data_bits.jl")
    include("cn0_estimation.jl")
    include("tracking_loop.jl")
    include("gpsl1.jl")
    include("gpsl5.jl")
end
