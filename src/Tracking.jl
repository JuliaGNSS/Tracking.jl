module Tracking
    using DocStringExtensions, GNSSSignals, DataStructures, LinearAlgebra, Statistics, StaticArrays, FFTW, Unitful
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
        databitbuffer::UInt
        num_databits::UInt
        num_integrationbits::UInt
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

    struct DataBits{T<:AbstractGNSSSystem}
        synchronisation_buffer::UInt
        num_bits_in_synchronisation_buffer::UInt
        first_found_after_num_prns::Int
        prompt_accumulator::Float64
        buffer::UInt
        num_bits_in_buffer::UInt
    end

    struct Initials
        carrier_doppler::typeof(1.0Hz)
        carrier_phase::Float64
        code_doppler::typeof(1.0Hz)
        code_phase::Float64
    end

    function Initials(res::TrackingResults)
        Initials(res.carrier_doppler, res.carrier_phase, res.code_doppler, res.code_phase)
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

    function TrackingResults(dopplers::Dopplers, phases::Phases, correlator_outputs, data_bits::DataBits, num_integrated_prns)
        TrackingResults(dopplers.carrier, phases.carrier, dopplers.code, phases.code, correlator_outputs, data_bits.buffer, data_bits.num_bits_in_buffer, num_integrated_prns)
    end

    function DataBits(system::T) where T <: AbstractGNSSSystem
        DataBits{T}(0, 0, -1, 0, 0, 0)
    end

    include("discriminators.jl")
    include("loop_filters.jl")
    include("gpsl1.jl")
    include("gpsl5.jl")
    include("tracking_loop.jl")
end
