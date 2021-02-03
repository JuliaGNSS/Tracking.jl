"""
$(SIGNATURES)

CarrierReplicaCPU that holds possible carrier representations.
The type is unknown, when TrackingState is initialized.
It is determined based on the signal, which is given in the track
function.
"""
struct CarrierReplicaCPU{
        CF32 <: StructArray{ComplexF32},
        CF64 <: StructArray{ComplexF64}
    }
    carrier_f32::CF32
    carrier_f64::CF64
end

function CarrierReplicaCPU()
    CarrierReplicaCPU(
        StructArray{Complex{Float32}}(undef, 0),
        StructArray{Complex{Float64}}(undef, 0)
    )
end

"""
$(SIGNATURES)

DownconvertedSignalCPU that holds possible downconverted signal representations.
The type is unknown, when TrackingState is initialized.
It is determined based on the signal, which is given in the track
function.
"""
struct DownconvertedSignalCPU{
        DSF32 <: StructArray{ComplexF32},
        DSF64 <: StructArray{ComplexF64}
    }
    downconverted_signal_f32::DSF32
    downconverted_signal_f64::DSF64
end

function DownconvertedSignalCPU(num_ants::NumAnts{1})
    DownconvertedSignalCPU(
        StructArray{Complex{Float32}}(undef, 0),
        StructArray{Complex{Float64}}(undef, 0)
    )
end

function DownconvertedSignalCPU(num_ants::NumAnts{N}) where N
    DownconvertedSignalCPU(
        StructArray{Complex{Float32}}(undef, 0),
        StructArray{Complex{Float64}}(undef, 0, N)
    )
end


"""
$(SIGNATURES)

TrackingState that holds the tracking state.
"""
struct TrackingState{
        S <: AbstractGNSS,
        C <: AbstractCorrelator,
        CALF <: AbstractLoopFilter,
        COLF <: AbstractLoopFilter,
        CN <: AbstractCN0Estimator,
        DS <: DownconvertedSignalCPU, # Union{DownconvertedSignalCPU, CuArray{ComplexF32}}
        CAR <: CarrierReplicaCPU, # Union{CarrierReplicaCPU, CuArray{ComplexF32}}
        COR <: Vector{Int8}, # Union{Vector{Int8}, CuArray{Float32}}
    }
    system::S
    init_carrier_doppler::typeof(1.0Hz)
    init_code_doppler::typeof(1.0Hz)
    carrier_doppler::typeof(1.0Hz)
    code_doppler::typeof(1.0Hz)
    carrier_phase::Float64
    code_phase::Float64
    correlator::C
    carrier_loop_filter::CALF
    code_loop_filter::COLF
    sc_bit_detector::SecondaryCodeOrBitDetector
    integrated_samples::Int
    prompt_accumulator::ComplexF64
    cn0_estimator::CN
    downconverted_signal::DS
    carrier::CAR
    code::COR
end

"""
$(SIGNATURES)

Convenient TrackingState constructor. Mandatory parameters are the GNSS system, the
carrier doppler `carrier_doppler` and the code phase `code_phase`. Optional parameters are
- Code doppler `code_doppler`, that defaults to carrier doppler multiplied with the code /
  center frequency ratio
- Carrier phase `carrier_phase`, that defaults to 0
- Carrier loop filter `carrier_loop_filter`, that defaults to the `ThirdOrderBilinarLF()`
- Code loop filter `code_loop_filter`, that defaults to the `SecondOrderBilinearLF()`
- Number of antenna elemants `num_ants`, that defaults to `NumAnts(1)`
- Correlator `correlator`, that defaults to `get_default_correlator(...)`
- Integrated samples `integrated_samples`, that defaults to 0
- Prompt accumulator for bit detection `prompt_accumulator`, that defaults to 0
- CN0 estimator `cn0_estimator`, that defaults to `MomentsCN0Estimator(20)`
"""
function TrackingState(
    system::S,
    carrier_doppler,
    code_phase;
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(system),
    carrier_phase = 0.0,
    carrier_loop_filter::CALF = ThirdOrderBilinearLF(),
    code_loop_filter::COLF = SecondOrderBilinearLF(),
    sc_bit_detector = SecondaryCodeOrBitDetector(),
    num_ants = NumAnts(1),
    correlator::C = get_default_correlator(system, num_ants),
    integrated_samples = 0,
    prompt_accumulator = zero(ComplexF64),
    cn0_estimator::CN = MomentsCN0Estimator(20)
) where {
    S <: AbstractGNSS,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator
}
    if found(sc_bit_detector)
        code_phase = mod(code_phase, get_code_length(system) *
            get_secondary_code_length(system))
    else
        code_phase = mod(code_phase, get_code_length(system))
    end
    downconverted_signal = DownconvertedSignalCPU(num_ants)
    carrier = CarrierReplicaCPU()
    code = Vector{Int8}(undef, 0)

    TrackingState{S, C, CALF, COLF, CN, typeof(downconverted_signal), typeof(carrier), typeof(code)}(
        system,
        carrier_doppler,
        code_doppler,
        carrier_doppler,
        code_doppler,
        carrier_phase / 2π,
        code_phase,
        correlator,
        carrier_loop_filter,
        code_loop_filter,
        sc_bit_detector,
        integrated_samples,
        prompt_accumulator,
        cn0_estimator,
        downconverted_signal,
        carrier,
        code
    )
end

@inline get_system(state::TrackingState) = state.system
@inline get_code_phase(state::TrackingState) = state.code_phase
@inline get_carrier_phase(state::TrackingState) = state.carrier_phase
@inline get_init_code_doppler(state::TrackingState) = state.init_code_doppler
@inline get_init_carrier_doppler(state::TrackingState) = state.init_carrier_doppler
@inline get_code_doppler(state::TrackingState) = state.code_doppler
@inline get_carrier_doppler(state::TrackingState) = state.carrier_doppler
@inline get_correlator(state::TrackingState) = state.correlator
@inline get_sc_bit_detector(state::TrackingState) = state.sc_bit_detector
@inline get_carrier_loop_filter(state::TrackingState) = state.carrier_loop_filter
@inline get_code_loop_filter(state::TrackingState) = state.code_loop_filter
@inline get_prompt_accumulator(state::TrackingState) = state.prompt_accumulator
@inline get_integrated_samples(state::TrackingState) = state.integrated_samples
@inline get_cn0_estimator(state::TrackingState) = state.cn0_estimator
@inline get_downconverted_signal(state::TrackingState) = state.downconverted_signal
@inline get_carrier(state::TrackingState) = state.carrier
@inline get_code(state::TrackingState) = state.code
