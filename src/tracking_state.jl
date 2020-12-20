"""
$(SIGNATURES)

TrackingState that holds the tracking state.
"""
struct TrackingState{
        S <: AbstractGNSSSystem{},
        C <: AbstractCorrelator,
        CALF <: AbstractLoopFilter,
        COLF <: AbstractLoopFilter,
        CN <: AbstractCN0Estimator,
        V <: Union{StructArray, CuArray},
    }
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
    downconverted_signal::V
    carrier::V
    code::V
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
    GNSS::AbstractGNSSSystem{T},
    carrier_doppler,
    code_phase;
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(GNSS),
    carrier_phase = 0.0,
    carrier_loop_filter::CALF = ThirdOrderBilinearLF(),
    code_loop_filter::COLF = SecondOrderBilinearLF(),
    sc_bit_detector = SecondaryCodeOrBitDetector(),
    num_ants = NumAnts(1),
    correlator::C = get_default_correlator(GNSS, num_ants),
    integrated_samples = 0,
    prompt_accumulator = zero(ComplexF64),
    cn0_estimator::CN = MomentsCN0Estimator(20),
    num_samples = 0
) where {
    T <: Array,
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator
}
    if found(sc_bit_detector)
        code_phase = mod(code_phase, get_code_length(GNSS) *
            get_secondary_code_length(GNSS))
    else
        code_phase = mod(code_phase, get_code_length(GNSS))
    end
    downconverted_signal = init_downconverted_signal(num_ants, num_samples)
    carrier = StructArray{Complex{Int16}}(undef, 0)
    code = Vector{Int16}(undef, 0)

    TrackingState{C, CALF, COLF, CN, typeof(downconverted_signal)}(
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
        code,
        num_samples,
    )
end

# GPU tracking state constructor
function TrackingState(
    GNSS::S,
    carrier_doppler,
    code_phase;
    code_doppler = carrier_doppler * get_code_center_frequency_ratio(GNSS),
    carrier_phase = 0.0,
    carrier_loop_filter::CALF = ThirdOrderBilinearLF(),
    code_loop_filter::COLF = SecondOrderBilinearLF(),
    sc_bit_detector = SecondaryCodeOrBitDetector(),
    num_ants = NumAnts(1),
    correlator::C = get_default_correlator(GNSS, num_ants),
    integrated_samples = 0,
    prompt_accumulator = zero(ComplexF64),
    cn0_estimator::CN = MomentsCN0Estimator(20),
    num_samples = 0
) where {
    T <: CuArray,
    S <: AbstractGNSSSystem{T},
    C <: AbstractCorrelator,
    CALF <: AbstractLoopFilter,
    COLF <: AbstractLoopFilter,
    CN <: AbstractCN0Estimator
}
    if found(sc_bit_detector)
        code_phase = mod(code_phase, get_code_length(GNSS) *
            get_secondary_code_length(GNSS))
    else
        code_phase = mod(code_phase, get_code_length(GNSS))
    end
    downconverted_signal = CuArray{Complex{Float32}}(undef, num_samples)
    carrier = CuArray{Complex{Float32}}(undef, num_samples)
    code = CuArray{Complex{Float32}}(undef, num_samples)

    TrackingState{S, C, CALF, COLF, CN, typeof(downconverted_signal)}(
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

# One dimensional signal initialization
function init_downconverted_signal(num_ants::NumAnts{1}, num_samples::Int)
    StructArray{Complex{Int16}}(undef, num_samples)
end

# N-dimensional signal initialization
function init_downconverted_signal(num_ants::NumAnts{N}, num_samples::Int) where N
    StructArray{Complex{Int16}}(undef, num_samples, N)
end

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
