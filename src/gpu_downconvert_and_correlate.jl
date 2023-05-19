struct GPUBuffers{T, DS <: CuArray{Complex{T}, 3}}
    downconverted_and_decoded_signal::DS
end

function GPUBuffers(
    ::Type{T},
    system::AbstractGNSS,
    state::SatState,
    num_samples,
) where T
	GPUBuffers(
        CuArray{ComplexF32}(undef, (num_samples, get_num_ants(state), get_num_accumulators(get_correlator(state))))
    )
end

function GPUBuffers(
    system::AbstractGNSS,
    state::SatState,
    num_samples,
)
	GPUBuffers(Float32, system, state, num_samples)
end

struct GPUDownconvertAndCorrelator{N, CB <: TupleLike{<:NTuple{N, <:AbstractGNSS}}, B <: TupleLike{<:NTuple{N, Vector{<:GPUBuffers}}}} <: AbstractDownconvertAndCorrelator
	code_buffers::CB
	buffers::B
	function GPUDownconvertAndCorrelator(code_buffers::TupleLike{<:NTuple{N, <:AbstractGNSS}}, buffers::TupleLike{<:NTuple{N, Vector{<:GPUBuffers}}}) where N
        new{N, typeof(code_buffers), typeof(buffers)}(code_buffers, buffers)
    end
end

function GPUDownconvertAndCorrelator(::Type{T}, system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}}, num_samples) where {T, N}
	GPUDownconvertAndCorrelator(
		map(convert_code_to_texture_memory, map(get_system, system_sats_states)),
		map(system_sats_states) do system_sats_state
            system = system_sats_state.system
            map(system_sats_state.states) do state
                GPUBuffers(T, system, state, num_samples)
            end
        end
	)
end

function GPUDownconvertAndCorrelator(
    system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}},
    num_samples::Integer,
) where N
	GPUDownconvertAndCorrelator(Float32, system_sats_states, num_samples)
end

import Adapt

Adapt.@adapt_structure GPSL1
Adapt.@adapt_structure GPSL5
Adapt.@adapt_structure GalileoE1B

recreate_system_with_texture(system::GPSL1, texture) = GPSL1(texture)
recreate_system_with_texture(system::GPSL5, texture) = GPSL5(texture)
recreate_system_with_texture(system::GalileoE1B, texture) = GalileoE1B(texture)

function convert_code_to_texture_memory(system::S) where S <: AbstractGNSS
	# Get only base code without secondary code, since otherwise code might be too
	# large for texture memory. Texture memory has a max size of 65536 in each
	# 2D dimension. GPSL5 would have a size of 102300 with secondary code.
	# Without secondary code GPSL5 has a code size of 10230.
	# The secondary code is multiplied in the kernel instead.
	# The same goes for any subcarrier code.
	codes = get_codes(system)[1:get_code_length(system), :]
	recreate_system_with_texture(system, CuTexture(
		CuTextureArray(CuArray(Float32.(codes))),
		address_mode = CUDA.ADDRESS_MODE_WRAP,
		interpolation = CUDA.NearestNeighbour(),
	))
end

function downconvert_and_correlate(
	downconvert_and_correlator::GPUDownconvertAndCorrelator,
	signal,
	sampling_frequency,
	intermediate_frequency,
	system_sats_states::TupleLike{<:NTuple{N, SystemSatsState}},
	params::TupleLike{<:NTuple{N, Vector{SampleParams}}},
) where N
	map(params, system_sats_states, downconvert_and_correlator.code_buffers, downconvert_and_correlator.buffers) do system_params, system_sats, system, buffers
		map(system_params, system_sats.states, buffers) do sat_params, sat_state, buffer
			if sat_params.signal_samples_to_integrate == 0
				return sat_state.correlator
			end
			carrier_frequency = sat_state.carrier_doppler + intermediate_frequency
			code_frequency = sat_state.code_doppler + get_code_frequency(system)
			downconvert_and_correlate!(
				system,
				signal,
				sat_state.correlator,
				sat_state.code_phase,
				sat_state.carrier_phase,
				code_frequency,
				carrier_frequency,
				sampling_frequency,
				sat_params.signal_start_sample,
				sat_params.signal_samples_to_integrate,
				sat_state.prn,
				buffer.downconverted_and_decoded_signal,
			)::typeof(sat_state.correlator)
		end
	end
end

function get_code(system::AbstractGNSS, phase, prn)
	get_code(system, get_modulation(system), phase, prn)
end

function get_code(system::AbstractGNSS, modulation::GNSSSignals.LOC, phase, prn)
	# Must add 0.5 because CUDA uses nearest neighbour instead of floor.
	system.codes[phase + 0.5f0, prn] * get_secondary_code(system, phase)
end

function get_code(system::AbstractGNSS, modulation::GNSSSignals.BOC, phase, prn)
	# Must add 0.5 because CUDA uses nearest neighbour instead of floor.
	system.codes[phase + 0.5f0, prn] * get_secondary_code(system, phase) *
	GNSSSignals.get_subcarrier_code(modulation, phase)
end

# This kernel currently assumes that we have more threads than number of samples
# to precess
# TODO: handle more samples than number of threads available
function downconvert_and_decode_prn_kernel!(
	downconverted_and_decoded_signal,
	signal,
	system::AbstractGNSS,
	prn::Int32,
	correlator_sample_shifts,
	num_samples::Int32,
	code_frequency,
	carrier_frequency,
	sampling_frequency,
	start_code_phase::Float32,
	start_carrier_phase::Float32,
	start_sample::Int32,
	num_ants::NumAnts{N},
) where {N}
	sample = ((blockIdx().x - 0x1) * blockDim().x + (threadIdx().x - 0x1))
	index = sample + 0x1
	if sample < num_samples
		carrier_wipe_off = cis(-Float32(2π) * (sample * carrier_frequency / sampling_frequency + start_carrier_phase))
		for sample_shift_index in eachindex(correlator_sample_shifts)
			sample_shift = correlator_sample_shifts[sample_shift_index]
			code = get_code(system, (sample + sample_shift) * code_frequency / sampling_frequency + start_code_phase, prn)
			for antenna_index = 0x1:N
				@inbounds downconverted_and_decoded_signal[index, antenna_index, sample_shift_index] = signal[sample + start_sample, antenna_index] * carrier_wipe_off * code
			end
		end
	end
    return
end

function downconvert_and_correlate!(
    code_buffer,
    signal,
    correlator::AbstractCorrelator{M},
    code_phase,
    carrier_phase,
    code_frequency,
    carrier_frequency,
    sampling_frequency,
    signal_start_sample,
    num_samples_left,
    prn,
	downconverted_and_decoded_signal,
) where M
	# Assume 1024 to be the max number of threads
	# TODO: Evaluate this at run time
	threads = min(size(signal, 1), 1024)
	blocks = cld(size(signal, 1), threads)
	num_correlators = size(downconverted_and_decoded_signal, 3)
	@cuda threads=threads blocks=blocks downconvert_and_decode_prn_kernel!(
		downconverted_and_decoded_signal,
		signal,
		code_buffer,
		Int32(prn),
		correlator.shifts,
		Int32(num_samples_left),
		Float32(code_frequency / Hz),
		Float32(carrier_frequency / Hz),
		Float32(sampling_frequency / Hz),
		Float32(code_phase),
		Float32(carrier_phase),
		Int32(signal_start_sample),
		NumAnts{M}()
	)
	correlated_signal = sum(view(downconverted_and_decoded_signal, 1:num_samples_left, :, :), dims = 1)
	result = reshape(Array(correlated_signal), M, num_correlators)
	gpu_add_to_accumulators(correlator, result)
end

function gpu_add_to_accumulators(correlator::AbstractCorrelator{1}, result)
	update_accumulator(correlator, SVector(map((a, b) -> a + b[1], get_accumulators(correlator), eachcol(result))))
end

function gpu_add_to_accumulators(correlator::AbstractCorrelator{M}, result) where M
	update_accumulator(correlator, map(+, get_accumulators(correlator), eachcol(result)))
end