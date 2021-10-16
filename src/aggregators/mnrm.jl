# A modified next reaction method 
struct MNRM <: AbstractAggregatorAlgorithm end
mutable struct MNRMJumpAggregation{T,S,F1,F2,RNG} <: AbstractSSAJumpAggregator
  next_jump::Int
  prev_jump::Int
  next_jump_time::T
  end_time::T
  react_rates::Vector{T}
  sum_rate::T
  ma_jumps::S
  rates::F1
  affects!::F2
  save_positions::Tuple{Bool,Bool}
  rng::RNG
  internal_waitingtimes::Vector{T}
end
MNRMJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; internal_waitingtimes, kwargs...) where {T,S,F1,F2,RNG} =
MNRMJumpAggregation{T,S,F1,F2,RNG}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, internal_waitingtimes)




  # creating the JumpAggregation structure (tuple-based constant jumps)
function aggregate(aggregator::MNRM, u, p, t, end_time, constant_jumps, ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using tuples
    rates, affects! = get_jump_info_tuples(constant_jumps)
    num_reactions = get_num_majumps(ma_jumps) + length(constant_jumps)
    internal_waitingtimes = zeros(typeof(t),num_reactions)
    build_jump_aggregation(MNRMJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; internal_waitingtimes = internal_waitingtimes, kwargs...)
end

function initialize!(p::MNRMJumpAggregation, integrator, u, params, t)
    num_reactions = get_num_majumps(p.ma_jumps) + length(p.rates)
    p.internal_waitingtimes = randexp(p.rng,num_reactions) # generate internal waiting time list for the initiation
    generate_jumps!(p, integrator, u, params, t)

    nothing
end

# calculate the next jump / jump time
function generate_jumps!(p::MNRMJumpAggregation, integrator, u, params, t)
    p.next_jump, ttnj = time_to_next_jump(p, u, params, t)
    @fastmath p.next_jump_time = t + ttnj
    update_internal_times!(p, ttnj, p.rng)
    nothing
end


@fastmath function time_to_next_jump(p::MNRMJumpAggregation{T,S,F1,F2,RNG}, u, params, t) where {T,S,F1 <: Tuple, F2 <: Tuple, RNG}
    
    @unpack internal_waitingtimes, react_rates = p
    fill_react_rates!(p, u, t)
    dt_list =@. copy(internal_waitingtimes/react_rates)
    ttnj, next_jump = findmin(dt_list)
    # sum_rate = 0. # sum_rate is not necessary in MNRM method
    next_jump, ttnj
end


# fill the propensity rates 
function fill_react_rates!(p::MNRMJumpAggregation,  u, t)
    react_rates = p.react_rates
  
    # mass action rates
    majumps   = p.ma_jumps
    num_majumps      = get_num_majumps(majumps)
    @inbounds for i in 1:num_majumps
        react_rates[i] = evalrxrate(u, i, majumps)
    end
    # constant jump rates
    rates = p.rates
    if !isempty(rates)
        idx = num_majumps + 1 
        @inbounds for rate in rates
            react_rates[idx] = rate(u, p, t)
            idx  += 1
        end
    end
end



function update_internal_times!(p::MNRMJumpAggregation,  ttnj, rng)
    p.internal_waitingtimes .-= ttnj*p.react_rates
    p.internal_waitingtimes[p.next_jump] = randexp(rng)
end


@inline function execute_jumps!(p::MNRMJumpAggregation, integrator, u, params, t)
    update_state!(p, integrator, u)
    nothing
end

# Change below

