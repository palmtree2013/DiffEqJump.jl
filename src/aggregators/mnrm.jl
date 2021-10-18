# A modified next reaction method 
struct MNRM <: AbstractAggregatorAlgorithm end
mutable struct MNRMJumpAggregation{T,S,F1,F2,RNG,DEPGR} <: AbstractSSAJumpAggregator
  next_jump::Int
  prev_jump::Int
  next_jump_time::T
  end_time::T
  cur_rates::Vector{T}
  sum_rate::T
  ma_jumps::S
  rates::F1
  affects!::F2
  save_positions::Tuple{Bool,Bool}
  rng::RNG
  internal_waitingtimes::Vector{T}
  dep_gr::DEPGR
end
function MNRMJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; num_specs, internal_waitingtimes, dep_graph=nothing, kwargs...) where {T,S,F1,F2,RNG}

    # a dependency graph is needed and must be provided if there are constant rate jumps
    if dep_graph === nothing
        if (get_num_majumps(maj) == 0) || !isempty(rs)
            error("To use ConstantRateJumps with the Next Reaction Method (NRM) algorithm a dependency graph must be supplied.")
        else
            dg = make_dependency_graph(num_specs, maj)
        end
    else
        dg = dep_graph

        # make sure each jump depends on itself
        add_self_dependencies!(dg)
    end
    # pq = MutableBinaryMinHeap{T}()
    MNRMJumpAggregation{T,S,F1,F2,RNG,typeof(dg)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, internal_waitingtimes, dg)
end

function aggregate(aggregator::MNRM, u, p, t, end_time, constant_jumps, ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using tuples
    rates, affects! = get_jump_info_tuples(constant_jumps)
    num_reactions = get_num_majumps(ma_jumps) + length(constant_jumps)
    internal_waitingtimes = zeros(typeof(t),num_reactions)
    build_jump_aggregation(MNRMJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; num_specs = length(u), internal_waitingtimes = internal_waitingtimes, kwargs...)
end

function initialize!(p::MNRMJumpAggregation, integrator, u, params, t)
    num_reactions = length(p.cur_rates)
    p.internal_waitingtimes = randexp(p.rng,num_reactions) # generate internal waiting time list for the initiation
    fill_rates_and_get_times!(p, u, params, t)
    nothing
end
@inline function execute_jumps!(p::MNRMJumpAggregation, integrator, u, params, t)
    update_state!(p, integrator, u)
    nothing
end
# calculate the next jump / jump time
function generate_jumps!(p::MNRMJumpAggregation, integrator, u, params, t)
    update_internal_times!(p, u, params, t)
    nothing
end



@inline function update_internal_times!(p::MNRMJumpAggregation, u, params, t)
    @unpack internal_waitingtimes, cur_rates, rates, ma_jumps, rng = p
    dep_rxs = p.dep_gr[copy(p.next_jump)]
    num_majumps = get_num_majumps(ma_jumps)
    @inbounds for rx in dep_rxs
        cur_rates[rx] = calculate_jump_rate(ma_jumps,num_majumps,rates,u,partial_path,t,rx)
    end
    pqdata = [internal_waitingtimes[i]/cur_rates[i] for i in 1:length(cur_rates)]
    # ttnj, next_jump = top_with_handle(MutableBinaryMinHeap(pqdata))
    ttnj, p.next_jump = findmin(pqdata)
    @fastmath p.next_jump_time = t + ttnj
    # update internal time
    @inbounds for rx in eachindex(cur_rates)
        internal_waitingtimes[rx] -= ttnj*cur_rates[rx]
    end
    internal_waitingtimes[p.next_jump] = randexp(rng)
    nothing
end


# fill the propensity rates 
function fill_rates_and_get_times!(p::MNRMJumpAggregation,  u, params,  t)
    @unpack internal_waitingtimes, cur_rates, rates, ma_jumps, rng = p
    # pqdata = [internal_waitingtimes[i]/cur_rates[i] for i in 1:length(cur_rates)]
    pqdata = Vector{typeof(t)}(undef,length(cur_rates))
    @inbounds for i in 1:get_num_majumps(ma_jumps)
        cur_rates[i] = evalrxrate(u, i, ma_jumps)
        pqdata[i] = internal_waitingtimes[i]/cur_rates[i]
    end
    # constant rates
    rates = p.rates
    idx   = get_num_majumps(ma_jumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        pqdata[idx] = internal_waitingtimes[i]/cur_rates[idx]
        idx += 1
    end
    ttnj, p.next_jump = findmin(pqdata)
    p.next_jump_time = t +ttnj
    @inbounds for rx in eachindex(cur_rates)
        internal_waitingtimes[rx] -= ttnj*cur_rates[rx]
    end
    internal_waitingtimes[p.next_jump] = randexp(rng)
end







# Change below

