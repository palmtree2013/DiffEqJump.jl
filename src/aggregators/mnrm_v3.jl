# A modified next reaction method 
struct MNRM <: AbstractAggregatorAlgorithm end
mutable struct MNRMJumpAggregation{T,S,F1,F2,RNG,DEPGR,PQ} <: AbstractSSAJumpAggregator
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
  dep_gr::DEPGR
  pq::PQ
end
function MNRMJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; num_specs, dep_graph=nothing, kwargs...) where {T,S,F1,F2,RNG}

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
    # if T == Float64
        # pq = MutableBinaryMinHeap{Float64, DataStructures.FasterForward}()
    # else
        pq = MutableBinaryMinHeap{T}()
    # end
    MNRMJumpAggregation{T,S,F1,F2,RNG,typeof(dg), typeof(pq)}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng, dg, pq)
end

function aggregate(aggregator::MNRM, u, p, t, end_time, constant_jumps, ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using tuples
    rates, affects! = get_jump_info_fwrappers(u, p, t, constant_jumps)
    build_jump_aggregation(MNRMJumpAggregation, u, p, t, end_time, ma_jumps, rates, affects!, save_positions, rng; num_specs = length(u), kwargs...)
end

function initialize!(p::MNRMJumpAggregation, integrator, u, params, t)
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
    time_to_next_jump!(p, params, t)
    nothing
end

@inline function update_internal_times!(p::MNRMJumpAggregation, u, params, t)
    @unpack pq, cur_rates, rates, ma_jumps, rng = p
    dep_rxs = p.dep_gr[p.next_jump]
    num_majumps = get_num_majumps(ma_jumps)
    @inbounds for rx in dep_rxs
        oldrate = cur_rates[rx]
        @inbounds cur_rates[rx] = calculate_jump_rate(ma_jumps,num_majumps,rates,u,params,t,rx)
        if rx != p.next_jump && oldrate > zero(oldrate)
            if cur_rates[rx] > zero(eltype(cur_rates))
                update!(pq, rx, pq*oldrate/cur_rates[rx])
            else
                update!(pq, rx, typemax(t))
            end
        else
            if cur_rates[rx] > zero(eltype(cur_rates))
                update!(pq, rx, randexp(rng)/cur_rates[rx])
            else
                update!(pq, rx, typemax(t))
            end
        end
    end
    nothing
end

function time_to_next_jump!(p, params, t)
    ttnj, p.next_jump = top_with_handle(p.pq)
    @fastmath p.next_jump_time = t + ttnj
    # update absolute time
    for rx in eachindex(p.cur_rates)
        update!(p.pq, rx, p.pq[rx] - ttnj)
    end
end

# fill the propensity rates 
# reevaulate all rates, recalculate all jump times, and reinit the priority queue
function fill_rates_and_get_times!(p::MNRMJumpAggregation, u, params, t)

    # mass action jumps
    majumps   = p.ma_jumps
    cur_rates = p.cur_rates
    pqdata = Vector{typeof(t)}(undef,length(cur_rates))
    @inbounds for i in 1:get_num_majumps(majumps)
        cur_rates[i] = evalrxrate(u, i, majumps)
        pqdata[i] = randexp(p.rng) / cur_rates[i]
    end

    # constant rates
    rates = p.rates
    idx   = get_num_majumps(majumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        pqdata[idx] = randexp(p.rng) / cur_rates[idx]
        idx += 1
    end

    # setup a new indexed priority queue to storing rx times
    p.pq = MutableBinaryMinHeap(pqdata)
    # ttnj, p.next_jump = findmin(p.pq)
    # @fastmath p.next_jump_time = t + ttnj
    time_to_next_jump!(p, params, t)
    nothing
end






# Change below

