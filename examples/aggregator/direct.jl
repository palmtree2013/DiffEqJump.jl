mutable struct DirectJumpAggregation{T,S,F1,F2,RNG} <: AbstractSSAJumpAggregator
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
end
DirectJumpAggregation(nj::Int, njt::T, et::T, crs::Vector{T}, sr::T, maj::S, rs::F1, affs!::F2, sps::Tuple{Bool,Bool}, rng::RNG; kwargs...) where {T,S,F1,F2,RNG} =
  DirectJumpAggregation{T,S,F1,F2,RNG}(nj, nj, njt, et, crs, sr, maj, rs, affs!, sps, rng)

# creating the JumpAggregation structure (tuple-based constant jumps)
function aggregate(aggregator::Direct, u, p, t, end_time, constant_jumps,
    ma_jumps, save_positions, rng; kwargs...)

    # handle constant jumps using tuples
    rates, affects! = get_jump_info_tuples(constant_jumps)

    build_jump_aggregation(DirectJumpAggregation, u, p, t, end_time, ma_jumps,
            rates, affects!, save_positions, rng; kwargs...)
end

function initialize!(p::DirectJumpAggregation, integrator, u, params, t)
    generate_jumps!(p, integrator, u, params, t)
    nothing
end
# calculate the next jump / jump time
function generate_jumps!(p::DirectJumpAggregation, integrator, u, params, t)
p.sum_rate, ttnj = time_to_next_jump(p, u, params, t)
@fastmath p.next_jump_time = t + ttnj
@inbounds p.next_jump = searchsortedfirst(p.cur_rates, rand(p.rng) * p.sum_rate) # 这里用到了 searchsortedfirst 是一个加速的function
nothing
end


@fastmath function time_to_next_jump(p::DirectJumpAggregation{T,S,F1,F2,RNG}, u, params, t) where {T,S,F1 <: Tuple, F2 <: Tuple, RNG}
    prev_rate = zero(t)
    new_rate  = zero(t)
    cur_rates = p.cur_rates
  
    # mass action rates
    majumps   = p.ma_jumps
    idx       = get_num_majumps(majumps)
    @inbounds for i in 1:idx
      new_rate     = evalrxrate(u, i, majumps)
      cur_rates[i] = new_rate + prev_rate
      prev_rate    = cur_rates[i]
    end
  
    # constant jump rates
    rates = p.rates
    if !isempty(rates)
      idx  += 1
      fill_cur_rates(u, params, t, cur_rates, idx, rates...)
      @inbounds for i in idx:length(cur_rates)
        cur_rates[i] = cur_rates[i] + prev_rate
        prev_rate    = cur_rates[i]
      end
    end
  
    @inbounds sum_rate = cur_rates[end]
    sum_rate, randexp(p.rng) / sum_rate
end

@inline function fill_cur_rates(u, p, t, cur_rates, idx, rate, rates...)
    @inbounds cur_rates[idx] = rate(u, p, t)
    idx += 1
    fill_cur_rates(u, p, t, cur_rates, idx, rates...)
end
  
@inline function fill_cur_rates(u, p, t, cur_rates, idx, rate)
@inbounds cur_rates[idx] = rate(u, p, t)
nothing
end
  
function get_jump_info_tuples(constant_jumps)
    if (constant_jumps !== nothing) && !isempty(constant_jumps)
      rates    = ((c.rate for c in constant_jumps)...,)
      affects! = ((c.affect! for c in constant_jumps)...,)
    else
      rates    = ()
      affects! = ()
    end
  
    rates, affects!
end
