"""
An aggregator interface for SSA-like algorithms.

### Required Fields
- `next_jump`          # the next jump to execute
- `prev_jump`          # the previous jump that was executed
- `next_jump_time`     # the time of the next jump
- `end_time`           # the time to stop a simulation
- `cur_rates`          # vector of current propensity values
- `sum_rate`           # sum of current propensity values
- `ma_jumps`           # any MassActionJumps for the system (scalar form)
- `rates`              # vector of rate functions for ConstantRateJumps
- `affects!`           # vector of affect functions for ConstantRateJumps
- `save_positions`     # tuple for whether to save the jumps before and/or after event
- `rng`                # random number generator

### Optional fields:
- `dep_gr`             # dependency graph, dep_gr[i] = indices of reactions that should
                       # be updated when rx i occurs.    
"""
abstract type AbstractSSAJumpAggregator <: AbstractJumpAggregator end

function (p::AbstractSSAJumpAggregator)(integrator)
    execute_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    register_next_jump_time!(p,integrator)
    generate_jumps!(p, integrator, integrator.u, integrator.p, integrator.t)
    nothing
end

"""
self_defined
"""
@inline function register_next_jump_time!(p::DirectJumpAggregation,integrator)
    integrator.tprev = integrator.t
    integrator.t = p.next_jump_time
    integrator.i +=1
    nothing
end

function build_jump_aggregation(jump_agg_type, u, p, t, end_time, ma_jumps, rates,
    affects!, save_positions, rng; kwargs...)

    # mass action jumps
    majumps = ma_jumps
    if majumps === nothing
        majumps = MassActionJump(Vector{typeof(t)}(),
        Vector{Vector{Pair{Int,eltype(u)}}}(),
        Vector{Vector{Pair{Int,eltype(u)}}}())
    end

    # current jump rates, allows mass action rates and constant jumps
    cur_rates = Vector{typeof(t)}(undef, get_num_majumps(majumps) + length(rates))

    sum_rate = zero(typeof(t))
    next_jump = 0
    next_jump_time = typemax(typeof(t))
    jump_agg_type(next_jump, next_jump_time, end_time, cur_rates, sum_rate,
    majumps, rates, affects!, save_positions, rng; kwargs...)
end

"""
fill_rates_and_sum!(p::AbstractSSAJumpAggregator, u, params, t)

Reevaluate all rates and their sum.
"""
function fill_rates_and_sum!(p::AbstractSSAJumpAggregator, u, params, t)
    sum_rate = zero(typeof(p.sum_rate))

    # mass action jumps
    majumps   = p.ma_jumps
    cur_rates = p.cur_rates
    @inbounds for i in 1:get_num_majumps(majumps)
        cur_rates[i] = evalrxrate(u, i, majumps)
        sum_rate    += cur_rates[i]
    end

    # constant rates
    rates = p.rates
    idx   = get_num_majumps(majumps) + 1
    @inbounds for rate in rates
        cur_rates[idx] = rate(u, params, t)
        sum_rate += cur_rates[idx]
        idx += 1
    end

    p.sum_rate = sum_rate
    nothing
end


@inline function update_state!(p::AbstractSSAJumpAggregator, integrator, u)
    @unpack ma_jumps, next_jump = p
    num_ma_rates = get_num_majumps(ma_jumps)
    if next_jump <= num_ma_rates # is next jump a mass action jump
        if u isa SVector
            integrator.u = executerx(u, next_jump, ma_jumps)
        else
            @inbounds executerx!(u, next_jump, ma_jumps)
        end
    else
        idx = next_jump - num_ma_rates
        @inbounds p.affects![idx](integrator)
    end

    # save jump that was just executed
    p.prev_jump = next_jump
    return integrator.u
end
