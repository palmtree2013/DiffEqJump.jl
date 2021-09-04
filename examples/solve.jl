mutable struct SSAIntegrator{uType,tType,P,S,SA} 
    u::uType
    t::tType
    tprev::tType
    p::P
    sol::S
    i::Int
    saveat::SA
    cur_saveat::Int
    end_time::tType
    save_end::Bool
end

mutable struct SSASolution{tType,uType}
    t::Vector{tType}
    u::Vector{uType}
end

function initialize_integrator(u0::uType,t0::tType,p::P,tspan::Tuple,saveat::SA;save_end::Bool=true) where {uType,tType,P,SA}
    _u0 = copy(u0)
    _t0 = copy(t0)
    end_time = tspan[2]
    ssa_sol = SSASolution([t0],[u0])
    if typeof(saveat) <: Number
        _saveat = tspan[1]:saveat:tspan[2]
    else
        _saveat = saveat
    end
    if _saveat !== nothing && !isempty(_saveat) && _saveat[1] == tspan[1]
        cur_saveat = 2
    else
        cur_saveat = 1
    end
    integrator = SSAIntegrator(_u0,_t0,_t0,p,ssa_sol,0,_saveat,cur_saveat,end_time,save_end)
    integrator
end
function initialize_dja(integrator,aggregatoralgo::AbstractAggregatorAlgorithm,jumpset::JumpSet,rng)
    _u0 = copy(integrator.u)
    _t0 = copy(integrator.t)
    end_time = copy(integrator.end_time)
    dja = aggregate(aggregatoralgo,_u0,p,_t0,end_time,jumpset.constant_jumps, jumpset.massaction_jump, save_positions, rng)
    initialize!(dja, integrator, _u0, p, _t0)
    return dja
end
function SSAsolve(u0,t0,p,tspan,jumpset,aggregatoralgo::AbstractAggregatorAlgorithm;saveat=nothing,rng=Random.seed!(rand(UInt64)),save_end::Bool=true)
    end_time = tspan[2]
    integrator = initialize_integrator(u0,t0,p,tspan,saveat,save_end=save_end)
    dja = initialize_dja(integrator,aggregatoralgo,jumpset,rng)
    if integrator.saveat == nothing || isempty(integrator.saveat)
        while integrator.t < integrator.end_time
            dja(integrator)
            push!(integrator.sol.t,integrator.t)
            push!(integrator.sol.u,copy(integrator.u))
        end
    end    
    if integrator.saveat !== nothing && !isempty(integrator.saveat)
    # Split to help prediction
        while integrator.t < end_time
            while integrator.cur_saveat <= length(integrator.saveat) &&
                integrator.saveat[integrator.cur_saveat] < integrator.t
                push!(integrator.sol.t,integrator.saveat[integrator.cur_saveat])
                push!(integrator.sol.u,copy(integrator.u))
                integrator.cur_saveat += 1
            end
            dja(integrator)
        end
    end
    if integrator.save_end && integrator.sol.t[end] != end_time
        push!(integrator.sol.t,end_time)
        push!(integrator.sol.u,copy(integrator.u))
    end
    return integrator.sol
end