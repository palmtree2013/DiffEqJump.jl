struct ConstantRateJump{F1,F2} <: AbstractJump
    rate::F1
    affect!::F2
end  


#TODO Simplify the code pmapper maybe
struct MassActionJump{T,S,U,V} <: AbstractJump
    scaled_rates::T
    reactant_stoch::S
    net_stoch::U
    param_mapper::V
  
    function MassActionJump{T,S,U,V}(rates::T, rs_in::S, ns::U, pmapper::V, scale_rates::Bool, useiszero::Bool, nocopy::Bool) where {T <: AbstractVector, S, U, V}
      sr  = nocopy ? rates : copy(rates)
      rs = nocopy ? rs_in : copy(rs_in)
      for i in eachindex(rs)
        if useiszero && (length(rs[i]) == 1) && iszero(rs[i][1][1])
          rs[i] = typeof(rs[i])()
        end
      end
  
      if scale_rates && !isempty(sr)
        scalerates!(sr, rs)
      end
      new(sr, rs, ns, pmapper)
    end
    function MassActionJump{Nothing,Vector{S},Vector{U},V}(::Nothing, rs_in::Vector{S}, ns::Vector{U}, pmapper::V, scale_rates::Bool, useiszero::Bool, nocopy::Bool) where {S<:AbstractVector, U<:AbstractVector, V}
      rs = nocopy ? rs_in : copy(rs_in)
      for i in eachindex(rs)
        if useiszero && (length(rs[i]) == 1) && iszero(rs[i][1][1])
          rs[i] = typeof(rs[i])()
        end
      end
      new(nothing, rs, ns, pmapper)
    end
    function MassActionJump{T,S,U,V}(rate::T, rs_in::S, ns::U, pmapper::V, scale_rates::Bool, useiszero::Bool, nocopy::Bool) where {T <: Number, S, U, V}
      rs = rs_in
      if useiszero && (length(rs) == 1) && iszero(rs[1][1])
        rs = typeof(rs)()
      end
      sr = scale_rates ? scalerate(rate, rs) : rate
      new(sr, rs, ns, pmapper)
    end
    function MassActionJump{Nothing,S,U,V}(::Nothing, rs_in::S, ns::U, pmapper::V, scale_rates::Bool, useiszero::Bool, nocopy::Bool) where {S, U, V}
      rs = rs_in
      if useiszero && (length(rs) == 1) && iszero(rs[1][1])
        rs = typeof(rs)()
      end
      new(nothing, rs, ns, pmapper)
    end
  
end
MassActionJump(usr::T, rs::S, ns::U, pmapper::V; scale_rates = true, useiszero = true, nocopy=false) where {T,S,U,V} = MassActionJump{T,S,U,V}(usr, rs, ns, pmapper, scale_rates, useiszero, nocopy)

MassActionJump(usr::T, rs, ns; scale_rates = true, useiszero = true, nocopy=false) where {T <: AbstractVector,S,U} = MassActionJump(usr, rs, ns, nothing; scale_rates=scale_rates, useiszero=useiszero, nocopy=nocopy)

MassActionJump(usr::T, rs, ns; scale_rates = true, useiszero = true, nocopy=false) where {T <: Number,S,U} = MassActionJump(usr, rs, ns, nothing; scale_rates=scale_rates, useiszero=useiszero, nocopy=nocopy)

# with parameter indices or mapping, multiple jump case
function MassActionJump(rs, ns; param_idxs=nothing, param_mapper=nothing, scale_rates = true, useiszero = true, nocopy=false)
  if param_mapper === nothing 
    (param_idxs === nothing) && error("If no parameter indices are given via param_idxs, an explicit parameter mapping must be passed in via param_mapper.")
    pmapper = MassActionJumpParamMapper(param_idxs)
  else
    (param_idxs !== nothing) && error("Only one of param_idxs and param_mapper should be passed.")
    pmapper = param_mapper
  end
                      
  MassActionJump(nothing, nocopy ? rs : copy(rs), ns, pmapper; scale_rates=scale_rates, 
                 useiszero=useiszero, nocopy=true)
end
@inline get_num_majumps(maj::MassActionJump) = length(maj.scaled_rates)
@inline get_num_majumps(maj::Nothing) = 0

function scalerates!(unscaled_rates::AbstractVector{U}, stochmat::AbstractVector{V}) where {U,S,T,W <: Pair{S,T}, V <: AbstractVector{W}}
    @inbounds for i in eachindex(unscaled_rates)
        coef = one(T)
        @inbounds for specstoch in stochmat[i]
            coef *= factorial(specstoch[2])
        end
        unscaled_rates[i] /= coef
    end
    nothing
end

function scalerate(unscaled_rate::U, stochmat::AbstractVector{Pair{S,T}}) where {U <: Number, S, T}
    coef = one(T)
    @inbounds for specstoch in stochmat
        coef *= factorial(specstoch[2])
    end
    unscaled_rate /= coef
end