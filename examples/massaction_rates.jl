@inline @fastmath function evalrxrate(speciesvec::AbstractVector{T}, rxidx::S, majump::MassActionJump{U,V,W,X})::R where {T,S,R,U <: AbstractVector{R},V,W,X}
    val = one(T)
    @inbounds for specstoch in majump.reactant_stoch[rxidx]
        specpop = speciesvec[specstoch[1]]
        val    *= specpop
        @inbounds for k = 2:specstoch[2]
            specpop -= one(specpop)
            val     *= specpop
        end
    end

    @inbounds return val * majump.scaled_rates[rxidx]
end

@inline @fastmath function executerx!(speciesvec::AbstractVector{T}, rxidx::S,
    majump::MassActionJump{U,V,W,X}) where {T,S,U,V,W,X}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
    speciesvec[specstoch[1]] += specstoch[2]
    end
    nothing
end