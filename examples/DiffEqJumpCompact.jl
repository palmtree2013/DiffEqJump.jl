module DiffEqJumpCompact
    using UnPack, StaticArrays
    using Random
    abstract type AbstractAggregatorAlgorithm end
    struct Direct <: AbstractAggregatorAlgorithm end

    abstract type AbstractJumpAggregator end
    abstract type AbstractJump end
    export SSAsolve, MassActionJump, ConstantRateJump
end