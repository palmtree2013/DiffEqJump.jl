using DiffEqJump, DiffEqBase, OrdinaryDiffEq
using Test
# using StableRNGs
# rng = StableRNG(12345)

# k10, G -> G*
# k01, G* -> G
# ρ,  G -> G + P
# d, P -> 0

u0 = [0,1,0]
reactstoch = [[1=>1],[2=>1],[1=>1],[3=>1]]
netstoch = [[1=>-1,2=>1],[1=>1,2=>-1],[3=>1],[3=>-1]]
prob = DiscreteProblem(u0,(0.0,250.0))
rates = [.1,0.25,20.,1.]
majumps = MassActionJump(rates,reactstoch,netstoch)
jump_prob_MNRM = JumpProblem(prob,MNRM(),majumps,save_positions=(false,false))
jump_prob_Direct = JumpProblem(prob,DirectCR(),majumps,save_positions=(false,false))
# JumpSet(majumps)
# DiffEqJump.using_params(majumps)

# sol = solve(jump_prob_MNRM,SSAStepper(),saveat=1.)
# sol.u

ens_sol_MNRM =@time solve(EnsembleProblem(jump_prob_MNRM), SSAStepper(), EnsembleThreads(); trajectories=1e4)
ens_sol_Direct =@time solve(EnsembleProblem(jump_prob_Direct), SSAStepper(), EnsembleThreads(); trajectories=1e4)

using DifferentialEquations.EnsembleAnalysis
end_P = componentwise_vectors_timepoint(ens_sol_MNRM,250.)[3]
end_P_D = componentwise_vectors_timepoint(ens_sol_Direct,250.)[3]

using Plots
histogram(end_P)
histogram!(end_P_D)
