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
rates = [.1,0.25,10.,1.2]
prob = DiscreteProblem(u0,(0.0,200.0))
majumps = MassActionJump(rates,reactstoch,netstoch)

jump_prob_MNRM = JumpProblem(prob,MNRM(),majumps,save_positions=(true,false))
jump_prob_NRM = JumpProblem(prob,NRM(),majumps,save_positions=(false,false))
jump_prob_Direct = JumpProblem(prob,Direct(),majumps,save_positions=(true,false))

integrator = DiffEqBase.__init(jump_prob_Direct,SSAStepper(),seed=1234)
aggregator = integrator.cb.initialize
# aggregator(integrator)
aggregator
# integrator.sol.t
integrator.tstop

sol =@time solve(jump_prob_Direct,SSAStepper())

sol_NRM =@time solve(jump_prob_NRM,SSAStepper(),saveat=1.)
sol.u
sol_NRM.u

ens_sol_MNRM =@time solve(EnsembleProblem(jump_prob_MNRM), SSAStepper(), EnsembleThreads(); trajectories=5e4)
ens_sol_NRM =@time solve(EnsembleProblem(jump_prob_NRM), SSAStepper(), EnsembleThreads(); trajectories=5e4)
ens_sol_Direct =@time solve(EnsembleProblem(jump_prob_Direct), SSAStepper(), EnsembleThreads(); trajectories=5e4)

using DifferentialEquations.EnsembleAnalysis
end_P = componentwise_vectors_timepoint(ens_sol_MNRM,200.)[3]
end_P_D = componentwise_vectors_timepoint(ens_sol_Direct,200.)[3]

using Plots
histogram(end_P, bins = 0:1:20, normalize = :pdf)
histogram!(end_P_D, bins = 0:1:20, normalize = :pdf)
