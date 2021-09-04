using DiffEqJump, Random
using OrdinaryDiffEq, Plots



rates1 = [0.1/1000.0]
reactant_stoich1 = [[1=>1,2=>1]]
net_stoich1 = [[1=>-1,2=>1]]
jump1_ = MassActionJump(rates1, reactant_stoich1, net_stoich1; scale_rates = false)

rate2 = (u,p,t) -> 0.01u[2]
affect! = function (integrator)
  integrator.u[2] -= 1
  integrator.u[3] += 1
end
jump2 = ConstantRateJump(rate2,affect!)

prob2 = DiscreteProblem([999.0,1.0,0.0],(0.0,250.0))
jump_prob2 = JumpProblem(prob2,Direct(),jump1_,jump2, save_positions=(false,false),saveat = 1, rng= Random.seed!(1234))
integrator = init(jump_prob2,FunctionMap())
# sol2 = solve(jump_prob2,SSAStepper(),saveat=1.)
# plot(sol2)
# jump_prob2.prob.p
# integrator.u|>typeof
integrator
prob2.f

integrator.saveat

prob2.u0

dja = jump_prob2.discrete_jump_aggregation
dja(integrator)
dja.next_jump
integrator


jump_prob2.discrete_jump_aggregation.next_jump

propertynames(prob2)
u0 = [999.0,1.0,0.0]
rng = Random.seed!(1234)
jumpset= JumpSet((),(jump2,),nothing,jump1_)

DiffEqJump.aggregate(Direct(),u0,prob2.p
,0.,250.,jumpset.constant_jumps,jumpset.massaction_jump,(false,false),rng)
# u0

jump_prob2|>propertynames
jump_prob2.discrete_jump_aggregation|>propertynames
jump_prob2.discrete_jump_aggregation.prev_jump
jump_prob2.discrete_jump_aggregation.rates
jump_prob2.discrete_jump_aggregation.affects!

jump_prob2.prob|>propertynames
jump_prob2.prob.f|propertynames
# integrator.p

integrator.f

u = [999.,1.,0.]

dja_ = jump_prob2.discrete_jump_aggregation
DiffEqJump.execute_jumps!(dja_, integrator, integrator.u, integrator.p, integrator.t)

dja_(integrator)
dja_(0,u,0,integrator)

propertynames(dja_)


integrator

integrator

DiffEqJump.generate_jumps!(dja_, integrator, integrator.u, integrator.p, integrator.t)

dja_

integrator

DiffEqJump.register_next_jump_time!(integrator, dja_, integrator.t)

integrator

DiffEqBase.add_tstop!(integrator,2)

integrator


DiffEqJump.initialize!(dja_,integrator,integrator.u,integrator.p,1)
dja_