using ModelingToolkit, Catalyst, Symbolics, Statistics, DataFrames
function create_neighbour(N;dims=1)
    dims_=Tuple(fill(N,dims))
    dict = Dict{CartesianIndex{dims},Vector{CartesianIndex{dims}}}()
    cart = CartesianIndices(dims_)
    for ind in cart
        a = Tuple(ind)
        neighbour = []
        for shift in [1,-1]
            for i in 1:dims
                b = Int.(zeros(dims))
                b[i] = shift
                c = CartesianIndex(a .+ Tuple(b))
                if c in cart
                    push!(neighbour,c)
                end
            end
        end        
        push!(dict,ind=>neighbour)
    end
    return dict
end

N = 128
dict_neighbour = @time create_neighbour(N,dims=1)
@parameters  D t
X = [Num(Symbolics.variable(:X,i; T = ModelingToolkit.FnType{Tuple{Any},Real}))(t) for i in 1:N]
begin # construct reaction network 
    rxs_tran = []
    for ind in CartesianIndices(X)
        for ind_neig in dict_neighbour[ind]
            push!(rxs_tran,Reaction(D,[X[ind]],[X[ind_neig]]))
        end
    end
    rxs = vcat(rxs_tran)
    @named rs = ReactionSystem(rxs,t,vec(X),[D])         
end

h = 1 / N
u0 = 10*ones(Int64, N)
tf = .01
methods = (Direct(),MNRM(),NRM(),DirectCR(),RSSA())
shortlabels = [string(leg)[12:end-2] for leg in methods]
prob    = DiscreteProblem(rs, u0, (0.0, tf), [1 / (h*h)])
ploth   = plot(reuse=false)
for (i,method) in enumerate(methods)
    println("Benchmarking method: ", method)
    jump_prob = JumpProblem(rs, prob, method, save_positions=(false, false))
    sol = solve(jump_prob, SSAStepper(), saveat=tf/1000.)
    plot!(ploth,sol.t,sol[Int(N//2),:],label=shortlabels[i], format=fmt)
end
plot!(ploth, title="Population at middle lattice site", xlabel="time",format=fmt)

function run_benchmark!(t, jump_prob, stepper)
    sol = solve(jump_prob, stepper)
    @inbounds for i in 1:length(t)
        t[i] = @elapsed (sol = solve(jump_prob, stepper))
    end
end


nsims = 50
benchmarks = Vector{Vector{Float64}}()
for method in methods
    jump_prob = JumpProblem(rs, prob, method, save_positions=(false, false))
    stepper = SSAStepper()
    t = Vector{Float64}(undef,nsims)
    run_benchmark!(t, jump_prob, stepper)
    push!(benchmarks, t)
end

medtimes = Vector{Float64}(undef,length(methods))
stdtimes = Vector{Float64}(undef,length(methods))
avgtimes = Vector{Float64}(undef,length(methods))
for i in 1:length(methods)
    medtimes[i] = median(benchmarks[i])
    avgtimes[i] = mean(benchmarks[i])
    stdtimes[i] = std(benchmarks[i])
end

df = DataFrame(names=shortlabels,medtimes=medtimes,relmedtimes=(medtimes/medtimes[1]),
                avgtimes=avgtimes, std=stdtimes, cv=stdtimes./avgtimes)


sa = [string(round(mt,digits=4),"s") for mt in df.medtimes]
bar(df.names,df.relmedtimes,legend=:false, fmt=fmt)
scatter!(df.names, .05 .+ df.relmedtimes, markeralpha=0, series_annotations=sa, fmt=fmt)
ylabel!("median relative to Direct")
title!("256 Site 1D Diffusion CTRW")