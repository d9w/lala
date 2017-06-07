using JSON
include("cmaes.jl")

defaults = JSON.parsefile("defaults.json")
shapedefaults = JSON.parsefile("shape_defaults.json")
ranges = JSON.parsefile("ranges.json");

seed = length(ARGS) > 0 ? parse(Int64, ARGS[1]) : 0
srand(seed)

deltas = 0.5 * ones(length(ranges))
cparams = rand(length(ranges))

function printcfg(a::Array{Float64}, s::Int64)
  for i in eachindex(a)
    if (a[i] < 0.0) || (a[i] > 1.0)
      a[i] = mod(a[i], 1.0)
    end
  end
  i = 1
  for (k, v) in ranges
    defaults[k] = v[1]+(v[2]-v[1])*a[i]
    i += 1
  end
  for (k, v) in shapedefaults
    defaults[k] = v
  end
  defaults["seed"] = s
  open("config.json", "w") do io
    JSON.print(io, defaults)
  end
end

function readlog(logname::String)
  res = readlines(pipeline(`cat $logname`, `grep "Reward"`))
  numres = Array{Float64}(length(res),length(split(chomp(split(res[1],"::")[2]),",")))
  for i in eachindex(res)
    numres[i,:] = map(x->parse(Float64,x), split(chomp(split(res[i],"::")[2]),","))
  end
# numres, -(numres[end,4]*1000+sum(numres[:,1]))
  numres, -sum(numres[:,1])
end

function fitness(a::Array{Float64})
  fitness = 0
  for i=1:1
    printcfg(a, seed+i)
    run(pipeline(`./console -f config.json`, stdout="out.log", stderr=string("err",i,".log")))
    fitness += readlog(string("err",i,".log"))[2]
  end
  fitness
end

pmin=cmaes(fitness, cparams, deltas, stopeval=1000)
