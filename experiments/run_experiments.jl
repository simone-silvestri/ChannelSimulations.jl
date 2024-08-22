using ChannelSimulations
using ChannelSimulations: default_closure, default_momentum_advection, default_tracer_advection
using Oceananigans
using Oceananigans.Units

MOM=parse(Int, get(ENV, "MOM", "0"))
TRA=parse(Int, get(ENV, "TRA", "0"))
EXP=parse(Int, get(ENV, "EXP", "0"))

if MOM == 0
  momentum_advection = default_momentum_advection
  closure = default_closure
  restart_file = "channel_checkpoint_13_iteration0.jld2"
else
  momentum_advection = Centered()
  closure = (default_closure, HorizontalScalarBiharmonicDiffusivity(; ν = 9e8))
  restart_file = nothing
end

if TRA == 0
  tracer_advection = default_tracer_advection
elseif TRA == 1
  tracer_advection = TracerAdvection(WENO(; order = 5), WENO(; order = 5), Centered())
elseif TRA == 2
  tracer_advection = TracerAdvection(WENO(; order = 9), WENO(; order = 9), Centered())
else	
  tracer_advection = UpwindBiased(; order = 3)
end

simulation = run_channel_simulation(; closure, tracer_advection, momentum_advection, testcase = EXP, restart_file) 
