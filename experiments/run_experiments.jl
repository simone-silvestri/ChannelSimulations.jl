using ChannelSimulations
using ChannelSimulations: default_closure, default_momentum_advection, default_tracer_advection
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: WENOVectorInvariant

MOM=parse(Int, get(ENV, "MOM", "0"))
TRA=parse(Int, get(ENV, "TRA", "0"))
EXP=parse(Int, get(ENV, "EXP", "0"))

χ = 0.0

if MOM == 0
  momentum_advection = default_momentum_advection
  closure = default_closure
  restart_file = nothing 
elseif MOM == 1 
  momentum_advection = VectorInvariant()
  closure = (default_closure, HorizontalScalarBiharmonicDiffusivity(; ν = 9e8))
  restart_file = nothing 
elseif MOM == 2
  momentum_advection = WENOVectorInvariant(; vorticity_order = 5)
  closure = default_closure
  restart_file = "channel_checkpoint_20_iteration3542400.jld2"
elseif MOM == 3
  momentum_advection = default_momentum_advection
  closure = default_closure
  restart_file = "channel_checkpoint_30_iteration2894400.jld2"
  χ = 0.05
end

if TRA == 0
  tracer_advection = default_tracer_advection
elseif TRA == 1
  tracer_advection = TracerAdvection(WENO(; order = 5), WENO(; order = 5), Centered())
elseif TRA == 2
  tracer_advection = TracerAdvection(WENO(; order = 9), WENO(; order = 9), Centered())
elseif TRA == 3	
  tracer_advection = UpwindBiased(; order = 3)
elseif TRA == 4	
  tracer_advection = TracerAdvection(UpwindBiased(; order = 3), UpwindBiased(; order = 3), Centered())
else 
  tracer_advection = WENO()
end

simulation = run_channel_simulation(; closure, tracer_advection, momentum_advection, testcase = EXP, restart_file, χ) 
