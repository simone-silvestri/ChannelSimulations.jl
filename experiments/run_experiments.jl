using ChannelSimulations
using ChannelSimulations: default_closure, 
                          default_momentum_advection, 
                          default_tracer_advection,
                          default_catke
using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: WENOVectorInvariant

CLO=parse(Int, get(ENV, "CLO", "0"))
MOM=parse(Int, get(ENV, "MOM", "0"))
TRA=parse(Int, get(ENV, "TRA", "0"))

EXP = string(clo) * string(mom) * string(tra)

χ = 0.05

if CLO == 0
  closure = default_closure
elseif CLO == 1
  closure = default_catke()
end  

if MOM == 0
  momentum_advection = default_momentum_advection
  restart_file = nothing 
elseif MOM == 1 
  momentum_advection = VectorInvariant()
  closure = (closure, HorizontalScalarBiharmonicDiffusivity(; ν = 9e8))
  restart_file = nothing 
elseif MOM == 2
  momentum_advection = WENOVectorInvariant(; vorticity_order = 5)
  restart_file = nothing 
elseif MOM == 3
  momentum_advection = default_momentum_advection
  restart_file = nothing 
end

if TRA == 0
  tracer_advection = default_tracer_advection
elseif TRA == 1
  tracer_advection = FluxFormAdvection(WENO(; order = 5), WENO(; order = 5), Centered())
elseif TRA == 2
  tracer_advection = FluxFormAdvection(WENO(; order = 9), WENO(; order = 9), Centered())
elseif TRA == 3	
  tracer_advection = UpwindBiased(; order = 3)
elseif TRA == 4	
  tracer_advection = FluxFormAdvection(UpwindBiased(; order = 3), UpwindBiased(; order = 3), Centered())
else 
  tracer_advection = WENO()
end

simulation = run_channel_simulation(; closure, tracer_advection, momentum_advection, testcase = EXP, restart_file, χ) 
