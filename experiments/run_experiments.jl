using ChannelSimulations
using ChannelSimulations: default_closure
using Oceananigans
using Oceananigans.Units

run_channel_simulation!()

tracer_advection = TracerAdvection(WENO(; order = 9), WENO(; order = 9), Centered()) 

run_channel_simulation!(; tracer_advection, testcase = "1")

ν = 5000^4 / 15days
closure = (default_closure, HorizontalScalarBiharmonicDiffusivity(; ν))

run_channel_simulation!(; closure, testcase = "2")

momentum_advection = VectorInvariant()

run_channel_simulation!(; closure, momentum_advection, testcase = "3")

tracer_advection = Centered()

closure = (default_closure, HorizontalScalarBiharmonicDiffusivity(; ν, κ = ν))

run_channel_simulation!(; closure, momentum_advection, tracer_advection, testcase = "4")
