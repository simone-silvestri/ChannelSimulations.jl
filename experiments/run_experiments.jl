using ChannelSimulations
using Oceananigans
using Oceananigans.Units

default_closure = ConvectiveAdjustmentVerticalDiffusivity(background_κz = 0,
                                                          convective_κz = 0.1,
                                                          background_νz = 1e-4,
                                                          convective_νz = 0.1)
ν = 9e8
momentum_advection = Centered()
closure = (default_closure, HorizontalScalarBiharmonicDiffusivity(; ν))

# Upwind 3
tracer_advection = UpwindBiased(; order = 3)
simulation = run_channel_simulation!(; closure, tracer_advection, momentum_advection, testcase = "mitgcm_upwind3_second")
# WENO 9 
tracer_advection = WENO(; order = 9)
simulation = run_channel_simulation!(; closure, tracer_advection, momentum_advection, testcase = "mitgcm_weno9_second")
