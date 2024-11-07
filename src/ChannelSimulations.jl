module ChannelSimulations

# Write your package code here.

export run_channel_simulation, run_spindown_simulation

using Printf
using Statistics

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures: FivePointHorizontalFilter

using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: VelocityStencil
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Grids: xnode, ynode, znode, ZStarVerticalCoordinate
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ZStarSpacingGrid, VelocityFields
using Oceananigans.Utils: ConsecutiveIterations
using KernelAbstractions: @kernel, @index

default_closure = ConvectiveAdjustmentVerticalDiffusivity(background_κz = 0,
                                                          convective_κz = 0.1,
                                                          background_νz = 3e-4,
                                                          convective_νz = 0.1)

default_momentum_advection = VectorInvariant(vertical_scheme   = WENO(),
                                             vorticity_scheme  = WENO(; order = 9),
                                             divergence_scheme = WENO())

default_tracer_advection = WENO(; order = 7)

include("VarianceDissipationComputation/VarianceDissipationComputation.jl")

using .VarianceDissipationComputation

include("channel_simulation.jl")
include("spindown_simulation.jl")

end
