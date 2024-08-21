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
using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.Advection: TracerAdvection
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ZStar, ZStarSpacingGrid, VelocityFields
using Oceananigans.Utils: ConsecutiveIterations
using KernelAbstractions: @kernel, @index

include("compute_dissipation.jl")
include("channel_simulation.jl")
include("spindown_simulation.jl")
include("visualize.jl")

end
