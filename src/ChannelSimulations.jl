module ChannelSimulations

# Write your package code here.

export run_channel_simulation, run_spindown_simulation

using Printf
using Statistics

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEEquation, CATKEVerticalDiffusivity

using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: VelocityStencil
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans.Grids
using Oceananigans.Models
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VelocityFields
using Oceananigans.Utils: ConsecutiveIterations
using KernelAbstractions: @kernel, @index

default_closure = ConvectiveAdjustmentVerticalDiffusivity(background_κz = 1e-5,
                                                          convective_κz = 0.0,
                                                          background_νz = 3e-4,
                                                          convective_νz = 0.1)

function default_catke()
    mixing_length = CATKEMixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
end

default_momentum_advection = WENOVectorInvariant()
default_tracer_advection = WENO(order = 7)

include("VarianceDissipations/VarianceDissipations.jl")

using .VarianceDissipations

include("channel_simulation.jl")
include("spindown_simulation.jl")
include("one_dimensional_simulation.jl")

end
