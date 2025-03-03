module VarianceDissipations
 
export VarianceDissipation, get_dissipation_fields

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Grids: architecture, topology
using Oceananigans.Utils
using Oceananigans.TimeSteppers
using Oceananigans.Fields
using Oceananigans.Fields: Field, VelocityFields
using Oceananigans.Operators
using Oceananigans.Operators: σⁿ, σ⁻
using Oceananigans.BoundaryConditions
using Oceananigans.TurbulenceClosures: viscosity,
                                       diffusivity, 
                                       ScalarDiffusivity, 
                                       ExplicitTimeDiscretization,
                                       ScalarBiharmonicDiffusivity,
                                       AbstractTurbulenceClosure,
                                       HorizontalFormulation,
                                       _diffusive_flux_x,
                                       _diffusive_flux_y,
                                       _diffusive_flux_z,
                                       diffusive_flux_x, 
                                       diffusive_flux_y, 
                                       diffusive_flux_z,
                                       ∂ⱼ_τ₁ⱼ, ∂ⱼ_τ₂ⱼ

using Oceananigans.Advection: _advective_tracer_flux_x, 
                              _advective_tracer_flux_y, 
                              _advective_tracer_flux_z

using Oceananigans.TimeSteppers: SplitRungeKutta3TimeStepper

using Oceananigans.Operators: volume
using KernelAbstractions: @kernel, @index

struct VarianceDissipation{P, K, A, D, S, G}
    advective_production :: P
    diffusive_production :: K
    advective_fluxes :: A
    diffusive_fluxes :: D
    previous_state :: S
    gradient_squared :: G
end

include("dissipation_utils.jl")

function VarianceDissipation(model; tracers = propertynames(model.tracers))
    tracers = tupleit(tracers)
    diffusivities = model.diffusivity_fields
    closure       = model.closure
    grid          = model.grid

    P    = NamedTuple{tracers}(tracer_fluxes(grid) for tracer in tracers)
    Fⁿ   = NamedTuple{tracers}(tracer_fluxes(grid) for tracer in tracers)
    Fⁿ⁻¹ = NamedTuple{tracers}(tracer_fluxes(grid) for tracer in tracers)
    
    K    = NamedTuple{tracers}(tracer_closure_dissipation(grid, diffusivities, closure, id) for id in eachindex(tracers))
    Vⁿ   = NamedTuple{tracers}(tracer_closure_dissipation(grid, diffusivities, closure, id) for id in eachindex(tracers))
    Vⁿ⁻¹ = NamedTuple{tracers}(tracer_closure_dissipation(grid, diffusivities, closure, id) for id in eachindex(tracers))    
    
    Uⁿ⁻¹ = VelocityFields(grid)
    Uⁿ   = VelocityFields(grid)
    
    cⁿ⁻¹ =  NamedTuple{tracers}(CenterField(grid) for tracer in tracers)

    previous_state   = merge(cⁿ⁻¹, (; Uⁿ⁻¹, Uⁿ))
    advective_fluxes = (; Fⁿ, Fⁿ⁻¹)
    diffusive_fluxes = (; Vⁿ, Vⁿ⁻¹)

    gradients = deepcopy(P)

    return VarianceDissipation(P, K, advective_fluxes, diffusive_fluxes, previous_state, gradients)
end

# Function to call in a callback
# Note: This works only if the callback is called with an IterationInterval(1), if not the
# previous fluxes and velocities will not be correct
# TODO: make sure that the correct velocities and fluxes are used even if 
# the callback is not called with an IterationInterval(1)
function (ϵ::VarianceDissipation)(model)

    # We first assemble values for Pⁿ⁻¹
    assemble_dissipation!(model, ϵ)

    # Then we update the fluxes to be used in the next time step
    update_fluxes!(model, ϵ)

    return nothing
end

const f = Face()
const c = Center()

include("get_dissipation_fields.jl")
include("advective_fluxes.jl")
include("diffusive_fluxes.jl")
include("update_fluxes.jl")
include("advective_dissipation.jl")
include("diffusive_dissipation.jl")
include("assemble_dissipation.jl")

end
