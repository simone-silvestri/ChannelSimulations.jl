using Oceananigans
using Oceananigans.Utils
using Oceananigans.TimeSteppers: implicit_step!
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

import Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: time_step_catke_equation! 

function time_step_catke_equation!(model)

    if model.closure isa Tuple
        catke_closure_idx = findfirst(clo -> clo isa CATKEVerticalDiffusivity, model.closure)
        
        # If there is no CATKE... do nothing!
        if isnothing(catke_closure_idx)
            return nothing
        end
        closure = model.closure[catke_closure_idx]

    # If there is no CATKE... do nothing!
    elseif !(model.closure isa CATKEVerticalDiffusivity)
        return nothing
    
    else
        closure = model.closure

    end

    # TODO: properly handle closure tuples
    e = model.tracers.e
    arch = model.architecture
    grid = model.grid
    Gⁿe = model.timestepper.Gⁿ.e
    G⁻e = model.timestepper.G⁻.e

    diffusivity_fields = model.diffusivity_fields
    κe = diffusivity_fields.κe
    Le = diffusivity_fields.Le
    previous_velocities = diffusivity_fields.previous_velocities
    tracer_index = findfirst(k -> k == :e, keys(model.tracers))
    implicit_solver = model.timestepper.implicit_solver

    Δt = model.clock.last_Δt
    Δτ = get_time_step(closure)

    if isnothing(Δτ)
        Δτ = Δt
        M = 1
    else
        M = ceil(Int, Δt / Δτ) # number of substeps
        Δτ = Δt / M
    end

    FT = eltype(grid)

    for m = 1:M # substep
        if m == 1 && M != 1
            χ = convert(FT, -0.5) # Euler step for the first substep
        else
            χ = model.timestepper.χ
        end

        # Compute the linear implicit component of the RHS (diffusivities, L)
        # and step forward
        launch!(arch, grid, :xyz,
                substep_turbulent_kinetic_energy!,
                κe, Le, grid, closure,
                model.velocities, previous_velocities, # try this soon: model.velocities, model.velocities,
                model.tracers, model.buoyancy, diffusivity_fields,
                Δτ, χ, Gⁿe, G⁻e)

        # Good idea?
        # previous_time = model.clock.time - Δt
        # previous_iteration = model.clock.iteration - 1
        # current_time = previous_time + m * Δτ
        # previous_clock = (; time=current_time, iteration=previous_iteration)

        implicit_step!(e, implicit_solver, closure,
                       model.diffusivity_fields, Val(tracer_index),
                       model.clock, Δτ)
    end

    return nothing
end

