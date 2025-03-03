using Oceananigans: fields

function update_fluxes!(model, dissipation)
    grid   = model.grid
    arch   = architecture(grid)
    params = KernelParameters(model.tracers[1])
    
    Uⁿ   = dissipation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation.previous_state.Uⁿ⁻¹ 
    U    = model.velocities

    launch!(arch, grid, params, _update_transport!, Uⁿ, Uⁿ⁻¹, grid, U)

    for (tracer_id, tracer_name) in enumerate(keys(dissipation.advective_production))
        update_fluxes!(dissipation, model, tracer_name, tracer_id)
    end

    return nothing
end

@inline parameters(grid, i) = topology(grid, i) == Flat ? UnitRange(1, 1) : UnitRange(-1, size(grid, i)+1)

function update_fluxes!(dissipation, model, tracer_name::Symbol, tracer_id)
    
    # Grab tracer properties
    cⁿ⁺¹ = model.tracers[tracer_name]
    cⁿ   = dissipation.previous_state[tracer_name]

    grid = model.grid
    arch = architecture(grid)
    Nx, Ny, Nz = size(grid)
    U = model.velocities

    # Include boundary regions
    Px = parameters(grid, 1)
    Py = parameters(grid, 2)
    Pz = parameters(grid, 3)

    params = KernelParameters(Px, Py, Pz)

    _update_advective_fluxes! = update_advective_fluxes_kernel(Val(tracer_name))
    _update_diffusive_fluxes! = update_diffusive_fluxes_kernel(Val(tracer_name))

    ####
    #### Update the advective fluxes and compute gradient squared
    ####

    Fⁿ   = dissipation.advective_fluxes.Fⁿ[tracer_name]
    Fⁿ⁻¹ = dissipation.advective_fluxes.Fⁿ⁻¹[tracer_name]
    Gⁿ   = dissipation.gradient_squared[tracer_name]
    advection = getadvection(model.advection, tracer_name)

    launch!(arch, grid, params, _update_advective_fluxes!, Gⁿ, Fⁿ, Fⁿ⁻¹, cⁿ, grid, advection, U, cⁿ⁺¹)

    ####
    #### Update the diffusive fluxes
    ####

    Vⁿ   = dissipation.diffusive_fluxes.Vⁿ[tracer_name]
    Vⁿ⁻¹ = dissipation.diffusive_fluxes.Vⁿ⁻¹[tracer_name]

    D = model.diffusivity_fields
    B = model.buoyancy
    clk  = model.clock
    clo  = model.closure
    model_fields = fields(model)

    launch!(arch, grid, params, _update_diffusive_fluxes!, Vⁿ, Vⁿ⁻¹, grid, clo, D, B, cⁿ⁺¹, Val(tracer_id), clk, model_fields)

    return nothing
end

@kernel function _update_transport!(Uⁿ, Uⁿ⁻¹, grid, U)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Uⁿ⁻¹.u[i, j, k] = Uⁿ.u[i, j, k]
        Uⁿ⁻¹.v[i, j, k] = Uⁿ.v[i, j, k]
        Uⁿ⁻¹.w[i, j, k] = Uⁿ.w[i, j, k]
          Uⁿ.u[i, j, k] = U.u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid) 
          Uⁿ.v[i, j, k] = U.v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) 
          Uⁿ.w[i, j, k] = U.w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid) 
    end
end

update_advective_fluxes_kernel(val_tracer_name) = _update_advective_tracer_fluxes!
update_advective_fluxes_kernel(::Val{:ζ})       = _update_advective_vorticity_fluxes!
update_diffusive_fluxes_kernel(val_tracer_name) = _update_diffusive_tracer_fluxes!
update_diffusive_fluxes_kernel(::Val{:ζ})       = _update_diffusive_vorticity_fluxes!

