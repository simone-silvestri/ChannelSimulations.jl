using Oceananigans.Models: AbstractModel

const AB2Model = AbstractModel{<:QuasiAdamsBashforth2TimeStepper}
const RK3Model = AbstractModel{<:SplitRungeKutta3TimeStepper}

function assemble_dissipation!(model, dissipation)

    for tracer_name in keys(dissipation.advective_production)
        assemble_dissipation!(dissipation, model, tracer_name)
    end

    return nothing
end

@inline c★(i, j, k, grid, cⁿ⁺¹, cⁿ) = @inbounds (cⁿ⁺¹[i, j, k] + cⁿ[i, j, k]) / 2
@inline c²(i, j, k, grid, cⁿ⁺¹, cⁿ) = @inbounds (cⁿ⁺¹[i, j, k] * cⁿ[i, j, k])

@inline ζ★(i, j, k, grid, uⁿ⁺¹, vⁿ⁺¹, ζⁿ) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ⁺¹, vⁿ⁺¹) + ζⁿ[i, j, k]) / 2
@inline ζ²(i, j, k, grid, uⁿ⁺¹, vⁿ⁺¹, ζⁿ) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ⁺¹, vⁿ⁺¹) * ζⁿ[i, j, k])

function assemble_dissipation!(dissipation, model::AB2Model, tracer_name::Symbol)
    
    grid = model.grid
    arch = architecture(grid)
    χ = model.timestepper.χ

    # General velocities
    Uⁿ⁺¹ = model.velocities
    Uⁿ   = dissipation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation.previous_state.Uⁿ⁻¹
    cⁿ⁺¹ = tracer_name == :ζ ? nothing : model.tracers[tracer_name]
    cⁿ   = dissipation.previous_state[tracer_name]

    _assemble_advective_dissipation! = assemble_advective_ab2_dissipation_kernel(Val(tracer_name))
    _assemble_diffusive_dissipation! = assemble_diffusive_ab2_dissipation_kernel(Val(tracer_name))

    ####
    #### Assemble the advective dissipation
    ####

    P    = dissipation.advective_production[tracer_name]
    Fⁿ   = dissipation.advective_fluxes.Fⁿ[tracer_name]
    Fⁿ⁻¹ = dissipation.advective_fluxes.Fⁿ⁻¹[tracer_name]

    launch!(arch, grid, :xyz, _assemble_advective_dissipation!, P, grid, χ, Fⁿ, Fⁿ⁻¹, Uⁿ⁺¹, Uⁿ, Uⁿ⁻¹, cⁿ⁺¹, cⁿ)

    ####
    #### Assemble the diffusive dissipation
    #### 

    K    = dissipation.diffusive_production[tracer_name]
    Vⁿ   = dissipation.diffusive_fluxes.Vⁿ[tracer_name]
    Vⁿ⁻¹ = dissipation.diffusive_fluxes.Vⁿ⁻¹[tracer_name]

    launch!(arch, grid, :xyz, _assemble_diffusive_dissipation!, K, grid, χ, Vⁿ, Vⁿ⁻¹, Uⁿ⁺¹, cⁿ⁺¹, cⁿ)

    return nothing
end


function assemble_dissipation!(dissipation, model::RK3Model, tracer_name::Symbol)
    
    grid = model.grid
    arch = architecture(grid)

    # General velocities
    Uⁿ   = dissipation.previous_state.Uⁿ
    cⁿ⁺¹ = tracer_name == :ζ ? nothing : model.tracers[tracer_name]
    cⁿ   = dissipation.previous_state[tracer_name]

    _assemble_advective_dissipation! = assemble_advective_rk3_dissipation_kernel(Val(tracer_name))
    _assemble_diffusive_dissipation! = assemble_diffusive_rk3_dissipation_kernel(Val(tracer_name))

    γ = if model.clock.stage == 1 
        nothing
    elseif model.clock.stage == 2
        model.timestepper.γ²
    else
        model.timestepper.γ³
    end

    ####
    #### Assemble the advective dissipation
    ####

    P    = dissipation.advective_production[tracer_name]
    Fⁿ   = dissipation.advective_fluxes.Fⁿ[tracer_name]
    Fⁿ⁻¹ = dissipation.advective_fluxes.Fⁿ⁻¹[tracer_name]

    launch!(arch, grid, :xyz, _assemble_advective_dissipation!, P, grid, γ, Fⁿ, Uⁿ, cⁿ⁺¹, cⁿ)

    ####
    #### Assemble the diffusive dissipation
    #### 

    K    = dissipation.diffusive_production[tracer_name]
    Vⁿ   = dissipation.diffusive_fluxes.Vⁿ[tracer_name]
    Vⁿ⁻¹ = dissipation.diffusive_fluxes.Vⁿ⁻¹[tracer_name]

    launch!(arch, grid, :xyz, _assemble_diffusive_dissipation!, K, grid, γ, Vⁿ, cⁿ⁺¹, cⁿ)

    return nothing
end

assemble_advective_ab2_dissipation_kernel(val_tracer_name) = _assemble_advective_ab2_tracer_dissipation!
assemble_advective_ab2_dissipation_kernel(::Val{:ζ})       = _assemble_advective_ab2_vorticity_dissipation!
assemble_diffusive_ab2_dissipation_kernel(val_tracer_name) = _assemble_diffusive_ab2_tracer_dissipation!
assemble_diffusive_ab2_dissipation_kernel(::Val{:ζ})       = _assemble_diffusive_ab2_vorticity_dissipation!

assemble_advective_rk3_dissipation_kernel(val_tracer_name) = _assemble_advective_rk3_tracer_dissipation!
assemble_advective_rk3_dissipation_kernel(::Val{:ζ})       = _assemble_advective_rk3_vorticity_dissipation!
assemble_diffusive_rk3_dissipation_kernel(val_tracer_name) = _assemble_diffusive_rk3_tracer_dissipation!
assemble_diffusive_rk3_dissipation_kernel(::Val{:ζ})       = _assemble_diffusive_rk3_vorticity_dissipation!
