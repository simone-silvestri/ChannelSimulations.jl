using Oceananigans.Grids: architecture
using Oceananigans.Utils
using Oceananigans.Fields: Field
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y, _advective_tracer_flux_z
using Oceananigans.Advection: horizontal_advection_U, horizontal_advection_V
using Oceananigans.Models.HydrostaticFreeSurfaceModels: vertical_scaling, previous_vertical_scaling

import Oceananigans.Utils: KernelParameters

const c = Center()
const f = Face()

const ZStarSimulation = Simulation{<:HydrostaticFreeSurfaceModel{<:Any, <:Any, <:Any, <:Any, <:ZStarSpacingGrid}}

@inline function KernelParameters(f::Field)
    sz = size(f.data)
    of = f.data.offsets
    return KernelParameters(sz, of)
end

function update_fluxes!(simulation)
    model = simulation.model
    grid = model.grid
    arch = architecture(grid)

    b = model.tracers.b
    F = model.advection.b
    M = model.advection.momentum
    U = model.velocities
    A = model.auxiliary_fields

    params = KernelParameters(b)

    launch!(arch, grid, params, _update_fluxes!, A, U, b, grid, F, M)
    
    fill_halo_regions!(A)

    return nothing
end

@kernel function _update_fluxes!(auxiliary_fields, velocities, b, grid, advection, momentum_advection)
    i, j, k = @index(Global, NTuple)
    u, v, w = velocities

    fˣⁿ⁻¹ = auxiliary_fields.fˣⁿ⁻¹
    fʸⁿ⁻¹ = auxiliary_fields.fʸⁿ⁻¹
    fᶻⁿ⁻¹ = auxiliary_fields.fᶻⁿ⁻¹
    fˣⁿ⁻² = auxiliary_fields.fˣⁿ⁻²
    fʸⁿ⁻² = auxiliary_fields.fʸⁿ⁻²
    fᶻⁿ⁻² = auxiliary_fields.fᶻⁿ⁻²

    ζˣⁿ⁻¹ = auxiliary_fields.ζˣⁿ⁻¹
    ζʸⁿ⁻¹ = auxiliary_fields.ζʸⁿ⁻¹
    ζˣⁿ⁻² = auxiliary_fields.ζˣⁿ⁻²
    ζʸⁿ⁻² = auxiliary_fields.ζʸⁿ⁻²

    Uⁿ⁻¹ = auxiliary_fields.Uⁿ⁻¹
    Vⁿ⁻¹ = auxiliary_fields.Vⁿ⁻¹
    Wⁿ⁻¹ = auxiliary_fields.Wⁿ⁻¹
    Uⁿ⁻² = auxiliary_fields.Uⁿ⁻² 
    Vⁿ⁻² = auxiliary_fields.Vⁿ⁻²
    Wⁿ⁻² = auxiliary_fields.Wⁿ⁻²
    bⁿ⁻¹ = auxiliary_fields.bⁿ⁻¹
    ζⁿ⁻¹ = auxiliary_fields.ζⁿ⁻¹

    @inbounds begin
        # Save previous advective fluxes
        fˣⁿ⁻²[i, j, k] = fˣⁿ⁻¹[i, j, k] 
        fʸⁿ⁻²[i, j, k] = fʸⁿ⁻¹[i, j, k] 
        fᶻⁿ⁻²[i, j, k] = fᶻⁿ⁻¹[i, j, k] 
        
        ζˣⁿ⁻²[i, j, k] = ζˣⁿ⁻¹[i, j, k]
        ζʸⁿ⁻²[i, j, k] = ζʸⁿ⁻¹[i, j, k]
        
        # Save previous transport and previous buoyancy
        Uⁿ⁻²[i, j, k] = Uⁿ⁻¹[i, j, k]
        Vⁿ⁻²[i, j, k] = Vⁿ⁻¹[i, j, k]
        Wⁿ⁻²[i, j, k] = Wⁿ⁻¹[i, j, k]
        Uⁿ⁻¹[i, j, k] = u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid) 
        Vⁿ⁻¹[i, j, k] = v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid) 
        Wⁿ⁻¹[i, j, k] = w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid) 
        bⁿ⁻¹[i, j, k] = b[i, j, k]
        ζⁿ⁻¹[i, j, k] = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v)
    end

    @inbounds begin
        # Calculate new advective fluxes
        fˣⁿ⁻¹[i, j, k] = _advective_tracer_flux_x(i, j, k, grid, advection, u, b) * vertical_scaling(i, j, k, grid, c, c, c)
        fʸⁿ⁻¹[i, j, k] = _advective_tracer_flux_y(i, j, k, grid, advection, v, b) * vertical_scaling(i, j, k, grid, c, c, c)
        fᶻⁿ⁻¹[i, j, k] = _advective_tracer_flux_z(i, j, k, grid, advection, w, b) * vertical_scaling(i, j, k, grid, c, c, c) 
        
        ζʸⁿ⁻¹[i, j, k] = Ayᶠᶜᶜ(i, j, k, grid) * horizontal_advection_U(i, j, k, grid, momentum_advection, u, v)
        ζˣⁿ⁻¹[i, j, k] = Axᶜᶠᶜ(i, j, k, grid) * horizontal_advection_V(i, j, k, grid, momentum_advection, u, v)
    end
end

function assemble_P_values!(simulation)
    model = simulation.model
    grid = model.grid
    arch = architecture(grid)

    b = model.tracers.b
    u = model.velocities.u
    v = model.velocities.v
    
    bⁿ⁻¹ = model.auxiliary_fields.bⁿ⁻¹  
    ζⁿ⁻¹ = model.auxiliary_fields.ζⁿ⁻¹

    Uⁿ⁻¹ = model.auxiliary_fields.Uⁿ⁻¹
    Vⁿ⁻¹ = model.auxiliary_fields.Vⁿ⁻¹
    Wⁿ⁻¹ = model.auxiliary_fields.Wⁿ⁻¹
    Uⁿ⁻² = model.auxiliary_fields.Uⁿ⁻²   
    Vⁿ⁻² = model.auxiliary_fields.Vⁿ⁻²  
    Wⁿ⁻² = model.auxiliary_fields.Wⁿ⁻²

    Pu   = model.auxiliary_fields.Pu
    Pv   = model.auxiliary_fields.Pv
    Pw   = model.auxiliary_fields.Pw

    Pζu   = model.auxiliary_fields.Pζu
    Pζv   = model.auxiliary_fields.Pζv

    ∂xb² = model.auxiliary_fields.∂xb²
    ∂yb² = model.auxiliary_fields.∂yb²
    ∂zb² = model.auxiliary_fields.∂zb²

    fˣⁿ⁻¹ = simulation.model.auxiliary_fields.fˣⁿ⁻¹
    fʸⁿ⁻¹ = simulation.model.auxiliary_fields.fʸⁿ⁻¹
    fᶻⁿ⁻¹ = simulation.model.auxiliary_fields.fᶻⁿ⁻¹

    fˣⁿ⁻² = simulation.model.auxiliary_fields.fˣⁿ⁻²
    fʸⁿ⁻² = simulation.model.auxiliary_fields.fʸⁿ⁻²
    fᶻⁿ⁻² = simulation.model.auxiliary_fields.fᶻⁿ⁻²

    ζˣⁿ⁻¹ = simulation.model.auxiliary_fields.ζˣⁿ⁻¹
    ζʸⁿ⁻¹ = simulation.model.auxiliary_fields.ζʸⁿ⁻¹

    ζˣⁿ⁻² = simulation.model.auxiliary_fields.ζˣⁿ⁻²
    ζʸⁿ⁻² = simulation.model.auxiliary_fields.ζʸⁿ⁻²

    ∂xζ² = simulation.model.auxiliary_fields.∂xζ²
    ∂yζ² = simulation.model.auxiliary_fields.∂yζ²

    χ = simulation.model.timestepper.χ

    launch!(arch, grid, :xyz, _compute_dissipation!, 
            Pu, Pv, Pw, Pζu, Pζv,
            ∂xb², ∂yb², ∂zb², 
            ∂xζ², ∂yζ²,
            grid, χ, 
            Uⁿ⁻¹, Vⁿ⁻¹, Wⁿ⁻¹, 
            Uⁿ⁻², Vⁿ⁻², Wⁿ⁻², 
            fˣⁿ⁻¹, fʸⁿ⁻¹, fᶻⁿ⁻¹, 
            fˣⁿ⁻², fʸⁿ⁻², fᶻⁿ⁻²,
            ζˣⁿ⁻¹, ζʸⁿ⁻¹, ζˣⁿ⁻², ζʸⁿ⁻²,
            b, bⁿ⁻¹, u, v, ζⁿ⁻¹)

    return nothing
end

@kernel function _compute_dissipation!(Pu, Pv, Pw, Pζu, Pζv,
                                       ∂xb², ∂yb², ∂zb², 
                                       ∂xζ², ∂yζ²,
                                       grid, χ, 
                                       Uⁿ⁻¹, Vⁿ⁻¹, Wⁿ⁻¹, 
                                       Uⁿ⁻², Vⁿ⁻², Wⁿ⁻², 
                                       fˣⁿ⁻¹, fʸⁿ⁻¹, fᶻⁿ⁻¹, 
                                       fˣⁿ⁻², fʸⁿ⁻², fᶻⁿ⁻²,
                                       ζˣⁿ⁻¹, ζʸⁿ⁻¹, ζˣⁿ⁻², ζʸⁿ⁻²,
                                       b, bⁿ⁻¹, u, v, ζⁿ⁻¹)
    
    i, j, k = @index(Global, NTuple)

    @inbounds Pζu[i, j, k] = compute_Pζᵁ(i, j, k, grid, Uⁿ⁻¹, Uⁿ⁻², χ, ζˣⁿ⁻¹, ζˣⁿ⁻², u, v, ζⁿ⁻¹)
    @inbounds Pζv[i, j, k] = compute_Pζⱽ(i, j, k, grid, Vⁿ⁻¹, Vⁿ⁻², χ, ζʸⁿ⁻¹, ζʸⁿ⁻², u, v, ζⁿ⁻¹)

    @inbounds Pu[i, j, k] = compute_Pᵁ(i, j, k, grid, Uⁿ⁻¹, Uⁿ⁻², χ, fˣⁿ⁻¹, fˣⁿ⁻², b, bⁿ⁻¹)
    @inbounds Pv[i, j, k] = compute_Pⱽ(i, j, k, grid, Vⁿ⁻¹, Vⁿ⁻², χ, fʸⁿ⁻¹, fʸⁿ⁻², b, bⁿ⁻¹)
    @inbounds Pw[i, j, k] = compute_Pᵂ(i, j, k, grid, Wⁿ⁻¹, Wⁿ⁻², χ, fᶻⁿ⁻¹, fᶻⁿ⁻², b, bⁿ⁻¹)

    @inbounds ∂xb²[i, j, k] = Axᶠᶜᶜ(i, j, k, grid) * δxᶠᶜᶜ(i, j, k, grid, bⁿ⁻¹)^2 / Δxᶠᶜᶜ(i, j, k, grid)
    @inbounds ∂yb²[i, j, k] = Ayᶜᶠᶜ(i, j, k, grid) * δyᶜᶠᶜ(i, j, k, grid, bⁿ⁻¹)^2 / Δyᶜᶠᶜ(i, j, k, grid)
    @inbounds ∂zb²[i, j, k] = Azᶜᶜᶠ(i, j, k, grid) * δzᶜᶜᶠ(i, j, k, grid, bⁿ⁻¹)^2 / Δzᶜᶜᶠ(i, j, k, grid)

    @inbounds ∂xζ²[i, j, k] = Axᶜᶠᶜ(i, j, k, grid) * δxᶜᶠᶜ(i, j, k, grid, ζⁿ⁻¹)^2 / Δxᶜᶠᶜ(i, j, k, grid)
    @inbounds ∂yζ²[i, j, k] = Ayᶠᶜᶜ(i, j, k, grid) * δyᶠᶜᶜ(i, j, k, grid, ζⁿ⁻¹)^2 / Δyᶠᶜᶜ(i, j, k, grid)
end

@inline b★(i, j, k, grid, bⁿ, bⁿ⁻¹) = @inbounds (bⁿ[i, j, k] + bⁿ⁻¹[i, j, k]) / 2
@inline b²(i, j, k, grid, b₁, b₂)   = @inbounds (b₁[i, j, k] * b₂[i, j, k])

@inline ζ★(i, j, k, grid, uⁿ, vⁿ, ζⁿ⁻¹) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ, vⁿ) + ζⁿ⁻¹[i, j, k]) / 2
@inline ζ²(i, j, k, grid, uⁿ, vⁿ, ζⁿ⁻¹) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ, vⁿ) * ζⁿ⁻¹[i, j, k])

@inline function compute_Pζᵁ(i, j, k, grid, Uⁿ⁻¹, Uⁿ⁻², χ, fˣⁿ⁻¹, fˣⁿ⁻², u, v, ζⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δˣb★ = δxᶜᶠᶜ(i, j, k, grid, ζ★, u, v, ζⁿ⁻¹)
    δˣb² = δxᶜᶠᶜ(i, j, k, grid, ζ², u, v, ζⁿ⁻¹)

    @inbounds begin
        𝒰ⁿ⁻¹ = C₁ * ℑxyᶜᶠᵃ(i, j, k, grid, Uⁿ⁻¹)
        𝒰ⁿ⁻² = C₂ * ℑxyᶜᶠᵃ(i, j, k, grid, Uⁿ⁻²)
        Fⁿ⁻¹ = C₁ * fˣⁿ⁻¹[i, j, k] 
        Fⁿ⁻² = C₂ * fˣⁿ⁻²[i, j, k] 
        𝒜x = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟x = 𝒰ⁿ⁻¹ - 𝒰ⁿ⁻²
    end

    return 2 * δˣb★ * 𝒜x - δˣb² * 𝒟x
end

@inline function compute_Pζⱽ(i, j, k, grid, Vⁿ⁻¹, Vⁿ⁻², χ, fʸⁿ⁻¹, fʸⁿ⁻², u, v, ζⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δʸb★ = δyᶠᶜᶜ(i, j, k, grid, ζ★, u, v, ζⁿ⁻¹)
    δʸb² = δyᶠᶜᶜ(i, j, k, grid, ζ², u, v, ζⁿ⁻¹)

    @inbounds begin
        𝒱ⁿ⁻¹ = C₁ * ℑxyᶠᶜᵃ(i, j, k, grid, Vⁿ⁻¹)
        𝒱ⁿ⁻² = C₂ * ℑxyᶠᶜᵃ(i, j, k, grid, Vⁿ⁻²)
        Fⁿ⁻¹ = C₁ * fʸⁿ⁻¹[i, j, k] 
        Fⁿ⁻² = C₂ * fʸⁿ⁻²[i, j, k] 
        𝒜y = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟y = 𝒱ⁿ⁻¹ - 𝒱ⁿ⁻²  
    end

    return 2 * δʸb★ * 𝒜y - δʸb² * 𝒟y
end

@inline function compute_Pᵁ(i, j, k, grid, Uⁿ⁻¹, Uⁿ⁻², χ, fˣⁿ⁻¹, fˣⁿ⁻², bⁿ, bⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δˣb★ = δxᶠᶜᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δˣb² = δxᶠᶜᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        𝒰ⁿ⁻¹ = C₁ *  Uⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        𝒰ⁿ⁻² = C₂ *  Uⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻¹ = C₁ * fˣⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fˣⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜x = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟x = 𝒰ⁿ⁻¹ - 𝒰ⁿ⁻²
    end

    return 2 * δˣb★ * 𝒜x - δˣb² * 𝒟x
end

@inline function compute_Pⱽ(i, j, k, grid, Vⁿ⁻¹, Vⁿ⁻², χ, fʸⁿ⁻¹, fʸⁿ⁻², bⁿ, bⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δʸb★ = δyᶜᶠᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δʸb² = δyᶜᶠᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        𝒱ⁿ⁻¹ = C₁ *  Vⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        𝒱ⁿ⁻² = C₂ *  Vⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻¹ = C₁ * fʸⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fʸⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜y = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟y = 𝒱ⁿ⁻¹ - 𝒱ⁿ⁻²  
    end

    return 2 * δʸb★ * 𝒜y - δʸb² * 𝒟y
end

@inline function compute_Pᵂ(i, j, k, grid, Wⁿ⁻¹, Wⁿ⁻², χ, fᶻⁿ⁻¹, fᶻⁿ⁻², bⁿ, bⁿ⁻¹)
   
    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δᶻb★ = δzᶜᶜᶠ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δᶻb² = δzᶜᶜᶠ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        𝒲ⁿ⁻¹ = C₁ *  Wⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        𝒲ⁿ⁻² = C₂ *  Wⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻¹ = C₁ * fᶻⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fᶻⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜z = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟z = 𝒲ⁿ⁻¹ - 𝒲ⁿ⁻²  
    end

    return 2 * δᶻb★ * 𝒜z - δᶻb² * 𝒟z
end

