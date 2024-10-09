using Oceananigans.Grids: architecture
using Oceananigans.Utils
using Oceananigans.Fields: Field
using Oceananigans.Operators
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y, _advective_tracer_flux_z
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
    U = model.velocities
    A = model.auxiliary_fields

    params = KernelParameters(b)

    launch!(arch, grid, params, _update_fluxes!, A, U, b, grid, F)
    
    return nothing
end

@kernel function _update_fluxes!(auxiliary_fields, velocities, b, grid, advection)
    i, j, k = @index(Global, NTuple)
    u, v, w = velocities

    fˣⁿ⁻¹ = auxiliary_fields.fˣⁿ⁻¹
    fʸⁿ⁻¹ = auxiliary_fields.fʸⁿ⁻¹
    fᶻⁿ⁻¹ = auxiliary_fields.fᶻⁿ⁻¹
    fˣⁿ⁻² = auxiliary_fields.fˣⁿ⁻²
    fʸⁿ⁻² = auxiliary_fields.fʸⁿ⁻²
    fᶻⁿ⁻² = auxiliary_fields.fᶻⁿ⁻²

    Uⁿ⁻¹ = auxiliary_fields.Uⁿ⁻¹
    Vⁿ⁻¹ = auxiliary_fields.Vⁿ⁻¹
    Wⁿ⁻¹ = auxiliary_fields.Wⁿ⁻¹
    bⁿ⁻¹ = auxiliary_fields.bⁿ⁻¹

    @inbounds begin
        # Save previous advective fluxes
        fˣⁿ⁻²[i, j, k] = fˣⁿ⁻¹[i, j, k] 
        fʸⁿ⁻²[i, j, k] = fʸⁿ⁻¹[i, j, k] 
        fᶻⁿ⁻²[i, j, k] = fᶻⁿ⁻¹[i, j, k] 
        
        # Save previous transport and previous buoyancy
        Uⁿ⁻¹[i, j, k] = u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid)
        Vⁿ⁻¹[i, j, k] = v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)
        Wⁿ⁻¹[i, j, k] = w[i, j, k] * Azᶜᶜᶠ(i, j, k, grid)
        bⁿ⁻¹[i, j, k] = b[i, j, k]
    end

    @inbounds begin
        # Calculate new advective fluxes
        fˣⁿ⁻¹[i, j, k] = _advective_tracer_flux_x(i, j, k, grid, advection, u, b) * vertical_scaling(i, j, k, grid, c, c, c)
        fʸⁿ⁻¹[i, j, k] = _advective_tracer_flux_y(i, j, k, grid, advection, v, b) * vertical_scaling(i, j, k, grid, c, c, c)
        fᶻⁿ⁻¹[i, j, k] = _advective_tracer_flux_z(i, j, k, grid, advection, w, b) * vertical_scaling(i, j, k, grid, c, c, c) 
    end
end

function assemble_P_values!(simulation)
    model = simulation.model
    grid = model.grid
    arch = architecture(grid)

    b = model.tracers.b
    bⁿ⁻¹ = model.auxiliary_fields.bⁿ⁻¹
    Uⁿ⁻¹ = model.auxiliary_fields.Uⁿ⁻¹
    Vⁿ⁻¹ = model.auxiliary_fields.Vⁿ⁻¹
    Wⁿ⁻¹ = model.auxiliary_fields.Wⁿ⁻¹

    Pu   = model.auxiliary_fields.Pu
    Pv   = model.auxiliary_fields.Pv
    Pw   = model.auxiliary_fields.Pw
    ∂xb² = model.auxiliary_fields.∂xb²
    ∂yb² = model.auxiliary_fields.∂yb²
    ∂zb² = model.auxiliary_fields.∂zb²

    fˣⁿ⁻¹ = simulation.model.auxiliary_fields.fˣⁿ⁻¹
    fʸⁿ⁻¹ = simulation.model.auxiliary_fields.fʸⁿ⁻¹
    fᶻⁿ⁻¹ = simulation.model.auxiliary_fields.fᶻⁿ⁻¹

    fˣⁿ⁻² = simulation.model.auxiliary_fields.fˣⁿ⁻²
    fʸⁿ⁻² = simulation.model.auxiliary_fields.fʸⁿ⁻²
    fᶻⁿ⁻² = simulation.model.auxiliary_fields.fᶻⁿ⁻²

    χ = simulation.model.timestepper.χ

    launch!(arch, grid, :xyz, _compute_dissipation!, 
            Pu, Pv, Pw, 
            ∂xb², ∂yb², ∂zb², 
            grid, χ, 
            Uⁿ⁻¹, Vⁿ⁻¹, Wⁿ⁻¹, 
            fˣⁿ⁻¹, fʸⁿ⁻¹, fᶻⁿ⁻¹, 
            fˣⁿ⁻², fʸⁿ⁻², fᶻⁿ⁻²,
            b, bⁿ⁻¹)

    return nothing
end

@kernel function _compute_dissipation!(Pu, Pv, Pw, 
                                       ∂xb², ∂yb², ∂zb², 
                                       grid, χ, 
                                       Uⁿ⁻¹, Vⁿ⁻¹, Wⁿ⁻¹, 
                                       fˣⁿ⁻¹, fʸⁿ⁻¹, fᶻⁿ⁻¹, 
                                       fˣⁿ⁻², fʸⁿ⁻², fᶻⁿ⁻²,
                                       b, bⁿ⁻¹)
    
    i, j, k = @index(Global, NTuple)

    @inbounds Pu[i, j, k] = compute_Pᵁ(i, j, k, grid, Uⁿ⁻¹, χ, fˣⁿ⁻¹, fˣⁿ⁻², b, bⁿ⁻¹)
    @inbounds Pv[i, j, k] = compute_Pⱽ(i, j, k, grid, Vⁿ⁻¹, χ, fʸⁿ⁻¹, fʸⁿ⁻², b, bⁿ⁻¹)
    @inbounds Pw[i, j, k] = compute_Pᵂ(i, j, k, grid, Wⁿ⁻¹, χ, fᶻⁿ⁻¹, fᶻⁿ⁻², b, bⁿ⁻¹)

    @inbounds ∂xb²[i, j, k] = Axᶠᶜᶜ(i, j, k, grid) * δxᶠᶜᶜ(i, j, k, grid, b)^2 / Δxᶠᶜᶜ(i, j, k, grid)
    @inbounds ∂yb²[i, j, k] = Ayᶜᶠᶜ(i, j, k, grid) * δyᶜᶠᶜ(i, j, k, grid, b)^2 / Δyᶜᶠᶜ(i, j, k, grid)
    @inbounds ∂zb²[i, j, k] = Azᶜᶜᶠ(i, j, k, grid) * δzᶜᶜᶠ(i, j, k, grid, b)^2 / Δzᶜᶜᶠ(i, j, k, grid)
end

@inline b★(i, j, k, grid, bⁿ, bⁿ⁻¹) = @inbounds (bⁿ[i, j, k] + bⁿ⁻¹[i, j, k]) / 2
@inline b²(i, j, k, grid, b₁, b₂)   = @inbounds (b₁[i, j, k] * b₂[i, j, k])

@inline function compute_Pᵁ(i, j, k, grid, U, χ, fˣⁿ⁻¹, fˣⁿ⁻², bⁿ, bⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δˣb★ = δxᶠᶜᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δˣb² = δxᶠᶜᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        Fⁿ⁻¹ = C₁ * fˣⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fˣⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜x = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟x = U[i, j, k] * δˣb²
    end

    return 2 * δˣb★ * 𝒜x - 𝒟x
end

@inline function compute_Pⱽ(i, j, k, grid, V, χ, fʸⁿ⁻¹, fʸⁿ⁻², bⁿ, bⁿ⁻¹)

    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δʸb★ = δyᶜᶠᶜ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δʸb² = δyᶜᶠᶜ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        Fⁿ⁻¹ = C₁ * fʸⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fʸⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜y = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟y = V[i, j, k] * δʸb²
    end

    return 2 * δʸb★ * 𝒜y - 𝒟y
end

@inline function compute_Pᵂ(i, j, k, grid, W, χ, fᶻⁿ⁻¹, fᶻⁿ⁻², bⁿ, bⁿ⁻¹)
   
    C₁ = convert(eltype(grid), 1.5 + χ)
    C₂ = convert(eltype(grid), 0.5 + χ)

    δᶻb★ = δzᶜᶜᶠ(i, j, k, grid, b★, bⁿ, bⁿ⁻¹)
    δᶻb² = δzᶜᶜᶠ(i, j, k, grid, b², bⁿ, bⁿ⁻¹)

    @inbounds begin
        Fⁿ⁻¹ = C₁ * fᶻⁿ⁻¹[i, j, k] / vertical_scaling(i, j, k, grid, c, c, c)
        Fⁿ⁻² = C₂ * fᶻⁿ⁻²[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, c)
        𝒜y = Fⁿ⁻¹ - Fⁿ⁻²
        𝒟y = W[i, j, k] * δᶻb²
    end

    return 2 * δᶻb★ * 𝒜y - 𝒟y
end

