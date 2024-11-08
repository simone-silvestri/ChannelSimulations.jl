@kernel function _assemble_diffusive_tracer_dissipation!(K, grid, χ, Vⁿ, Vⁿ⁻¹, Uⁿ⁺¹, cⁿ⁺¹, cⁿ)
    i, j, k = @index(Global, NTuple)
    compute_diffusive_tracer_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, cⁿ⁺¹, cⁿ)
end

@inline function compute_diffusive_tracer_dissipation!(K::Tuple, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, cⁿ⁺¹, cⁿ)
    for n in eachindex(K)
        compute_diffusive_tracer_dissipation!(K[n], i, j, k, grid, Vⁿ[n], Vⁿ⁻¹[n], χ, cⁿ⁺¹, cⁿ)
    end
end

@inline compute_diffusive_tracer_dissipation!(::Nothing, args...) = nothing

@inline function compute_diffusive_tracer_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, cⁿ⁺¹, cⁿ)
    C₁  = 3//2 + χ
    C₂  = 1//2 + χ

    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)

    @inbounds begin
        fx₁ = C₁ * Vⁿ.x[i, j, k] / vertical_scaling(i, j, k, grid, f, c, c)
        fy₁ = C₁ * Vⁿ.y[i, j, k] / vertical_scaling(i, j, k, grid, c, f, c)
        fz₁ = C₁ * Vⁿ.z[i, j, k] / vertical_scaling(i, j, k, grid, c, c, f)

        fx₂ = C₂ * Vⁿ⁻¹.x[i, j, k] / previous_vertical_scaling(i, j, k, grid, f, c, c)
        fy₂ = C₂ * Vⁿ⁻¹.y[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, f, c)
        fz₂ = C₂ * Vⁿ⁻¹.z[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, f)    
    
        K.x[i, j, k] = 2 * δˣc★ * (fx₁ - fx₂)
        K.y[i, j, k] = 2 * δʸc★ * (fy₁ - fy₂)
        K.z[i, j, k] = 2 * δᶻc★ * (fz₁ - fz₂)
    end
end

@kernel function _assemble_diffusive_vorticity_dissipation!(K, grid, χ, Vⁿ, Vⁿ⁻¹, Uⁿ⁺¹, cⁿ⁺¹, ζⁿ)
    i, j, k = @index(Global, NTuple)
    compute_diffusive_vorticity_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, Uⁿ⁺¹, ζⁿ)
end

@inline function compute_diffusive_vorticity_dissipation!(K::Tuple, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, Uⁿ⁺¹, ζⁿ)
    for n in eachindex(K)
        compute_diffusive_vorticity_dissipation!(K[n], i, j, k, grid, Vⁿ[n], Vⁿ⁻¹[n], χ, Uⁿ⁺¹, ζⁿ)
    end
end

@inline compute_diffusive_vorticity_dissipation!(::Nothing, args...) = nothing

@inline function compute_diffusive_vorticity_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, Uⁿ⁺¹, ζⁿ)
    C₁  = convert(eltype(grid), 1.5 + χ)
    C₂  = convert(eltype(grid), 0.5 + χ)
    
    δˣζ★ = δxᶠᶜᶜ(i, j, k, grid, ζ★, Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)
    δʸζ★ = δyᶜᶠᶜ(i, j, k, grid, ζ★, Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)
    
    @inbounds begin
        K.x[i, j, k] = 2 * δˣζ★ * (C₁ * Vⁿ.x[i, j, k] - C₂ * Vⁿ⁻¹.x[i, j, k])
        K.y[i, j, k] = 2 * δʸζ★ * (C₁ * Vⁿ.y[i, j, k] - C₂ * Vⁿ⁻¹.y[i, j, k])
    end
end