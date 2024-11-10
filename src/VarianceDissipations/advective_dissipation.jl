# TODO: This is only for AB2, figure out how to generalize this for other timesteppers for example RK3
@kernel _assemble_advective_tracer_dissipation!(::Nothing, args...) = nothing

@kernel function _assemble_advective_tracer_dissipation!(P, grid, χ::FT, Fⁿ, Fⁿ⁻¹, Uⁿ⁺¹, Uⁿ, Uⁿ⁻¹, cⁿ⁺¹, cⁿ) where FT
    i, j, k = @index(Global, NTuple)

    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δˣc² = δxᶠᶜᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δʸc² = δyᶜᶠᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc² = δzᶜᶜᶠ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    C₁  = convert(FT, 3/2) + χ
    C₂  = convert(FT, 1/2) + χ

    @inbounds begin
        u₁ = C₁ * Uⁿ.u[i, j, k] / vertical_scaling(i, j, k, grid, f, c, c)
        v₁ = C₁ * Uⁿ.v[i, j, k] / vertical_scaling(i, j, k, grid, c, f, c)
        w₁ = C₁ * Uⁿ.w[i, j, k] / vertical_scaling(i, j, k, grid, c, c, f)

        u₂ = C₂ * Uⁿ⁻¹.u[i, j, k] / previous_vertical_scaling(i, j, k, grid, f, c, c)
        v₂ = C₂ * Uⁿ⁻¹.v[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, f, c)
        w₂ = C₂ * Uⁿ⁻¹.w[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, f)

        fx₁ = C₁ * Fⁿ.x[i, j, k] / vertical_scaling(i, j, k, grid, f, c, c)
        fy₁ = C₁ * Fⁿ.y[i, j, k] / vertical_scaling(i, j, k, grid, c, f, c)
        fz₁ = C₁ * Fⁿ.z[i, j, k] / vertical_scaling(i, j, k, grid, c, c, f)

        fx₂ = C₂ * Fⁿ⁻¹.x[i, j, k] / previous_vertical_scaling(i, j, k, grid, f, c, c)
        fy₂ = C₂ * Fⁿ⁻¹.y[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, f, c)
        fz₂ = C₂ * Fⁿ⁻¹.z[i, j, k] / previous_vertical_scaling(i, j, k, grid, c, c, f)

        P.x[i, j, k] = 2 * δˣc★ * (fx₁ - fx₂) - δˣc² * (u₁ - u₂)
        P.y[i, j, k] = 2 * δʸc★ * (fy₁ - fy₂) - δʸc² * (v₁ - v₂)
        P.z[i, j, k] = 2 * δᶻc★ * (fz₁ - fz₂) - δᶻc² * (w₁ - w₂)
    end
end

@kernel _assemble_advective_vorticity_dissipation!(::Nothing, args...) = nothing

@kernel function _assemble_advective_vorticity_dissipation!(P, grid, χ::FT, Fⁿ, Fⁿ⁻¹, Uⁿ⁺¹, Uⁿ, Uⁿ⁻¹, c, ζⁿ) where FT
    i, j, k = @index(Global, NTuple)

    δˣζ★ = δxᶠᶜᶜ(i, j, k, grid, ζ★, Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)
    δˣζ² = δxᶠᶜᶜ(i, j, k, grid, ζ², Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)

    δʸζ★ = δyᶜᶠᶜ(i, j, k, grid, ζ★, Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)
    δʸζ² = δyᶜᶠᶜ(i, j, k, grid, ζ², Uⁿ⁺¹.u, Uⁿ⁺¹.v, ζⁿ)

    C₁  = convert(FT, 3/2) + χ
    C₂  = convert(FT, 1/2) + χ

    @inbounds begin
        u₁ = C₁ * ℑxyᶜᶠᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, Uⁿ.u)   / Δyᶠᶜᶜ(i, j, k, grid)
        u₂ = C₂ * ℑxyᶜᶠᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, Uⁿ⁻¹.u) / Δyᶠᶜᶜ(i, j, k, grid)
        v₁ = C₁ * ℑxyᶠᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, Uⁿ.v)   / Δxᶜᶠᶜ(i, j, k, grid)
        v₂ = C₂ * ℑxyᶠᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, Uⁿ⁻¹.v) / Δxᶜᶠᶜ(i, j, k, grid)

        fx₁ = C₁ * Fⁿ.x[i, j, k]
        fy₁ = C₁ * Fⁿ.y[i, j, k]

        fx₂ = C₂ * Fⁿ⁻¹.x[i, j, k]
        fy₂ = C₂ * Fⁿ⁻¹.y[i, j, k]

        P.x[i, j, k] = 2 * δˣζ★ * (fx₁ - fx₂) - δˣζ² * (u₁ - u₂)
        P.y[i, j, k] = 2 * δʸζ★ * (fy₁ - fy₂) - δʸζ² * (v₁ - v₂)
    end
end
