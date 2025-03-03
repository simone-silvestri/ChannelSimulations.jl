# TODO: This is only for AB2, figure out how to generalize this for other timesteppers for example RK3
@kernel _assemble_advective_tracer_dissipation!(::Nothing, args...) = nothing

@kernel function _assemble_ab2_advective_tracer_dissipation!(P, grid, χ::FT, Fⁿ, Fⁿ⁻¹, Uⁿ⁺¹, Uⁿ, Uⁿ⁻¹, cⁿ⁺¹, cⁿ) where FT
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
        u₁ = C₁ * Uⁿ.u[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        v₁ = C₁ * Uⁿ.v[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        w₁ = C₁ * Uⁿ.w[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        u₂ = C₂ * Uⁿ⁻¹.u[i, j, k] / previous_vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        v₂ = C₂ * Uⁿ⁻¹.v[i, j, k] / previous_vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        w₂ = C₂ * Uⁿ⁻¹.w[i, j, k] / previous_vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        fx₁ = C₁ * Fⁿ.x[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        fy₁ = C₁ * Fⁿ.y[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        fz₁ = C₁ * Fⁿ.z[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        fx₂ = C₂ * Fⁿ⁻¹.x[i, j, k] / previous_vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        fy₂ = C₂ * Fⁿ⁻¹.y[i, j, k] / previous_vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        fz₂ = C₂ * Fⁿ⁻¹.z[i, j, k] / previous_vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        P.x[i, j, k] = 2 * δˣc★ * (fx₁ - fx₂) - δˣc² * (u₁ - u₂)
        P.y[i, j, k] = 2 * δʸc★ * (fy₁ - fy₂) - δʸc² * (v₁ - v₂)
        P.z[i, j, k] = 2 * δᶻc★ * (fz₁ - fz₂) - δᶻc² * (w₁ - w₂)
    end
end

@kernel function _assemble_rk3_advective_tracer_dissipation!(P, grid, γ, Fⁿ, Uⁿ, cⁿ⁺¹, cⁿ) 
    i, j, k = @index(Global, NTuple)

    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δˣc² = δxᶠᶜᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δʸc² = δyᶜᶠᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc² = δzᶜᶜᶠ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)

    @inbounds begin
        u₁ = Uⁿ.u[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        v₁ = Uⁿ.v[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        w₁ = Uⁿ.w[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        fx₁ = Fⁿ.x[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        fy₁ = Fⁿ.y[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        fz₁ = Fⁿ.z[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        P.x[i, j, k] += 2 * δˣc★ * fx₁ - δˣc² * u₁
        P.y[i, j, k] += 2 * δʸc★ * fy₁ - δʸc² * v₁
        P.z[i, j, k] += 2 * δᶻc★ * fz₁ - δᶻc² * w₁

        P.x[i, j, k] *= γ
        P.y[i, j, k] *= γ
        P.z[i, j, k] *= γ
    end
end

@kernel function _assemble_rk3_advective_tracer_dissipation!(P, grid, ::Nothing, Fⁿ, Uⁿ, cⁿ⁺¹, cⁿ) 
    i, j, k = @index(Global, NTuple)

    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δˣc² = δxᶠᶜᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δʸc² = δyᶜᶠᶜ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)
    
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc² = δzᶜᶜᶠ(i, j, k, grid, c², cⁿ⁺¹, cⁿ)

    @inbounds begin
        u₁ = Uⁿ.u[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        v₁ = Uⁿ.v[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        w₁ = Uⁿ.w[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        fx₁ = Fⁿ.x[i, j, k] / vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        fy₁ = Fⁿ.y[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        fz₁ = Fⁿ.z[i, j, k] / vertical_scaling(i, j, k, grid, Center(), Center(), Face())

        P.x[i, j, k] = 2 * δˣc★ * fx₁ - δˣc² * u₁
        P.y[i, j, k] = 2 * δʸc★ * fy₁ - δʸc² * v₁
        P.z[i, j, k] = 2 * δᶻc★ * fz₁ - δᶻc² * w₁
    end
end