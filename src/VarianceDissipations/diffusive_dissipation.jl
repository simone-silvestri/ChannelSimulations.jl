@kernel function _assemble_diffusive_ab2_tracer_dissipation!(K, grid, χ, Vⁿ, Vⁿ⁻¹, Uⁿ⁺¹, cⁿ⁺¹, cⁿ)
    i, j, k = @index(Global, NTuple)
    compute_diffusive_ab2_tracer_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, cⁿ⁺¹, cⁿ)
end

@inline function compute_diffusive_ab2_tracer_dissipation!(K::Tuple, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ, cⁿ⁺¹, cⁿ)
    for n in eachindex(K)
        compute_diffusive_ab2_tracer_dissipation!(K[n], i, j, k, grid, Vⁿ[n], Vⁿ⁻¹[n], χ, cⁿ⁺¹, cⁿ)
    end
end

@inline compute_diffusive_ab2_tracer_dissipation!(::Nothing, args...) = nothing

@inline function compute_diffusive_ab2_tracer_dissipation!(K, i, j, k, grid, Vⁿ, Vⁿ⁻¹, χ::FT, cⁿ⁺¹, cⁿ) where FT
    C₁  = convert(FT, 3/2) + χ
    C₂  = convert(FT, 1/2) + χ

    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)

    @inbounds begin
        fx₁ = C₁ * Vⁿ.x[i, j, k] / σⁿ(i, j, k, grid, f, c, c)
        fy₁ = C₁ * Vⁿ.y[i, j, k] / σⁿ(i, j, k, grid, c, f, c)
        fz₁ = C₁ * Vⁿ.z[i, j, k] / σⁿ(i, j, k, grid, c, c, f)

        fx₂ = C₂ * Vⁿ⁻¹.x[i, j, k] / σ⁻(i, j, k, grid, f, c, c)
        fy₂ = C₂ * Vⁿ⁻¹.y[i, j, k] / σ⁻(i, j, k, grid, c, f, c)
        fz₂ = C₂ * Vⁿ⁻¹.z[i, j, k] / σ⁻(i, j, k, grid, c, c, f)
    
        K.x[i, j, k] = 2 * δˣc★ * (fx₁ - fx₂)
        K.y[i, j, k] = 2 * δʸc★ * (fy₁ - fy₂)
        K.z[i, j, k] = 2 * δᶻc★ * (fz₁ - fz₂)
    end
end

@kernel function _assemble_diffusive_rk3_tracer_dissipation!(K, grid, γ, Vⁿ, cⁿ⁺¹, cⁿ)
    i, j, k = @index(Global, NTuple)
    compute_diffusive_rk3_tracer_dissipation!(K, i, j, k, grid, Vⁿ, γ, cⁿ⁺¹, cⁿ)
end

@inline function compute_diffusive_rk3_tracer_dissipation!(K::Tuple, i, j, k, grid, Vⁿ, γ, cⁿ⁺¹, cⁿ)
    for n in eachindex(K)
        compute_diffusive_ab2_tracer_dissipation!(K[n], i, j, k, grid, Vⁿ[n], γ, cⁿ⁺¹, cⁿ)
    end
end

@inline compute_diffusive_rk3_tracer_dissipation!(::Nothing, args...) = nothing
@inline compute_diffusive_rk3_tracer_dissipation!(::Nothing, i, j, k, grid, Vⁿ, ::Nothing, cⁿ⁺¹, cⁿ) = nothing

@inline function compute_diffusive_rk3_tracer_dissipation!(K, i, j, k, grid, Vⁿ, ::Nothing, cⁿ⁺¹, cⁿ)
    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)

    @inbounds begin
        fx₁ = Vⁿ.x[i, j, k] / σⁿ(i, j, k, grid, f, c, c)
        fy₁ = Vⁿ.y[i, j, k] / σⁿ(i, j, k, grid, c, f, c)
        fz₁ = Vⁿ.z[i, j, k] / σⁿ(i, j, k, grid, c, c, f)

        K.x[i, j, k] = 2 * δˣc★ * fx₁ 
        K.y[i, j, k] = 2 * δʸc★ * fy₁ 
        K.z[i, j, k] = 2 * δᶻc★ * fz₁ 
    end
end

@inline function compute_diffusive_rk3_tracer_dissipation!(K, i, j, k, grid, Vⁿ, γ, cⁿ⁺¹, cⁿ)
    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)    
    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)
    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ⁺¹, cⁿ)

    @inbounds begin
        fx₁ = Vⁿ.x[i, j, k] / σⁿ(i, j, k, grid, f, c, c)
        fy₁ = Vⁿ.y[i, j, k] / σⁿ(i, j, k, grid, c, f, c)
        fz₁ = Vⁿ.z[i, j, k] / σⁿ(i, j, k, grid, c, c, f)

        K.x[i, j, k] += 2 * δˣc★ * fx₁ 
        K.y[i, j, k] += 2 * δʸc★ * fy₁ 
        K.z[i, j, k] += 2 * δᶻc★ * fz₁ 

        K.x[i, j, k] *= γ
        K.y[i, j, k] *= γ
        K.z[i, j, k] *= γ
    end
end