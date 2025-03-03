@kernel function _update_diffusive_tracer_fluxes!(Vⁿ, Vⁿ⁻¹, grid, closure, diffusivity, bouyancy, cⁿ, tracer_id, clk, model_fields) 
    i, j, k = @index(Global, NTuple)
    compute_diffusive_tracer_fluxes!(Vⁿ, Vⁿ⁻¹, i, j, k, grid, closure, diffusivity, bouyancy, cⁿ, tracer_id, clk, model_fields)
end

@inline function compute_diffusive_tracer_fluxes!(Vⁿ, Vⁿ⁻¹, i, j, k, grid, closure::Tuple, K, args...) 
    for n in eachindex(closure)
        compute_diffusive_tracer_fluxes!(Vⁿ[n], Vⁿ⁻¹[n], i, j, k, grid, closure[n], K[n], args...)
    end
end

@inline compute_diffusive_tracer_fluxes!(::Nothing, args...) = nothing

@inline function compute_diffusive_tracer_fluxes!(Vⁿ, Vⁿ⁻¹, i, j, k, grid, clo, K, b, cⁿ, c_id, clk, fields)
    @inbounds begin
        Vⁿ⁻¹.x[i, j, k] = Vⁿ.x[i, j, k]
        Vⁿ⁻¹.y[i, j, k] = Vⁿ.y[i, j, k]
        Vⁿ⁻¹.z[i, j, k] = Vⁿ.z[i, j, k]

        ETD = ExplicitTimeDiscretization()

        Vⁿ.x[i, j, k] = diffusive_flux_x(i, j, k, grid, ETD, clo, K, c_id, cⁿ, clk, fields, b) * Axᶠᶜᶜ(i, j, k, grid) * σⁿ(i, j, k, grid, f, c, c)
        Vⁿ.y[i, j, k] = diffusive_flux_y(i, j, k, grid, ETD, clo, K, c_id, cⁿ, clk, fields, b) * Ayᶜᶠᶜ(i, j, k, grid) * σⁿ(i, j, k, grid, c, f, c)
        Vⁿ.z[i, j, k] = diffusive_flux_z(i, j, k, grid, ETD, clo, K, c_id, cⁿ, clk, fields, b) * Azᶜᶜᶠ(i, j, k, grid) * σⁿ(i, j, k, grid, c, c, f)
    end
end