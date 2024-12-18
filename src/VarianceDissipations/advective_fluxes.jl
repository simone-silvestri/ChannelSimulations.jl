@kernel function _update_advective_tracer_fluxes!(Gⁿ, Fⁿ, Fⁿ⁻¹, cⁿ⁻¹, grid, advection, U, cⁿ)
    i, j, k = @index(Global, NTuple)
    
    @inbounds begin
        # Save previous advective fluxes
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        Fⁿ⁻¹.z[i, j, k] = Fⁿ.z[i, j, k]
        
        cⁿ⁻¹[i, j, k] = cⁿ[i, j, k]

        # Calculate new advective fluxes
        Fⁿ.x[i, j, k] = _advective_tracer_flux_x(i, j, k, grid, advection, U.u, cⁿ) * vertical_scaling(i, j, k, grid, Face(), Center(), Center())
        Fⁿ.y[i, j, k] = _advective_tracer_flux_y(i, j, k, grid, advection, U.v, cⁿ) * vertical_scaling(i, j, k, grid, Center(), Face(), Center())
        Fⁿ.z[i, j, k] = _advective_tracer_flux_z(i, j, k, grid, advection, U.w, cⁿ) * vertical_scaling(i, j, k, grid, Center(), Center(), Face())
        
        Gⁿ.x[i, j, k] = Axᶠᶜᶜ(i, j, k, grid) * δxᶠᶜᶜ(i, j, k, grid, cⁿ)^2 / Δxᶠᶜᶜ(i, j, k, grid)
        Gⁿ.y[i, j, k] = Ayᶜᶠᶜ(i, j, k, grid) * δyᶜᶠᶜ(i, j, k, grid, cⁿ)^2 / Δyᶜᶠᶜ(i, j, k, grid)
        Gⁿ.z[i, j, k] = Azᶜᶜᶠ(i, j, k, grid) * δzᶜᶜᶠ(i, j, k, grid, cⁿ)^2 / Δzᶜᶜᶠ(i, j, k, grid)
    end
end

@kernel function _update_advective_vorticity_fluxes!(Gⁿ, Fⁿ, Fⁿ⁻¹, ζⁿ⁻¹, grid, advection, U, cⁿ)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        # Save previous advective fluxes
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        
        ζⁿ⁻¹[i, j, k] = ζ₃ᶠᶠᶜ(i, j, k, grid, U.u, U.v)

        # Calculate new advective fluxes
        Fⁿ.x[i, j, k] =   horizontal_advection_V(i, j, k, grid, advection, U.u, U.v) * Axᶜᶠᶜ(i, j, k, grid) 
        Fⁿ.y[i, j, k] = - horizontal_advection_U(i, j, k, grid, advection, U.u, U.v) * Ayᶠᶜᶜ(i, j, k, grid)

        Gⁿ.x[i, j, k] = Axᶜᶠᶜ(i, j, k, grid) * δxᶜᶠᶜ(i, j, k, grid, ζ₃ᶠᶠᶜ, U.u, U.v)^2 / Δxᶜᶠᶜ(i, j, k, grid)
        Gⁿ.y[i, j, k] = Ayᶠᶜᶜ(i, j, k, grid) * δyᶠᶜᶜ(i, j, k, grid, ζ₃ᶠᶠᶜ, U.u, U.v)^2 / Δyᶠᶜᶜ(i, j, k, grid)
    end
end