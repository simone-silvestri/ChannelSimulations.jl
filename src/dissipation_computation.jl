using Oceananigans: instantiated_location
using Oceananigans.Grids: architecture
using Oceananigans.Utils
using Oceananigans.TimeSteppers
using Oceananigans.Fields
using Oceananigans.Fields: Field, VelocityFields
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y, _advective_tracer_flux_z
using Oceananigans.Operators: volume
using KernelAbstractions: @kernel, @index

import Oceananigans.Utils: KernelParameters

struct TracerVarianceDissipation{P, K, A, D, S, G}
    advective_production :: P
    diffusive_production :: K
    advective_fluxes :: A
    diffusive_fluxes :: D
    previous_state :: S
    gradient_squared :: G
end

# Function to call in a callback
# Note: This works only if the callback is called with an IterationInterval(1), if not the
# previous fluxes and velocities will not be correct
# TODO: make sure that the correct velocities and fluxes are used even if 
# the callback is not called with an IterationInterval(1)
function (dc::TracerVarianceDissipation)(simulation)
    # We first assemble values for Pⁿ⁻¹
    assemble_P_values!(simulation, dc)

    # Then we update the fluxes to be used in the next time step
    update_fluxes!(simulation, dc)

    return nothing
end

@inline function KernelParameters(f::Field)
    sz = size(f.data)
    of = f.data.offsets
    return KernelParameters(sz, of)
end

function fluxes_fields(grid)
    x = XFaceField(grid)
    y = YFaceField(grid)
    z = ZFaceField(grid)

    return (; x, y, z)
end

function get_dissipation_fields(t::TracerVarianceDissipation)
    f = NamedTuple()
    
    for name in keys(t.advective_production)
        f = merge(f, get_dissipation_fields(t, name))
    end

    return f
end

function get_dissipation_fields(t::TracerVarianceDissipation, tracer_name) 
    P = t.advective_production[tracer_name]
    G = t.gradient_squared[tracer_name]

    dirs = tracer_name == :ζ ? (:x, :y, :z) : (:x, :y)

    prod_names = Tuple(Symbol(:P, tracer_name, dir) for dir in dirs)
    grad_names = Tuple(Symbol(:G, tracer_name, dir) for dir in dirs)

    prod = Tuple(getproperty(P, dir) for dir in dirs)
    grad = Tuple(getproperty(G, dir) for dir in dirs)

    return NamedTuple{tuple(prod_names..., grad_names...)}(tuple(prod..., grad...))
end

function TracerVarianceDissipation(model; tracers = propertynames(model.tracers), include_vorticity = true)
        
    if !(model.timestepper isa QuasiAdamsBashforth2TimeStepper)
        throw(ArgumentError("DissipationComputation requires a QuasiAdamsBashforth2TimeStepper"))
    end

    tracers = tupleit(tracers)

    grid = model.grid
    P    = NamedTuple{tracers}(fluxes_fields(grid) for tracer in tracers)
    Fⁿ   = NamedTuple{tracers}(fluxes_fields(grid) for tracer in tracers)
    Fⁿ⁻¹ = NamedTuple{tracers}(fluxes_fields(grid) for tracer in tracers)
    Uⁿ⁻¹ = VelocityFields(grid)
    Uⁿ   = VelocityFields(grid)
    
    cⁿ⁻¹ =  NamedTuple{tracers}(CenterField(grid) for tracer in tracers)

    if include_vorticity
        Fζⁿ   = (; x = YFaceField(grid), y = XFaceField(grid))
        Fζⁿ⁻¹ = (; x = YFaceField(grid), y = XFaceField(grid))
        Pζ    = (; x = ZFaceField(grid), y = ZFaceField(grid))
        ζⁿ⁻¹  = Field{Face, Face, Center}(grid)

        P    = merge(P,    (; ζ = Pζ))
        Fⁿ   = merge(Fⁿ,   (; ζ = Fζⁿ))
        Fⁿ⁻¹ = merge(Fⁿ⁻¹, (; ζ = Fζⁿ⁻¹))
        cⁿ⁻¹ = merge(cⁿ⁻¹, (; ζ = ζⁿ⁻¹))
    end

    previous_state = merge(cⁿ⁻¹, (; Uⁿ⁻¹, Uⁿ))
    advective_fluxes = (; Fⁿ, Fⁿ⁻¹)

    gradients = deepcopy(P)

    return TracerVarianceDissipation(P, advective_fluxes, previous_state, gradients)
end

@inline getadvection(advection, tracer_name) = advection
@inline getadvection(advection::NamedTuple, tracer_name) = @inbounds ifelse(tracer_name == :ζ, advection.momentum, advection[tracer_name])

function update_fluxes!(simulation, dissipation_computation)
    model = simulation.model

    grid = model.grid
    arch = architecture(grid)

    params = KernelParameters(model.tracers[1])
    
    Uⁿ   = dissipation_computation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation_computation.previous_state.Uⁿ⁻¹ 

    U = model.velocities

    launch!(architecture(grid), grid, params, _update_transport!, Uⁿ, Uⁿ⁻¹, grid, U)

    for tracer_name in keys(dissipation_computation.advective_production)
        if tracer_name != :ζ
            c    = model.tracers[tracer_name]
            cⁿ⁻¹ = dissipation_computation.previous_state[tracer_name]
            Fⁿ   = dissipation_computation.advective_fluxes.Fⁿ[tracer_name]
            Fⁿ⁻¹ = dissipation_computation.advective_fluxes.Fⁿ⁻¹[tracer_name]
            A    = getadvection(model.advection, tracer_name)

            launch!(arch, grid, params, _update_fluxes!, Fⁿ, Fⁿ⁻¹, cⁿ⁻¹, grid, A, U, c)
        end
    end

    if :ζ ∈ keys(dissipation_computation.advective_production)
        ζⁿ⁻¹     = dissipation_computation.previous_state[:ζ]
        Fⁿ       = dissipation_computation.advective_fluxes.Fⁿ[:ζ]
        Fⁿ⁻¹     = dissipation_computation.advective_fluxes.Fⁿ⁻¹[:ζ]
        A        = model.advection.momentum

        launch!(arch, grid, params, _update_vorticity_fluxes!, Fⁿ, Fⁿ⁻¹, ζⁿ⁻¹, grid, A, U)
    end

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

@kernel function _update_fluxes!( Fⁿ, Fⁿ⁻¹, cⁿ⁻¹, grid, advection, U, c)
    i, j, k = @index(Global, NTuple)
    u, v, w = U

    @inbounds begin
        # Save previous advective fluxes
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        Fⁿ⁻¹.z[i, j, k] = Fⁿ.z[i, j, k]
        
        cⁿ⁻¹[i, j, k] = c[i, j, k]

        # Calculate new advective fluxes
        Fⁿ.x[i, j, k] = _advective_tracer_flux_x(i, j, k, grid, advection, u, c) 
        Fⁿ.y[i, j, k] = _advective_tracer_flux_y(i, j, k, grid, advection, v, c) 
        Fⁿ.z[i, j, k] = _advective_tracer_flux_z(i, j, k, grid, advection, w, c) 
    end
end

@kernel function _update_vorticity_fluxes!( Fⁿ, Fⁿ⁻¹, ζⁿ⁻¹, grid, advection, U)
    i, j, k = @index(Global, NTuple)
    u, v, w = U

    @inbounds begin
        # Save previous advective fluxes
        Fⁿ⁻¹.x[i, j, k] = Fⁿ.x[i, j, k]
        Fⁿ⁻¹.y[i, j, k] = Fⁿ.y[i, j, k]
        
        ζⁿ⁻¹[i, j, k] = ζ₃ᶠᶠᶜ(i, j, k, grid, u, v)

        # Calculate new advective fluxes
        Fⁿ.x[i, j, k] =   horizontal_advection_V(i, j, k, grid, advection, u, v) * Axᶜᶠᶜ(i, j, k, grid)
        Fⁿ.y[i, j, k] = - horizontal_advection_U(i, j, k, grid, advection, u, v) * Ayᶠᶜᶜ(i, j, k, grid)
    end
end

function assemble_P_values!(simulation, dissipation_computation)
    model = simulation.model
    grid = model.grid
    arch = architecture(grid)

    χ = simulation.model.timestepper.χ

    # General velocities
    Uⁿ   = dissipation_computation.previous_state.Uⁿ
    Uⁿ⁻¹ = dissipation_computation.previous_state.Uⁿ⁻¹

    for tracer_name in keys(dissipation_computation.advective_production)
        if tracer_name != :ζ
            c = model.tracers[tracer_name]
            cⁿ⁻¹ = dissipation_computation.previous_state[tracer_name]
            P = dissipation_computation.advective_production[tracer_name]
            G = dissipation_computation.gradient_squared[tracer_name]
            Fⁿ = dissipation_computation.advective_fluxes.Fⁿ[tracer_name]
            Fⁿ⁻¹ = dissipation_computation.advective_fluxes.Fⁿ⁻¹[tracer_name]

            launch!(arch, grid, :xyz, _compute_dissipation!, P, G, grid, χ, 
                                                             Fⁿ, Fⁿ⁻¹, 
                                                             Uⁿ, Uⁿ⁻¹, 
                                                             c, cⁿ⁻¹)
        end
    end

    if :ζ ∈ keys(dissipation_computation.advective_production)
        ζⁿ⁻¹ = dissipation_computation.previous_state[:ζ]
        Fⁿ   = dissipation_computation.advective_fluxes.Fⁿ[:ζ]
        Fⁿ⁻¹ = dissipation_computation.advective_fluxes.Fⁿ⁻¹[:ζ]
        P    = dissipation_computation.advective_production[:ζ]
        G    = dissipation_computation.gradient_squared[:ζ]

        launch!(arch, grid, :xyz, _compute_enstrophy_dissipation!, P, grid, χ, 
                                                                   Fⁿ, Fⁿ⁻¹, 
                                                                   Uⁿ, Uⁿ⁻¹, 
                                                                   ζⁿ⁻¹)
    end

    return nothing
end

@inline c★(i, j, k, grid, cⁿ, cⁿ⁻¹) = @inbounds (cⁿ[i, j, k] + cⁿ⁻¹[i, j, k]) / 2
@inline c²(i, j, k, grid, c₁, c₂)   = @inbounds (c₁[i, j, k] * c₂[i, j, k])

@kernel function _compute_dissipation!(P,
                                       grid, χ, 
                                       Fⁿ, Fⁿ⁻¹, 
                                       Uⁿ, Uⁿ⁻¹, 
                                       cⁿ, cⁿ⁻¹)
    
    i, j, k = @index(Global, NTuple)


    δˣc★ = δxᶠᶜᶜ(i, j, k, grid, c★, cⁿ, cⁿ⁻¹)
    δˣc² = δxᶠᶜᶜ(i, j, k, grid, c², cⁿ, cⁿ⁻¹)

    δʸc★ = δyᶜᶠᶜ(i, j, k, grid, c★, cⁿ, cⁿ⁻¹)
    δʸc² = δyᶜᶠᶜ(i, j, k, grid, c², cⁿ, cⁿ⁻¹)

    δᶻc★ = δzᶜᶜᶠ(i, j, k, grid, c★, cⁿ, cⁿ⁻¹)
    δᶻc² = δzᶜᶜᶠ(i, j, k, grid, c², cⁿ, cⁿ⁻¹)

    @inbounds P.x[i, j, k] = compute_dissipation(i, j, k, grid, χ, Fⁿ.x, Fⁿ⁻¹.x, Uⁿ.u, Uⁿ⁻¹.u, δˣc★, δˣc²)
    @inbounds P.y[i, j, k] = compute_dissipation(i, j, k, grid, χ, Fⁿ.y, Fⁿ⁻¹.y, Uⁿ.v, Uⁿ⁻¹.v, δʸc★, δʸc²)
    @inbounds P.z[i, j, k] = compute_dissipation(i, j, k, grid, χ, Fⁿ.z, Fⁿ⁻¹.z, Uⁿ.w, Uⁿ⁻¹.w, δᶻc★, δᶻc²)
end

@inline c★(i, j, k, grid, cⁿ, cⁿ⁻¹) = @inbounds (cⁿ[i, j, k] + cⁿ⁻¹[i, j, k]) / 2
@inline c²(i, j, k, grid, c₁, c₂)   = @inbounds (c₁[i, j, k] * c₂[i, j, k])

@inline ζ★(i, j, k, grid, uⁿ, vⁿ, ζⁿ⁻¹) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ, vⁿ) + ζⁿ⁻¹[i, j, k]) / 2
@inline ζ²(i, j, k, grid, uⁿ, vⁿ, ζⁿ⁻¹) = @inbounds (ζ₃ᶠᶠᶜ(i, j, k, grid, uⁿ, vⁿ) * ζⁿ⁻¹[i, j, k])

@kernel function _compute_enstrophy_dissipation!(P,
                                                 grid, χ, 
                                                 Fⁿ, Fⁿ⁻¹, 
                                                 Uⁿ, Uⁿ⁻¹, 
                                                 ζⁿ⁻¹)
    
    i, j, k = @index(Global, NTuple)

    δˣζ★ = δxᶠᶜᶜ(i, j, k, grid, ζ★, U.uⁿ, U.vⁿ, ζⁿ⁻¹)
    δˣζ² = δxᶠᶜᶜ(i, j, k, grid, ζ², U.uⁿ, U.vⁿ, ζⁿ⁻¹)

    δʸζ★ = δyᶜᶠᶜ(i, j, k, grid, ζ★, U.uⁿ, U.vⁿ, ζⁿ⁻¹)
    δʸζ² = δyᶜᶠᶜ(i, j, k, grid, ζ², U.uⁿ, U.vⁿ, ζⁿ⁻¹)

    @inbounds P.x[i, j, k] = compute_enstrophy_dissipation(i, j, k, grid, χ, Fⁿ.x, Fⁿ⁻¹.x, ℑxyᶜᶠᵃ, Uⁿ.u, Uⁿ⁻¹.u, δˣζ★, δˣζ²)
    @inbounds P.y[i, j, k] = compute_enstrophy_dissipation(i, j, k, grid, χ, Fⁿ.y, Fⁿ⁻¹.y, ℑxyᶜᶠᵃ, Uⁿ.v, Uⁿ⁻¹.v, δʸζ★, δʸζ²)
end

@inline function compute_dissipation(i, j, k, grid, χ, fⁿ, fⁿ⁻¹, Uⁿ, Uⁿ⁻¹, δc★, δc²)

    C₁  = convert(eltype(grid), 1.5 + χ)
    C₂  = convert(eltype(grid), 0.5 + χ)
    loc = instantiated_location(Uⁿ)

    @inbounds begin
        𝒰ⁿ   = C₁ * Uⁿ[i, j, k] 
        𝒰ⁿ⁻¹ = C₂ * Uⁿ⁻¹[i, j, k] 
        Fⁿ   = C₁ * fⁿ[i, j, k] 
        Fⁿ⁻¹ = C₂ * fⁿ⁻¹[i, j, k] 
        A = Fⁿ - Fⁿ⁻¹
        D = 𝒰ⁿ - 𝒰ⁿ⁻¹
    end
    
    return (2 * δc★ * A - δc² * D) / volume(i, j, k, grid, loc...)
end 

@inline function enstrophy_dissipation(i, j, k, grid, χ, fⁿ, fⁿ⁻¹, ℑxy, Uⁿ, Uⁿ⁻¹, δc★, δc²)

    C₁  = convert(eltype(grid), 1.5 + χ)
    C₂  = convert(eltype(grid), 0.5 + χ)
    loc = instantiated_location(Uⁿ)

    @inbounds begin
        𝒰ⁿ   = C₁ * ℑxy(i, j, k, grid, Uⁿ)
        𝒰ⁿ⁻¹ = C₂ * ℑxy(i, j, k, grid, Uⁿ⁻¹)
        Fⁿ   = C₁ * fⁿ[i, j, k] 
        Fⁿ⁻¹ = C₂ * fⁿ⁻¹[i, j, k] 
        A = Fⁿ - Fⁿ⁻¹
        D = 𝒰ⁿ - 𝒰ⁿ⁻¹
    end
    
    return (2 * δc★ * A - δc² * D) / volume(i, j, k, grid, loc...)
end 
