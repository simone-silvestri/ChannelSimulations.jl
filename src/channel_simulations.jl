import Oceananigans.OutputReaders: extract_field_time_series

@inline extract_field_time_series(a1, a2...) = ()
@inline extract_field_time_series(a1...) = ()

const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]

@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(j, grid, Center())
    Q = ifelse(y > p.y_shutoff, zero(grid), p.Qᵇ * cos(3π * y / p.Ly))
    return Q
end

@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))

@inline mask(y, p) = max(0, (y - p.Ly + p.Lsponge) / p.Lsponge)

@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    z = znode(k, grid, Center())
    y = ynode(j, grid, Center())

    target_b = initial_buoyancy(z, p)
    
    b = @inbounds model_fields.b[i, grid.Ny, k]

    return mask(y, p) / timescale * (target_b - b)
end

@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(j, grid, Center())
    return - p.τ * sin(π * y / p.Ly)
end

@inline u_drag(i, j, k, grid, clock, model_fields, p) = @inbounds ifelse(k == 1, - p.μ * model_fields.u[i, j, 1] / Δzᶜᶜᶜ(i, j, 1, grid), zero(grid))
@inline v_drag(i, j, k, grid, clock, model_fields, p) = @inbounds ifelse(k == 1, - p.μ * model_fields.v[i, j, 1] / Δzᶜᶜᶜ(i, j, 1, grid), zero(grid))

default_closure = ConvectiveAdjustmentVerticalDiffusivity(background_κz = 1e-5,
						          convective_κz = 0.1,
					                  background_νz = 3e-4,
						          convective_νz = 0.1)

default_momentum_advection = VectorInvariant(vertical_scheme   = WENO(),
                                             vorticity_scheme  = WENO(; order = 9),
                                             divergence_scheme = WENO())

@info "Building a model..."

default_tracer_advection = TracerAdvection(WENO(; order = 7), WENO(; order = 7), Centered()) 

function run_channel_simulation!(; momentum_advection = default_momentum_advection, 
                                     tracer_advection = default_tracer_advection, 
                                              closure = default_closure,
                                                zstar = true,
 					 initial_file = "tIni_80y_90L.bin",
                                             testcase = "0")
    # Architecture
    arch = GPU()

    # number of grid points
    Nx = 200
    Ny = 400
    Nz = 90

    # Ninty levels spacing
    Δz = [10.0 * ones(6)...,
        11.25, 12.625, 14.125, 15.8125, 17.75, 19.9375, 22.375, 25.125, 28.125, 31.625, 35.5, 39.75,
        42.0 * ones(56)...,
        39.75, 35.5, 31.625, 28.125, 25.125, 22.375, 19.9375, 17.75, 15.8125, 14.125, 12.625, 11.25,
        10.0 * ones(4)...]

    z_faces = zeros(Nz+1)
    for k in Nz : -1 : 1
        z_faces[k] = z_faces[k+1] - Δz[Nz - k + 1]
    end

    grid = RectilinearGrid(arch,
                        topology = (Periodic, Bounded, Bounded),
                        size = (Nx, Ny, Nz),
                        halo = (6, 6, 6),
                        x = (0, Lx),
                        y = (0, Ly),
                        z = z_faces)

    @info "Built a grid: $grid."

    #####
    ##### Boundary conditions
    #####

    α  = 2e-4     # [K⁻¹] thermal expansion coefficient 
    g  = 9.8061   # [m/s²] gravitational constant
    cᵖ = 3994.0   # [J/K]  heat capacity
    ρ  = 999.8    # [kg/m³] reference density

    parameters = (
        Ly  = grid.Ly,
        Lz  = grid.Lz,
        Δy  = grid.Δyᵃᶜᵃ,
        Qᵇ  = 10 / (ρ * cᵖ) * α * g, # buoyancy flux magnitude [m² s⁻³]    
        y_shutoff = 5 / 6 * grid.Ly, # shutoff location for buoyancy flux [m] 
        τ  = 0.1 / ρ,                # surface kinematic wind stress [m² s⁻²]
        μ  = 1.102e-3,               # bottom drag damping time-scale [ms⁻¹]
        Lsponge = 9 / 10 * Ly,       # sponge region for buoyancy restoring [m]
        ν  = 3e-4,                   # viscosity for "no-slip" lateral boundary conditions
        ΔB = 8 * α * g,              # surface vertical buoyancy gradient [s⁻²]
        H  =  grid.Lz,               # domain depth [m]
        h  = 1000.0,                 # exponential decay scale of stable stratification [m]
        λt = 7.0days                 # relaxation time scale [s]
    )

    # Initial condition from MITgcm
    Tinit = Array{Float64}(undef, Nx*Ny*Nz)
    read!(initial_file, Tinit)
    Tinit = bswap.(Tinit) |> Array{Float64}
    Tinit = reshape(Tinit, Nx, Ny, Nz)
    binit = reverse(Tinit, dims = 3) .* α .* g

    buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, discrete_form = true, parameters = parameters)

    buoyancy_restoring = Forcing(buoyancy_relaxation; discrete_form = true, parameters)

    u_stress_bc = FluxBoundaryCondition(u_stress; discrete_form = true, parameters)

    u_drag_forcing = Forcing(u_drag; discrete_form = true, parameters)
    v_drag_forcing = Forcing(v_drag; discrete_form = true, parameters)

    b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)
    u_bcs = FieldBoundaryConditions(top = u_stress_bc, bottom = ValueBoundaryCondition(0))
    v_bcs = FieldBoundaryConditions(bottom = ValueBoundaryCondition(0))

    #####
    ##### Coriolis
    #####

    coriolis = BetaPlane(f₀ = -1e-4, β = 1e-11)

    free_surface = SplitExplicitFreeSurface(grid; substeps = 90)

    bⁿ⁻¹ = CenterField(grid)
    Uⁿ⁻¹ = VelocityFields(grid)
    χ    = VelocityFields(grid)
    ∂b²  = VelocityFields(grid)
    ℱ    = VelocityFields(grid)

    auxiliary_fields = (; bⁿ⁻¹, 
                        uⁿ⁻¹  = Uⁿ⁻¹.u,
                        vⁿ⁻¹  = Uⁿ⁻¹.v,
                        wⁿ⁻¹  = Uⁿ⁻¹.w,
                        fˣⁿ⁻¹ = ℱ.u,
                        fʸⁿ⁻¹ = ℱ.v,
                        fᶻⁿ⁻¹ = ℱ.w,
                        χu    = χ.u,
                        χv    = χ.v,
                        χw    = χ.w,
                        ∂xb²  = ∂b².u,
                        ∂yb²  = ∂b².v,
                        ∂zb²  = ∂b².w)

    generalized_vertical_coordinate = zstar ?  ZStar() : nothing

    model = HydrostaticFreeSurfaceModel(; grid,
                                        free_surface,
                                        momentum_advection,
                                        tracer_advection,
                                        buoyancy = BuoyancyTracer(),
                                        coriolis,
                                        generalized_vertical_coordinate,
                                        closure,
                                        tracers = :b,
                                        forcing = (; b = buoyancy_restoring, u = u_drag_forcing, v = v_drag_forcing),
                                        auxiliary_fields,
                                        boundary_conditions = (b = b_bcs, u = u_bcs, v = v_bcs))

    @info "Built $model."

    model.timestepper.χ = 0.0

    #####
    ##### Initial conditions
    #####

    # resting initial condition
    bᵢ(x, y, z) = parameters.ΔB * (exp(z / parameters.h) - exp(-grid.Lz / parameters.h)) / 
                                (1 - exp(-grid.Lz / parameters.h)) * (1 + cos(20π * x / Lx) / 100)

    set!(model, b = binit) 

    #####
    ##### Simulation building
    #####

    Δt₀ = 1minutes

    # 50 years of simulation
    simulation = Simulation(model; Δt = Δt₀, stop_time = 200days)
    conjure_time_step_wizard!(simulation; cfl = 0.2, max_Δt = 5minutes, max_change = 1.1)

    # add progress callback
    wall_clock = [time_ns()]

    function print_progress(sim)
        @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u, w): (%6.3e, %6.3e) m/s, next Δt: %s\n",
            100 * (sim.model.clock.time / sim.stop_time),
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(1e-9 * (time_ns() - wall_clock[1])),
            maximum(abs, interior(sim.model.velocities.u)),
            maximum(abs, interior(sim.model.velocities.w)),
            prettytime(sim.Δt))

        wall_clock[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(20))

    run!(simulation)

    # Remove wizard now
    delete!(simulation.callbacks, :time_step_wizard)

    model.clock.time = 0
    model.clock.iteration = 0

    simulation.stop_time = 18000days
    simulation.Δt = 6minutes

    simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                            schedule = TimeInterval(18000days),
                                                            prefix = "abernathey_channel_" * string(testcase),
                                                            overwrite_existing = true)

    #####
    ##### Diagnostics
    #####

    simulation.callbacks[:compute_diagnostics] = Callback(compute_χ_values,  IterationInterval(1))
    simulation.callbacks[:update_velocities]   = Callback(update_velocities, IterationInterval(1))

    u, v, w = model.velocities
    b = model.tracers.b
    outputs = (; u, v, w, b)

    grid_variables   = zstar ? (; sⁿ = model.grid.Δzᵃᵃᶠ.sⁿ, ∂t_∂s = model.grid.Δzᵃᵃᶠ.∂t_∂s) : NamedTuple()
    snapshot_outputs = merge(model.velocities,  model.tracers)
    snapshot_outputs = merge(snapshot_outputs,  grid_variables, model.auxiliary_fields)
    average_outputs  = merge(snapshot_outputs,  model.auxiliary_fields)

    #####
    ##### Build checkpointer and output writer
    #####

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, snapshot_outputs, 
                                                            schedule = ConsecutiveIterations(TimeInterval(360days)),
                                                            filename = "abernathey_channel_snapshots_" * string(testcase),
                                                            overwrite_existing = true)

    simulation.output_writers[:averages] = JLD2OutputWriter(model, average_outputs, 
                                                            schedule = AveragedTimeInterval(5 * 360days),
                                                            filename = "abernathey_channel_averages_" * string(testcase),
                                                            overwrite_existing = true)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end
