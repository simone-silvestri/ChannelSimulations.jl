using ChannelSimulations.VarianceDissipations

const Lx = 1000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]

@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(i, j, 1, grid, Center(), Center(), Center())
    Q = ifelse(y > p.y_shutoff, zero(grid), p.Qᵇ * cos(3π * y / p.Ly))
    return Q
end

@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))

@inline mask(y, p) = max(0, (y - p.Ly + p.Lsponge) / p.Lsponge)

@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    z = znode(i, j, k, grid, Center(), Center(), Center())
    y = ynode(i, j, k, grid, Center(), Center(), Center())

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

default_bottom_height = (x, y) -> y < 1000kilometers ?  5.600000000000001e-15 * y^3 - 8.4e-9 * y^2 - 200 : -3000.0

function default_grid(arch, zstar, bottom_height)

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

    z_faces = zstar ? MutableVerticalDiscretization(z_faces) : z_faces

    grid = RectilinearGrid(arch,
                        topology = (Periodic, Bounded, Bounded),
                        size = (Nx, Ny, Nz),
                        halo = (6, 6, 6),
                        x = (0, Lx),
                        y = (0, Ly),
                        z = z_faces)

    @info "Built a grid: $grid."

    return isnothing(bottom_height) ? grid : ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
end

hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

function run_channel_simulation(; momentum_advection = default_momentum_advection, 
                                    tracer_advection = default_tracer_advection, 
                                             closure = default_closure,
                                               zstar = false,
                                        restart_file = nothing,
                                                arch = CPU(),
                                       bottom_height = nothing,
                                         timestepper = :SplitRungeKutta3,
                                                grid = default_grid(arch, zstar, bottom_height),
 					                    initial_file = "tIni_80y_90L.bin",
                                            testcase = "0")

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
        μ  = 1.1e-3,               # bottom drag damping time-scale [ms⁻¹]
        Lsponge = 9 / 10 * Ly,       # sponge region for buoyancy restoring [m]
        ν  = 3e-4,                   # viscosity for "no-slip" lateral boundary conditions
        ΔB = 8 * α * g,              # surface vertical buoyancy gradient [s⁻²]
        H  =  grid.Lz,               # domain depth [m]
        h  = 1000.0,                 # exponential decay scale of stable stratification [m]
        λt = 7.0days                 # relaxation time scale [s]
    )

    buoyancy_flux_bc   = FluxBoundaryCondition(buoyancy_flux, discrete_form = true, parameters = parameters)
    buoyancy_restoring = Forcing(buoyancy_relaxation; discrete_form = true, parameters)
    u_stress_bc        = FluxBoundaryCondition(u_stress; discrete_form = true, parameters)

    # Drag is added as a forcing to allow both bottom drag _and_ a no-slip BC
    u_drag_forcing = Forcing(u_drag; discrete_form = true, parameters)
    v_drag_forcing = Forcing(v_drag; discrete_form = true, parameters)

    b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)
    u_bcs = FieldBoundaryConditions(bottom = ValueBoundaryCondition(0), top = u_stress_bc)
    v_bcs = FieldBoundaryConditions(bottom = ValueBoundaryCondition(0))

    #####
    ##### Coriolis
    #####

    coriolis = BetaPlane(f₀ = -1e-4, β = 1e-11)
    free_surface = SplitExplicitFreeSurface(grid; substeps = 90)

    vertical_coordinate = zstar ? ZStar() : ZCoordinate()

    tracers = hasclosure(closure, CATKEVerticalDiffusivity) ? (:b, :e) : (:b, )

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection,
                                          tracer_advection,
                                          buoyancy = BuoyancyTracer(),
                                          coriolis,
                                          closure,
                                          tracers,
                                          timestepper,
                                          vertical_coordinate,
                                          forcing = (; b = buoyancy_restoring, u = u_drag_forcing, v = v_drag_forcing),
                                          boundary_conditions = (b = b_bcs, u = u_bcs, v = v_bcs))

    @info "Built $model."

    #####
    ##### Initial conditions
    #####

    Nx, Ny, Nz = size(grid)

    if  restart_file isa String # Initialize from spinned up solution
        set!(model, restart_file)

    else # resting initial condition
      
      binit = if initial_file isa String
          
          # Initial condition from MITgcm
          Tinit = Array{Float64}(undef, Nx*Ny*Nz)
          read!(initial_file, Tinit)
          Tinit = bswap.(Tinit) |> Array{Float64}
          Tinit = reshape(Tinit, Nx, Ny, Nz)
          binit = reverse(Tinit, dims = 3) .* α .* g

      else
          (x, y, z) -> initial_buoyancy(z, parameters)
      end
          

      set!(model, b = binit) 

    end

    #####
    ##### Simulation building
    #####

    Δt₀ = 2minutes

    # 50 years of simulation
    simulation = Simulation(model; Δt = Δt₀, stop_time = 150days)

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

    # Fuck the spin up!
    if !(restart_file isa String) # Spin up!        
        run!(simulation)

        # Reset time step and simulation time
        model.clock.time = 0
        model.clock.iteration = 0
    end

    simulation.stop_time = 14400days

    if timestepper == :SplitRungeKutta3
       simulation.Δt = 15minutes
    else
       simulation.Δt = 5minutes
    end

    simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                            schedule = TimeInterval(1800days),
                                                            prefix = "channel_checkpoint_" * string(testcase),
                                                            overwrite_existing = true)
    #####
    ##### Diagnostics
    #####
    
    ϵ = Oceananigans.Simulations.VarianceDissipation(:b, grid)
    simulation.callbacks[:compute_variance] = Callback(ϵ, IterationInterval(1))

    @info "added the tracer variance diagnostic"

    f = Oceananigans.Simulations.VarianceDissipationComputations.flatten_dissipation_fields(ϵ)
    b = model.tracers.b

    Gbx = ∂x(b)^2
    Gby = ∂y(b)^2
    Gbz = ∂z(b)^2

    g = (; Gbx, Gby, Gbz)

    snapshot_outputs = merge(model.velocities, model.tracers, f, g)
    average_outputs  = merge(snapshot_outputs, f, g)

    #####
    ##### Build checkpointer and output writer
    #####

    if isnothing(restart_file)
        overwrite_existing = true
    else
        overwrite_existing = false
    end

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, snapshot_outputs; 
                                                            schedule = ConsecutiveIterations(TimeInterval(360days)),
                                                            filename = "snapshots_" * string(testcase),
                                                            overwrite_existing)

    simulation.output_writers[:averages] = JLD2OutputWriter(model, average_outputs; 
                                                            schedule = AveragedTimeInterval(5 * 360days),
                                                            filename = "averages_" * string(testcase),
                                                            overwrite_existing)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end
