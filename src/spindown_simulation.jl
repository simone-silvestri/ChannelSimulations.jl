
function run_spindown_simulation(; momentum_advection = default_momentum_advection, 
                                     tracer_advection = default_tracer_advection, 
                                              closure = default_closure,
                                                zstar = true,
 					                     restart_file = nothing,
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
    ##### Coriolis
    #####

    coriolis = BetaPlane(f₀ = -1e-4, β = 1e-11)

    free_surface = SplitExplicitFreeSurface(grid; substeps = 90)

    generalized_vertical_coordinate = zstar ?  ZStar() : nothing

    model = HydrostaticFreeSurfaceModel(; grid,
                                        free_surface,
                                        momentum_advection,
                                        tracer_advection,
                                        buoyancy = BuoyancyTracer(),
                                        coriolis,
                                        generalized_vertical_coordinate,
                                        closure,
                                        tracers = (:b, :c))
    @info "Built $model."

    model.timestepper.χ = 0.0

    #####
    ##### Initial conditions
    #####

    isnothing(restart_file) && throw("Restart file cannot be a 'Nothing'!!!")
    set!(model, restart_file)

    # tracer initial condition gaussian centered at 
    # - 1500 meters depth 
    # - Lx / 2 in x 
    # - Ly / 2 in y
    # spread of 100 km in x and y and 100 meters in z
    x₀ = Lx / 2
    y₀ = Ly / 2
    z₀ = Lz / 2
    σˣ = σʸ = 100kilometers
    σᶻ = 100
    cᵢ(x, y, z) = exp( - (x - x₀)^2 / σˣ^2 - (y - y₀)^2 / σʸ^2 - (z - z₀)^2 / σᶻ^2)

    set!(model, c = cᵢ) 

    #####
    ##### Simulation building
    #####

    Δt = 6minutes

    # 50 years of simulation
    simulation = Simulation(model; Δt, stop_time = 360days)

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

    #####
    ##### Diagnostics
    #####

    grid_variables   = zstar ? (; sⁿ = model.grid.Δzᵃᵃᶠ.sⁿ, ∂t_∂s = model.grid.Δzᵃᵃᶠ.∂t_∂s) : NamedTuple()
    snapshot_outputs = merge(model.velocities,  model.tracers)
    snapshot_outputs = merge(snapshot_outputs,  grid_variables, model.auxiliary_fields)

    #####
    ##### Build checkpointer and output writer
    #####

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, snapshot_outputs, 
                                                            schedule = ConsecutiveIterations(TimeInterval(60days)),
                                                            filename = "spindown_snapshots_" * string(testcase),
                                                            overwrite_existing = true)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end
