using Statistics
using JLD2
using Printf
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.ImmersedBoundaries
using SeawaterPolynomials.TEOS10

const r_faces = [-6000.0,             -5561.590629795381, 
                 -5154.859192242641,  -4777.516750368039, 
                 -4427.439757803566,  -4102.658108238392, 
                 -3801.344048375254,  -3521.801891998017, 
                 -3262.458477264797,  -3021.854313523849, 
                 -2798.635367829655,  -2591.545444936741, 
                 -2399.4191178885562, -2221.175169417374, 
                 -2055.810507245753,  -1902.394519047109, 
                 -1760.06383529716,   -1628.0174705434547, 
                 -1505.5123157498406, -1391.8589563483943, 
                 -1286.417792464344,  -1188.5954394800156, 
                 -1097.841388681472,  -1013.6449091951929, 
                 -935.5321737800015,  -863.0635922992519, 
                 -795.8313378670168,  -733.4570517463264, 
                 -675.5897140834629,  -621.9036684955659, 
                 -572.0967893946566,  -525.8887817344348, 
                 -483.01960361144364, -443.2480028435804, 
                 -406.3501592903355,  -372.11842527424255, 
                 -340.3601570150808,  -310.89663050056856, 
                 -283.56203569246134, -258.20254340780957, 
                 -234.67543962412464, -212.84832233663514, 
                 -192.5983564478374,  -173.8115824961299, 
                 -156.3822753333063,  -140.21234914177717, 
                 -125.21080544317253, -111.29322099191894, 
                 -98.38127267184086,  -86.40229672207751, 
                 -75.28887981179821,  -64.97847966243299, 
                 -55.41307308241844,  -46.538829433725674, 
                 -38.30580769255844,  -30.667675399389157, 
                 -23.581447916685597, -17.007246526964995, 
                 -10.908074009839579, -5.249606435082352, 
                 0.0]

@inline ϕ²(i, j, k, grid, 𝒰)       = @inbounds 𝒰[i, j, k]^2
@inline speedᶠᶜᶜ(i, j, k, grid, f) = @inbounds sqrt(f.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², f.v))
@inline speedᶜᶠᶜ(i, j, k, grid, f) = @inbounds sqrt(f.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², f.u))

function run_near_global_simulation(; momentum_advection = default_momentum_advection, 
                                        tracer_advection = default_tracer_advection, 
                                                 closure = default_closure,
                                                   zstar = true,
                                            restart_file = nothing,
                                              bathymetry = "bathymetry1440x600.jld2",
                                               init_file = "initial_conditions.jld2",
                                                    arch = CPU(),
                                             timestepper = :SplitRungeKutta3,
                                                testcase = "0")
        
    # 0.25 degree resolution
    Nx = 1440
    Ny = 600
    Nz = 60

    z_faces = zstar ? MutableVerticalDiscretization(r_faces) : r_faces

    # Remember the convention!! On the surface a negative flux increases a positive decreases
    bathymetry = jldopen(bathymetry)["bathymetry"]
    bathymetry = on_architecture(arch, bathymetry)

    # A spherical domain
    grid = LatitudeLongitudeGrid(arch,
                                size = (Nx, Ny, Nz),
                                longitude = (-180, 180),
                                latitude = (-75, 75),
                                halo = (7, 7, 7),
                                z = z_faces)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry); active_cells_map = true)

    # Quadratic Drag force:
    @inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * speedᶠᶜᶜ(i, j, 1, grid, fields) * fields.u[i, j, 1]
    @inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * speedᶜᶠᶜ(i, j, 1, grid, fields) * fields.v[i, j, 1]

    @inline u_immersed_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * speedᶠᶜᶜ(i, j, 1, grid, fields) * fields.u[i, j, 1]
    @inline v_immersed_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * speedᶜᶠᶜ(i, j, 1, grid, fields) * fields.v[i, j, 1]

    # Quadratic bottom drag:
    μ = 0.003 # ms⁻¹

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form=true, parameters=μ)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(u_immersed_drag, discrete_form=true, parameters=μ))
    v_immersed_drag_bc = ImmersedBoundaryCondition(bottom = FluxBoundaryCondition(v_immersed_drag, discrete_form=true, parameters=μ))
    
    u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_drag_bc)
    v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_drag_bc)

    free_surface = SplitExplicitFreeSurface(grid; substeps=80)

    buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())

    model = HydrostaticFreeSurfaceModel(; grid,
                                          free_surface,
                                          momentum_advection,
                                          coriolis = HydrostaticSphericalCoriolis(),
                                          buoyancy,
                                          tracer_advection,
                                          tracers = (:T, :S),
                                          closure,
                                          timestepper,
                                          boundary_conditions = (u=u_bcs, v=v_bcs))

    #####
    ##### Initial condition:
    #####

    u, v, w = model.velocities
    T = model.tracers.T
    S = model.tracers.S

    if  restart_file isa String # Initialize from spinned up solution
        set!(model, restart_file)

    else # resting initial condition
        T_init = jldopen(init_file)["T"]
        S_init = jldopen(init_file)["S"]
    
        set!(model, T=T_init, S=S_init)

    end

    @info "model initialized"

    Δt₀ = 1minutes

    # First 90 days of adjustment with a lower timestep
    simulation = Simulation(model; Δt = Δt₀, stop_time = 90days, align_time_step = false)

    wall_time = Ref(time_ns())

    function progress(sim)
        u = sim.model.velocities.u
        w = sim.model.velocities.w
        v = sim.model.velocities.v

        @info @sprintf("Time: % 12s, iteration: %d, max(|u|, |w|): %.2e, %.2e ms⁻¹, wall time: %s",
                        prettytime(sim.model.clock.time),
                        sim.model.clock.iteration, maximum(abs, u), maximum(abs, w),
                        prettytime((time_ns() - wall_time[]) * 1e-9))

        wall_time[] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    #####
    ##### Simulation setup
    #####

    # Fuck the spin up!
    if !(restart_file isa String) # Spin up!    
        
        if timestepper == :SplitRungeKutta3 # Add wizard
            conjure_time_step_wizard!(simulation; cfl = 0.2, max_Δt = 5minutes, max_change = 1.1)    
        end

        simulation.output_writers[:first_checkpointer] = Checkpointer(model,
                                                                      schedule = TimeInterval(90days),
                                                                      prefix = "restart_nearglobal" * string(testcase),
                                                                      overwrite_existing = true)
        run!(simulation)

        if timestepper == :SplitRungeKutta3 # Remove wizard
            delete!(simulation.callbacks, :time_step_wizard)
        end

        # Remove first checkpoint
        delete!(simulation.output_writers, :first_checkpointer)

        # Reset time step and simulation time
        model.clock.time = 0
        model.clock.iteration = 0
    end
    
    simulation.stop_time = 1800days

    if timestepper == :SplitRungeKutta3
       simulation.Δt = 18minutes
    else
       simulation.Δt = 6minutes
    end
    
    #####
    ##### Diagnostics
    #####
    
    ϵT = Oceananigans.Simulations.VarianceDissipation(:T, grid)
    ϵS = Oceananigans.Simulations.VarianceDissipation(:S, grid)
    simulation.callbacks[:compute_T_variance] = Callback(ϵT, IterationInterval(1))
    simulation.callbacks[:compute_S_variance] = Callback(ϵS, IterationInterval(1))

    @info "added the tracer variance diagnostic"

    fT = Oceananigans.Simulations.VarianceDissipationComputations.flatten_dissipation_fields(ϵT)
    fS = Oceananigans.Simulations.VarianceDissipationComputations.flatten_dissipation_fields(ϵS)
    
    GTx = ∂x(T)^2
    GTy = ∂y(T)^2
    GTz = ∂z(T)^2

    GSx = ∂x(S)^2
    GSy = ∂y(S)^2
    GSz = ∂z(S)^2

    g = (; GTx, GTy, GTz, GSx, GSy, GSz)

    snapshot_outputs = merge(model.velocities, model.tracers, fT, fS, g)
    average_outputs  = merge(snapshot_outputs, fT, fS, g)

    #####
    ##### Build checkpointer and output writer
    #####

    if isnothing(restart_file)
        overwrite_existing = true
    else
        overwrite_existing = false
    end

    simulation.output_writers[:snapshots] = JLD2Writer(model, snapshot_outputs; 
                                                       schedule = ConsecutiveIterations(TimeInterval(100days)),
                                                       filename = "snapshots_nearglobal_" * string(testcase),
                                                       overwrite_existing)

    simulation.output_writers[:averages] = JLD2Writer(model, average_outputs; 
                                                      schedule = AveragedTimeInterval(360days),
                                                      filename = "averages_nearglobal_" * string(testcase),
                                                      overwrite_existing)

    @info "Running the simulation..."

    run!(simulation)

    return simulation
end
