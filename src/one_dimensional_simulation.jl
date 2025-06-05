using Oceananigans.Models.VarianceDissipationComputations: VarianceDissipation

using Oceananigans

# the initial condition
@inline G(x, β, z) = exp(-β*(x - z)^2)
@inline F(x, α, a) = √(max(1 - α^2*(x-a)^2, 0.0))

Z = -0.7
δ = 0.005
β = log(2)/(36*δ^2)
a = 0.5
α = 10

@inline function c₀(x) 
    if x <= -0.6 && x >= -0.8
        return 1/6*(G(x, β, Z-δ) + 4*G(x, β, Z) + G(x, β, Z+δ))
    elseif x <= -0.2 && x >= -0.4
        return 1.0
    elseif x <= 0.2 && x >= 0.0
        return 1.0 - abs(10 * (x - 0.1))
    elseif x <= 0.6 && x >= 0.4
        return 1/6*(F(x, α, a-δ) + 4*F(x, α, a) + F(x, α, a+δ))
    else
        return 0.0
    end
end

function run_one_dimensional_simulation(N; advection = WENO(order=7), timestepper = :SplitRungeKutta3)
    # 1D grid constructions
    grid = RectilinearGrid(size=N, x=(-1, 1), halo=7, topology = (Periodic, Flat, Flat))


    Δt_max = 0.2 * minimum_xspacing(grid)
    c_real = CenterField(grid)
    set!(c_real, c₀)

    model = HydrostaticFreeSurfaceModel(; grid, 
                                          velocities=PrescribedVelocityFields(u=1),
                                          timestepper,
                                          tracer_advection=advection,
                                          tracers=:c)

    set!(model, c=c₀)

    variance_dissipation = VarianceDissipation(model)

    simulation = Simulation(model, Δt=Δt_max, stop_time=10)

    simulation.callbacks[:compute_variance] = Callback(variance_dissipation, IterationInterval(1), callsite=Oceananigans.TendencyCallsite())
    @info "added the tracer variance diagnostic"

    snapshot_outputs = merge(model.tracers, get_dissipation_fields(variance_dissipation))

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, snapshot_outputs; 
                                                            schedule = ConsecutiveIterations(TimeInterval(0.1)),
                                                            filename = "snapshots_one_dimensional",
                                                            overwrite_existing = true)

    run!(simulation)

    return simulation
end
