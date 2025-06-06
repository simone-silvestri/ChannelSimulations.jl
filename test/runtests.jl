using ChannelSimulations
using Oceananigans.Models.VarianceDissipationComputations: VarianceDissipation
using Test
using Oceananigans

@testset "ChannelSimulations.jl" begin
    grid  = RectilinearGrid(size = (10, 10, 10), extent = (100, 100, 100))
    
    hclosure = HorizontalScalarDiffusivity(κ = 1.0)
    vclosure = VerticalScalarDiffusivity(κ = 1e-4)

    closure = (hclosure, vclosure)
    tracer_advection = UpwindBiased(order = 1)

    velocities = PrescribedVelocityFields(u = 0.1, v = 0.1, w = 0.0)

    model = HydrostaticFreeSurfaceModel(; grid, 
                                          closure, 
                                          tracers = :b, 
                                          buoyancy = nothing,
                                          velocities,
                                          tracer_advection)
    
    set!(model, b = (x, y, z) -> rand())
    ϵ = VarianceDissipation(:b, grid)

    simulation = Simulation(model; Δt = 1, stop_iteration = 2)

    @testset "VarianceDissipation runs"  begin
        @test begin
            ϵ(model)
            true
        end
    end

    simulation.callbacks[:dissipation] = Callback(ϵ, IterationInterval(1))

    @testset "Simulation runs" begin
        @test begin
            run!(simulation)
            true
        end
    end
end
