using ChannelSimulations
using Oceananigans.Simulations.VarianceDissipationComputations: VarianceDissipation
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
    ϵ = VarianceDissipation(model)

    @testset "VarianceDissipation construction"  begin
        @test :b ∈ propertynames(ϵ.diffusive_production)
        @test :b ∈ propertynames(ϵ.advective_production)
    end

    simulation = Simulation(model; Δt = 1, stop_iteration = 2)

    @testset "VarianceDissipation runs"  begin
        @test begin
            ϵ(model)
            true
        end
    end

    simulation.callbacks[:dissipation] = Callback(ε, IterationInterval(1))

    run!(simulation)

    # @testset "VarianceDissipation correctness"  begin
    #     Gbx = ϵ.gradient_squared[:b].x
    #     Gby = ϵ.gradient_squared[:b].y
    #     Gbz = ϵ.gradient_squared[:b].z

    #     Khbx = ϵ.diffusive_production[:b][1].x
    #     Khby = ϵ.diffusive_production[:b][1].y
    #     Khbz = ϵ.diffusive_production[:b][1].z

    #     Kzbx = ϵ.diffusive_production[:b][2].x
    #     Kzby = ϵ.diffusive_production[:b][2].y
    #     Kzbz = ϵ.diffusive_production[:b][2].z

    #     κhx = mean(filter(!isnan, interior(Khbx) ./ interior(Gbx) ./ 2))
    #     κhy = mean(filter(!isnan, interior(Khby) ./ interior(Gby) ./ 2))
    #     κhz = mean(filter(!isnan, interior(Khbz) ./ interior(Gbz) ./ 2))
    #     κzx = mean(filter(!isnan, interior(Kzbx) ./ interior(Gbx) ./ 2))
    #     κzy = mean(filter(!isnan, interior(Kzby) ./ interior(Gby) ./ 2))
    #     κzz = mean(filter(!isnan, interior(Kzbz) ./ interior(Gbz) ./ 2))

    #     @test κhx ≈ 1.0 
    #     @test κhy ≈ 1.0
    #     @test κhz ≈ 0.0

    #     @test κzx ≈ 0.0
    #     @test κzy ≈ 0.0
    #     @test κzz ≈ 1e-4

    #     Abx = ϵ.advective_production[:b].x
    #     Aby = ϵ.advective_production[:b].y
    #     Abz = ϵ.advective_production[:b].z

    #     κAx = mean(filter(!isnan, interior(Abx) ./ interior(Gbx) ./ 2))
    #     κAy = mean(filter(!isnan, interior(Aby) ./ interior(Gby) ./ 2))
    #     κAz = mean(filter(!isnan, interior(Abz) ./ interior(Gbz) ./ 2))

    #     @test κhx ≈ 1.0 
    #     @test κhy ≈ 1.0
    #     @test κhz ≈ 0.0
    # end

end
