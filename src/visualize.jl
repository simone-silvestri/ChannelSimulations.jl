using CairoMakie, SixelTerm, JLD2
using Statistics: mean

function visualize_field(; file = jldopen("abernathey_channel_averages.jld2"),
                    tracer_name = :b, 
                      iteration = keys(file["timeseries/t"])[end],
                      with_halo = 6)
    
    Abx = file["timeseries/A" * string(tracer_name) *  "x/" * iteration] 
    Aby = file["timeseries/A" * string(tracer_name) *  "y/" * iteration] 
    Abz = file["timeseries/A" * string(tracer_name) *  "z/" * iteration]  
    
    Dbx = file["timeseries/D" * string(tracer_name) *  "x/" * iteration] 
    Dby = file["timeseries/D" * string(tracer_name) *  "y/" * iteration] 
    Dbz = file["timeseries/D" * string(tracer_name) *  "z/" * iteration]  
    
    Gbx = file["timeseries/G" * string(tracer_name) *  "x/" * iteration]
    Gby = file["timeseries/G" * string(tracer_name) *  "y/" * iteration]
    Gbz = file["timeseries/G" * string(tracer_name) *  "z/" * iteration]

    if !isnothing(with_halo)
        H = with_halo

        Abx = Abx[H+1:end-H, H+1:end-H, H+1:end-H]     
        Aby = Aby[H+1:end-H, H+1:end-H, H+1:end-H]  
        Abz = Abz[H+1:end-H, H+1:end-H, H+1:end-H]  
        Dbx = Dbx[H+1:end-H, H+1:end-H, H+1:end-H]  
        Dby = Dby[H+1:end-H, H+1:end-H, H+1:end-H]  
        Dbz = Dbz[H+1:end-H, H+1:end-H, H+1:end-H]  
        Gbx = Gbx[H+1:end-H, H+1:end-H, H+1:end-H]  
        Gby = Gby[H+1:end-H, H+1:end-H, H+1:end-H]  
        Gbz = Gbz[H+1:end-H, H+1:end-H, H+1:end-H]  
    end

    cxm = mean(Abx, dims = 1)[1, :, :]
    cym = mean(Aby, dims = 1)[1, :, :]
    czm = mean(Abz, dims = 1)[1, :, :]

    dxm = mean(Dbx, dims = 1)[1, :, :]
    dym = mean(Dby, dims = 1)[1, :, :]
    dzm = mean(Dbz, dims = 1)[1, :, :]

    bxm = mean(Gbx, dims = 1)[1, :, :]
    bym = mean(Gby, dims = 1)[1, :, :]
    bzm = mean(Gbz, dims = 1)[1, :, :]

    κx = - cxm ./ bxm ./ 2
    κy = - cym ./ bym ./ 2
    κz = - czm ./ bzm ./ 2
    κdz = - dzm ./ bzm ./ 2

    cxc = deepcopy(cxm)
    cyc = (cym[1:end-1, :] .+ cym[2:end, :]) / 2
    czc = (czm[:, 1:end-1] .+ czm[:, 2:end]) / 2

    bxc = deepcopy(bxm)
    byc = (bym[1:end-1, :] .+ bym[2:end, :]) / 2
    bzc = (bzm[:, 1:end-1] .+ bzm[:, 2:end]) / 2

    κi = - (cxc .+ cyc .+ czc) ./ (bxc .+ byc .+ bzc) ./ 2

    κix = - cxc ./ bzc ./ 2
    κiy = - cyc ./ bzc ./ 2 
    κiz = - czc ./ bzc ./ 2 

    κx[isnan.(κx)] .= 0
    κy[isnan.(κy)] .= 0
    κz[isnan.(κz)] .= 0
    κi[isnan.(κi)] .= 0

    Nz = 90
    Δz = [10.0 * ones(6)...,
        11.25, 12.625, 14.125, 15.8125, 17.75, 19.9375, 22.375, 25.125, 28.125, 31.625, 35.5, 39.75,
        42.0 * ones(56)...,
        39.75, 35.5, 31.625, 28.125, 25.125, 22.375, 19.9375, 17.75, 15.8125, 14.125, 12.625, 11.25,
        10.0 * ones(4)...]

    z_faces = zeros(Nz+1)
    for k in Nz : -1 : 1
        z_faces[k] = z_faces[k+1] - Δz[Nz - k + 1]
    end

    zF = z_faces
    zC = (zF[1:end-1] .+ zF[2:end]) ./ 2

    Δh = 5000

    irange = 1:390

    κixm = mean(κix[irange, :], dims = 1)[1, :] 
    κiym = mean(κiy[irange, :], dims = 1)[1, :] 
    κizm = mean(κiz[irange, :], dims = 1)[1, :] 
    κxm = mean(κx[irange, :], dims = 1)[1, :] 
    κym = mean(κy[irange, :], dims = 1)[1, :] 
    κzm = mean(κz[irange, :], dims = 1)[1, :] 
    κdm = mean(κdz[irange, :], dims = 1)[1, :] 
    κim = mean(κi[irange, :], dims = 1)[1, :] 

    fig = Figure(); ax = Axis(fig[1, 1], ylabel = "z [m]", xlabel = "Diffusivity [m²/s]")

    lines!(ax, κixm, zC, color = :blue,  linewidth = 2, label = "κx Sx² = - ⟨Px⟩ / ⟨(∂zb)²⟩") 
    lines!(ax, κiym, zC, color = :green, linewidth = 2, label = "κy Sy² = - ⟨Py⟩ / ⟨(∂zb)²⟩") 
    lines!(ax, κzm,  zF, color = :red,   linewidth = 2, label = "κz     = - ⟨Pz⟩ / ⟨(∂zb)²⟩") 

    vlines!(ax, -1e-5; linestyle = :dash, linewidth = 0.5, color = :grey)
    vlines!(ax,  1e-5; linestyle = :dash, linewidth = 0.5, color = :grey)
    axislegend(ax, position = :rc) 
  #  xlims!(ax, (-2e-5, 3.5e-5))

    fig2 = Figure(); ax = Axis(fig2[1, 1], ylabel = "z [m]", xlabel = "log(abs(diffusivity))")

    lines!(ax, log10.(abs.(κixm)), zC, color = :blue,  linewidth = 2, label = "abs(κx Sx²)") 
    lines!(ax, log10.(abs.(κiym)), zC, color = :green, linewidth = 2, label = "abs(κy Sy²)") 
    lines!(ax, log10.(abs.(κzm)),  zF, color = :red,   linewidth = 2, label = "abs(κz)") 
    lines!(ax, log10.(abs.(κixm .+ κiym .+ κizm)), zC, color = :black, linewidth = 2, label = "abs(total)", linestyle = :dash) 

    vlines!(ax, -5; linestyle = :dash, linewidth = 0.5, color = :grey)
    axislegend(ax, position = :lc) 

#     xlims!(ax, (-10, -3))

    return fig, fig2, (; κx, κy, κz, κixm, κiym, κdm, κzm, κxm, κym, cxm, cym, czm, bxm, bym, bzm, zC, zF, dxm, dym, dzm)
end

function plot_cases(nt1, nt2, nt3, nt4, var = :κixm; var2 = nothing)

    fig = Figure(); ax = Axis(fig[1, 1], ylabel = "z [m]", xlabel = "Diffusivity [m²/s]")

    t1 = getproperty(nt1, var)
    t2 = getproperty(nt2, var)
    t3 = getproperty(nt3, var)
    t4 = getproperty(nt4, var)

    if length(size(t1)) == 2
        t1 = mean(t1, dims=1)[1, :]
        t2 = mean(t2, dims=1)[1, :]
        t3 = mean(t3, dims=1)[1, :]
        t4 = mean(t4, dims=1)[1, :]
    end
   
    if !isnothing(var2)
        t21 = getproperty(nt1, var)
        t22 = getproperty(nt2, var)
        t23 = getproperty(nt3, var)
        t24 = getproperty(nt4, var)

        if length(size(t21)) == 2
            t21 = mean(t21, dims=1)[1, :]
            t22 = mean(t22, dims=1)[1, :]
            t23 = mean(t23, dims=1)[1, :]
            t24 = mean(t24, dims=1)[1, :]
        end
   
        t1 .+= t21
        t2 .+= t22
        t3 .+= t23
        t4 .+= t24
    end

    z = if length(t1) == 91
        :zF
    else
        :zC
    end

    lines!(ax, t1, getproperty(nt1, z), color = :blue,   linewidth = 2, label = "nt1") 
    lines!(ax, t2, getproperty(nt2, z), color = :green,  linewidth = 2, label = "nt2") 
    lines!(ax, t3, getproperty(nt3, z), color = :red,    linewidth = 2, label = "nt3") 
    lines!(ax, t4, getproperty(nt4, z), color = :purple, linewidth = 2, label = "nt4") 

    vlines!(ax, -1e-5; linestyle = :dash, linewidth = 0.5, color = :grey)
    vlines!(ax,  1e-5; linestyle = :dash, linewidth = 0.5, color = :grey)
    axislegend(ax, position = :rc) 

    fig2 = Figure(); ax = Axis(fig2[1, 1], ylabel = "z [m]", xlabel = "log(abs(diffusivity))")
    
    lines!(ax, log10.(abs.(t1)), getproperty(nt1, z), color = :blue,   linewidth = 2, label = "nt1") 
    lines!(ax, log10.(abs.(t2)), getproperty(nt2, z), color = :green,  linewidth = 2, label = "nt2") 
    lines!(ax, log10.(abs.(t3)), getproperty(nt3, z), color = :red,    linewidth = 2, label = "nt3") 
    lines!(ax, log10.(abs.(t4)), getproperty(nt4, z), color = :purple, linewidth = 2, label = "nt4") 

    vlines!(ax, -5; linestyle = :dash, linewidth = 0.5, color = :grey)
    axislegend(ax, position = :lc) 

    return fig, fig2
end
