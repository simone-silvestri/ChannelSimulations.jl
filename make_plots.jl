using Statistics

include("comparison.jl")
include("src/visualize.jl")

var = read_variable("varDiag.0001036800.data", 12);
_, _, k = visualize_field(file = jldopen("abernathey_channel_averages_mitgcm_weno9_second.jld2"))

Py = k.cym
Px = k.cxm 
Pz = k.czm

By = k.bym
Bx = k.bxm
Bz = k.bzm
 
kx = Px ./ 2 ./ Bx
ky = Py ./ 2 ./ By
kz = Pz ./ 2 ./ Bz

Bzi = (Bz[:, 1:end-1] .+ Bz[:, 2:end]) ./ 2
kxi = Px ./ 2 ./ Bzi 
kyi = (Py[1:end-1, :] .+ Py[2:end, :]) ./ 4 ./ Bzi 
kzi = (kz[:, 1:end-1] .+ kz[:, 2:end]) ./ 2

Δz = abs.(k.zF[1:end-1] .- k.zF[2:end])

Pxm = mean(var[:, :, :, 5] , dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2  
Pym = mean(var[:, :, :, 6] , dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2 
Pzm = mean(var[:, :, :, 7] , dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2 
Bxm = mean(var[:, :, :, 9] , dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2 ./ 5000 
Bym = mean(var[:, :, :, 10], dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2 ./ 5000 
Bzm = mean(var[:, :, :, 11], dims = 1)[1, :, :] .* 4e-8 .* 9.8061^2 

for k in 1:90
    Bzm[:, k] ./= Δz[k]
end

kxm = Pxm ./ 2 ./ Bxm
kym = Pym ./ 2 ./ Bym
kzm = Pzm ./ 2 ./ Bzm

kxim = Pxm ./ 2 ./ Bzm
kyim = Pym ./ 2 ./ Bzm

fig = Figure(size = (1200, 600), fontsize=8); axis = Axis(fig[1, 1], ylabel = "z [m]", xscale = log10, xlabel = "[s⁻⁴]")
lines!(axis, abs.(mean(Pxm[3:end-3, :], dims = 1)[1, :]), k.zC, label = "Px, MITgcm", color = :blue)
lines!(axis, abs.(mean(Px[3:end-3, :] , dims = 1)[1, :]), k.zC, label = "Px, Ocean", color = :blue, linestyle = :dash)
lines!(axis, abs.(mean(Pym[3:end-3, :], dims = 1)[1, :]), k.zC, label = "Py, MITgcm", color = :deepskyblue1)
lines!(axis, abs.(mean(Py[3:end-3, :] , dims = 1)[1, :]), k.zC, label = "Py, Ocean", color = :deepskyblue1, linestyle = :dash)
lines!(axis, abs.(mean(Bxm[3:end-3, :], dims = 1)[1, :]), k.zC, label = "∂ˣb², MITgcm", color = :red)
lines!(axis, abs.(mean(Bx[3:end-3, :] , dims = 1)[1, :]), k.zC, label = "∂ˣb², Ocean", color = :red, linestyle = :dash)
lines!(axis, abs.(mean(Bym[3:end-3, :], dims = 1)[1, :]), k.zC, label = "∂ʸb², MITgcm", color = :orange)
lines!(axis, abs.(mean(By[3:end-3, :] , dims = 1)[1, :]), k.zC, label = "∂ʸb², Ocean", color = :orange, linestyle = :dash)
lines!(axis, abs.(mean(Bzm[3:end-3, :], dims = 1)[1, :]), k.zC, label = "∂ᶻb², MITgcm", color = :purple)
lines!(axis, abs.(mean(Bzi[3:end-3, :] , dims = 1)[1, :]), k.zC, label = "∂ᶻb², Ocean", color = :purple, linestyle = :dash)
axislegend(axis, position = :lt, framecolor = :transparent)

axis = Axis(fig[1, 2], ylabel = "z [m]", xlabel = "[m²s⁻¹]")
lines!(axis, - mean(kxm[3:end-3, :], dims = 1)[1, :], k.zC, label = "κx, MITgcm", color = :red)
lines!(axis, - mean(kx[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κx, Ocean", color = :red, linestyle = :dash) 
lines!(axis, - mean(kym[3:end-3, :], dims = 1)[1, :], k.zC, label = "κy, MITgcm",  color = :blue)
lines!(axis, - mean(ky[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κy, Ocean",  color = :blue, linestyle = :dash)
xlims!(axis, (-1, 10))
axislegend(axis, position = :rc)

kt  = kxi  .+ kyi  .+ kzi
ktm = kxim .+ kyim .+ kzm 

axis = Axis(fig[1, 3], ylabel = "z [m]", xlabel = "[m²s⁻¹]")
lines!(axis, - mean(kxim[3:end-3, :], dims = 1)[1, :], k.zC, label = "κxi, MITgcm", color = :red)
lines!(axis, - mean(kxi[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κxi, Ocean", color = :red, linestyle = :dash) 
lines!(axis, - mean(kyim[3:end-3, :], dims = 1)[1, :], k.zC, label = "κyi, MITgcm",  color = :blue)
lines!(axis, - mean(kyi[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κyi, Ocean",  color = :blue, linestyle = :dash)
lines!(axis, - mean(kzm[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κz, MITgcm",  color = :green)
lines!(axis, - mean(kzi[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κz, Ocean",  color = :green, linestyle = :dash)
lines!(axis, - mean(ktm[3:end-3, :],  dims = 1)[1, :], k.zC, label = "κi, MITgcm",  color = :black)
lines!(axis, - mean(kt[3:end-3, :],   dims = 1)[1, :], k.zC, label = "κi, Ocean",  color = :black, linestyle = :dash)
vlines!(axis,   1e-5, linewidth = 0.2, linestyle = :dash, color = :grey)
vlines!(axis, - 1e-6, linewidth = 0.2, linestyle = :dash, color = :grey)
xlims!(axis, (-1e-5, 1e-4))
axislegend(axis, position = :rc)


