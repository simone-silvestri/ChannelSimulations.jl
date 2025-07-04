### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ec938fce-60b1-11ef-3e29-f5396e753494
begin
	using Pkg
	Pkg.add("CairoMakie")
	Pkg.add("HypertextLiteral")
	Pkg.add("PlutoUI")
	using Statistics
	using PlutoUI
	using CairoMakie
	using JLD2
	using HypertextLiteral
	nothing
end

# ╔═╡ cdb48125-a977-4e3a-91ec-1861792dd9d1
@bind screenWidth @htl("""
	<div>
	<script>
		var div = currentScript.parentElement
		div.value = screen.width
	</script>
	</div>
""")

# ╔═╡ 1c578cdd-cd83-499b-98a0-1292b566a0f8
begin
	cellWidth= min(1000, screenWidth*0.90)
	@htl("""
		<style>
			pluto-notebook {
				margin: auto;
				width: $(cellWidth)px;
			}
		</style>
	""")
end


# ╔═╡ bc52c23b-5573-4931-93a5-5b598e8cc83b
# Some constants
begin
	α  = 2e-4
	g  = 9.81
	Δh = 5000

    # Ninty levels spacing
    Δz = [10.0 * ones(6)...,
          11.25, 12.625, 14.125, 15.8125, 17.75, 19.9375, 22.375, 25.125, 28.125, 31.625, 35.5, 39.75,
          42.0 * ones(56)...,
          39.75, 35.5, 31.625, 28.125, 25.125, 22.375, 19.9375, 17.75, 15.8125, 14.125, 12.625, 11.25,
          10.0 * ones(4)...]

	
    zF = zeros(91)
    for k in 90 : -1 : 1
        zF[k] = zF[k+1] - Δz[90 - k + 1]
    end
	
	zC = (zF[2:end] + zF[1:end-1]) / 2

	nothing
end

# ╔═╡ 37f329a1-c664-4cc9-acab-c6c5d357c7da
begin
	function instantaneous_heatmaps(v1, v2, v3, colormap, title; fcrange = nothing, fxlim = nothing, labpos = :lt)

		Ly = range(0, 2, length=size(v1, 2))
		Lx = range(0, 1, length=size(v1, 1))

		fig = Figure(size = (1100, 600), fontsize = 15)
		ax = Axis(fig[1, 1], xlabel = "y [km]", ylabel = "x [km]", yticks = ([0, 0.5, 1, 1.5, 2], ["0", "0.5", "1", "1.5", "2"]))
		hm = isnothing(fcrange) ? heatmap!(ax, Lx, Ly, v1; colormap) : heatmap!(ax, Lx, Ly, v1; colormap, colorrange = fcrange) 
				
		ax  = Axis(fig[1, 2], xlabel = "y [km]", ylabel = "", yticks = ([0, 0.5, 1, 15, 2], ["", "", "", "", ""]))
		hm = isnothing(fcrange) ? heatmap!(ax, Lx, Ly, v2; colormap) : heatmap!(ax, Lx, Ly, v2; colormap, colorrange = fcrange) 
	
		ax  = Axis(fig[1, 3], xlabel = "y [km]", ylabel = "", yticks = ([0, 0.5, 1, 15, 2], ["", "", "", "", ""]))
		hm = isnothing(fcrange) ? heatmap!(ax, Lx, Ly, v3; colormap) : heatmap!(ax, Lx, Ly, v3; colormap, colorrange = fcrange) 
		cb = Colorbar(fig[0, 1:3], hm, vertical = false, label = title, width = Relative(2/3))

		save(title * ".png", fig)
		
		return fig
	end
	
	# a utility function for plotting:
	function plot_and_save(v1, v2, v3, vl, title; fxlim = nothing, labpos = :lt)

		lim = 0
		
		v11 = v1[1:end-lim, :]
		v21 = v2[1:end-lim, :]
		v31 = v3[1:end-lim, :]
		
		L = range(0, 2000, length=size(v11, 1))
		
		fig = Figure(size = (400, 400), fontsize = 16)
		D1 = mean(v11, dims = 1)[1, :]
		D2 = mean(v21, dims = 1)[1, :]
		D3 = mean(v31, dims = 1)[1, :]
		
		ax = Axis(fig[1, 1]; xlabel = "Average " * title, ylabel = "Depth [km]")
		lines!(ax, D1, zC ./ 1e3, label = "Case 1")
		lines!(ax, D2, zC ./ 1e3, linestyle = :dash, label = "Case 2")
		lines!(ax, D3, zC ./ 1e3, linestyle = :dashdot, label = "Case 3")
		axislegend(; framecolor = :transparent, position = labpos, backgroundcolor = :transparent)
		if !isnothing(fxlim)
			xlims!(ax, fxlim)
		end

		vlines!(ax, vl; linestyle = :dash, color = :grey, linewidth = 0.75)
		save(title * "onlyplot.png", fig)

		return fig
	end

	# a utility function for plotting:
	function plot_heatmaps(v1, v2, v3, colorrange, colormap, title; fcrange = nothing, fxlim = nothing, labpos = :lt)

		lim = 0
		
		v11 = v1[1:end-lim, :]
		v21 = v2[1:end-lim, :]
		v31 = v3[1:end-lim, :]
		
		L = range(0, 2000, length=size(v11, 1))
		
		fig = Figure(size = (1200, 400), fontsize = 10)
		ax = Axis(fig[1, 1], xlabel = "y [km]", ylabel = "Depth [km]", yticks = ([-3, -2, -1, 0], ["-3", "-2", "-1", "0"]))
		hm = isnothing(fcrange) ? heatmap!(ax, L, zC./ 1e3, v11; colormap) : heatmap!(ax, L, zC ./ 1e3, v11; colormap, colorrange = fcrange) 
		
		ax  = Axis(fig[1, 2], xlabel = "y [km]", ylabel = "", yticks = ([-3, -2, -1, 0], ["", "", "", ""]))
		hm = isnothing(fcrange) ? heatmap!(ax, L, zC ./ 1e3, v21; colormap) : heatmap!(ax, L, zC./ 1e3, v21; colormap, colorrange = fcrange) 
		
		ax  = Axis(fig[1, 3], xlabel = "y [km]", ylabel = "", yticks = ([-3, -2, -1, 0], ["", "", "", ""]))
		hm = isnothing(fcrange) ? heatmap!(ax, L, zC ./ 1e3, v31; colormap) : heatmap!(ax, L, zC./ 1e3, v31; colormap, colorrange = fcrange) 
		cb = Colorbar(fig[0, 1:3], hm, vertical = false, label = title, width = Relative(2/3))
	
		D1 = mean(v11, dims = 1)[1, :]
		D2 = mean(v21, dims = 1)[1, :]
		D3 = mean(v31, dims = 1)[1, :]
		
		ax = Axis(fig[1, 4]; xlabel = "Meridional - mean " * title, ylabel = "Depth [km]")
		lines!(ax, D1, zC ./ 1e3, label = "Case 1")
		lines!(ax, D2, zC ./ 1e3, linestyle = :dash, label = "Case 2")
		lines!(ax, D3, zC ./ 1e3, linestyle = :dashdot, label = "Case 3")
		axislegend(; framecolor = :transparent, position = labpos, backgroundcolor = :transparent)
		if !isnothing(fxlim)
			xlims!(ax, fxlim)
		end
		
		save(title * ".png", fig)

		return fig
	end

	function plot_contours(v1, v2, v3, colorrange, colormap, title; levels = 10, fxlim = nothing, labpos = :lt)

		lim = 0
		
		v11 = v1[1:end-lim, :]
		v21 = v2[1:end-lim, :]
		v31 = v3[1:end-lim, :]
		
		L = range(0, 2000, length=size(v11, 1))
				
		fig = Figure(size = (1200, 400), fontsize = 10)
		ax = Axis(fig[1, 1], xlabel = "y [km]", ylabel = "Depth [km]", yticks = ([-3, -2, -1, 0], ["-3", "-2", "-1", "0"]))
		hm = contourf!(ax, L, zC ./ 1e3, v11; colormap, levels) 
		xlims!(ax, (0, L[end]))
		ylims!(ax, (-3, 0))
		
		ax  = Axis(fig[1, 2], xlabel = "y [km]", ylabel = "", yticks = ([-3, -2, -1, 0], ["", "", "", ""]))
		hm = contourf!(ax, L, zC./ 1e3, v21; colormap, levels) 
		xlims!(ax, (0, L[end]))
		ylims!(ax, (-3, 0))
		
		ax  = Axis(fig[1, 3], xlabel = "y [km]", ylabel = "", yticks = ([-3, -2, -1, 0], ["", "", "", ""]))
		hm = contourf!(ax, L, zC./ 1e3, v31; colormap, levels) 
		xlims!(ax, (0, L[end]))
		ylims!(ax, (-3, 0))
		cb = Colorbar(fig[0, 1:3], hm, vertical = false, label = title, width = Relative(2/3))
	
		D1 = mean(v11, dims = 1)[1, :]
		D2 = mean(v21, dims = 1)[1, :]
		D3 = mean(v31, dims = 1)[1, :]

		ax = Axis(fig[1, 4]; xlabel = "Meridional - mean " * title, ylabel = "Depth [km]")
		lines!(ax, D1, zC ./ 1e3, label = "Case 1")
		lines!(ax, D2, zC ./ 1e3, linestyle = :dash, label = "Case 2")
		lines!(ax, D3, zC ./ 1e3, linestyle = :dashdot, label = "Case 3")
		axislegend(; framecolor = :transparent, position = labpos, backgroundcolor = :transparent)
		if !isnothing(fxlim)
			xlims!(ax, fxlim)
		end
		
		save(title * ".png", fig)

		return fig
	end
	
	# A utility function to read the variables from MITgcm
	function read_variable(file, Nvar; Nx = 200, Ny = 400, Nz = 90)
     	var = Array{Float32}(undef, Nx*Ny*Nz*Nvar)
     	read!(file, var)
     	var = bswap.(var) |> Array{Float64}
     	var = reshape(var, Nx, Ny, Nz, Nvar)
     	var = reverse(var, dims = 3)

    	return var
	end
end

# ╔═╡ be071ec1-0067-41de-b13f-2c185b117fd1
md"""
# Simulation:

Re-entrant channel with prescribed zonal wind stress and heat-flux and restore to an exponentially stratified profile at the northern boundary.
The grid spacing is 5 kilometers in the horizontal and variable vertical spacing with 90 levels, finer near the top and the bottom.
All the simulations are run for 40 years starting from the same initial condition. All outputs are averaged over five years.

6 cases are run: 
There is a choice between momentum case 0 and 1 and tracer case 0, 1, 2, and 3

### Momentum choice:
- 0 -> WENO vector invariant, 9th order for vorticity, 5th for the rest
- 1 -> Centered advection with biharmonic viscosity and ``\nu = 9e8``

### Tracer choice:
- 0 -> WENO 7th order in all 3 directions
- 1 -> WENO 5th order in ``x`` and ``y`` and Centered in ``z``
- 2 -> WENO 9th order in ``x`` and ``y`` and Centered in ``z``
- 3 -> Upwind 3rd order in all 3 directions

Case 1:
- momentum $(@bind mom_case1 Select([1, 0])) 
- tracer   $(@bind tra_case1 Select([3, 0, 1, 2, 4]))

Case 2:
- momentum $(@bind mom_case2 Select([1, 0]))
- tracer   $(@bind tra_case2 Select([2, 0, 1, 3, 4]))

Case 3:
- momentum $(@bind mom_case3 Select([0, 1, 2, 3]))
- tracer   $(@bind tra_case3 Select([2, 0, 1, 3, 4]))

"""

# ╔═╡ b40130ad-e33f-4b99-8394-2bdcaebc8567
begin
	mom1  = string(mom_case1)
	tra1  = string(tra_case1)
	case1 = parse(Int, mom1 * tra1)
	
	mom2  = string(mom_case2)
	tra2  = string(tra_case2)
	case2 = parse(Int, mom2 * tra2)

	mom3  = string(mom_case3)
	tra3  = string(tra_case3)
	case3 = parse(Int, mom3 * tra3)
	
	tracer_advection1 = if tra1 == "0"
		"WENO 7th order in all 3 directions"
	elseif tra1 == "1"
		"WENO 5th order in x and y and Centered in z"
	elseif tra1 == "2"
		"WENO 9th order in x and y and Centered in z"
	else
		"Upwind third order"
	end

	tracer_advection2 = if tra2 == "0"
		"WENO 7th order in all 3 directions"
	elseif tra2 == "1"
		"WENO 5th order in x and y and Centered in z"
	elseif tra2 == "2"
		"WENO 9th order in x and y and Centered in z"
	else
		"Upwind third order"
	end

	tracer_advection3 = if tra3 == "0"
		"WENO 7th order in all 3 directions"
	elseif tra3 == "1"
		"WENO 5th order in x and y and Centered in z"
	elseif tra3 == "2"
		"WENO 9th order in x and y and Centered in z"
	else
		"Upwind third order"
	end
	
	momentum1 = if mom1 == "0"
		"WENO vector invariant, 9th order for vorticity, 5th for the rest"
	else
		"Centered advection with biharmonic viscosity and ν = 9e8"
	end
	
	momentum2 = if mom2 == "0"
		"WENO vector invariant, 9th order for vorticity, 5th for the rest"
	else
		"Centered advection with biharmonic viscosity and ν = 9e8"
	end
	
	momentum3 = if mom3 == "0"
		"WENO vector invariant, 9th order for vorticity, 5th for the rest"
	else
		"Centered advection with biharmonic viscosity and ν = 9e8"
	end
	
	md"""
	# Numerical Details

	### Case 1

	- momentum: $(momentum1)
	- closure: convective adjustment with ``\kappa = \nu = 0.1`` and background viscosity of ``3e-4``
	- tracer: $(tracer_advection1)
	
	### Case 2

	- momentum: $(momentum2)
	- closure: convective adjustment with ``\kappa = \nu = 0.1`` and background viscosity of ``3e-4``
	- tracer: $(tracer_advection2)
	
	### Case 2

	- momentum: $(momentum3)
	- closure: convective adjustment with ``\kappa = \nu = 0.1`` and background viscosity of ``3e-4``
	- tracer: $(tracer_advection3)
	"""
end

# ╔═╡ 01b9f5ae-ec5e-4190-a795-43f6d23744cd
begin
	lev = 90
	
	file1s = jldopen("../channel_snapshots_" * string(case1) * ".jld2")
	iter1s = keys(file1s["timeseries/t"])[end]
	bo1s   = file1s["timeseries/b/" * iter1s][:, :, lev] ./ 2e-4 ./ 9.80665
	uo1s   = file1s["timeseries/u/" * iter1s][:, :, lev]
	vo1s   = file1s["timeseries/v/" * iter1s][:, :, lev]
	wo1s   = file1s["timeseries/w/" * iter1s][:, :, 70] .* 1000
	
	file2s = jldopen("../channel_snapshots_" * string(case2) * ".jld2")
	iter2s = keys(file2s["timeseries/t"])[end]
	bo2s   = file2s["timeseries/b/" * iter2s][:, :, lev] ./ 2e-4 ./ 9.80665
	uo2s   = file2s["timeseries/u/" * iter2s][:, :, lev]
	vo2s   = file2s["timeseries/v/" * iter2s][:, :, lev]
	wo2s   = file2s["timeseries/w/" * iter2s][:, :, 70] .* 1000

	file3s = jldopen("../channel_snapshots_" * string(case3) * ".jld2")
	iter3s = keys(file3s["timeseries/t"])[end]
	bo3s   = file3s["timeseries/b/" * iter3s][:, :, 90] ./ 2e-4 ./ 9.80665
	uo3s   = file3s["timeseries/u/" * iter3s][:, :, 90]
	vo3s   = file3s["timeseries/v/" * iter3s][:, :, 90]
	wo3s   = file3s["timeseries/w/" * iter3s][:, :, 70] .* 1000

	figbi = instantaneous_heatmaps(bo1s, bo2s, bo3s, :magma,  "Temperature"; fcrange = (2, 7.8))
	figui = instantaneous_heatmaps(uo1s, uo2s, uo3s, :viridis, "u-vel"; fcrange=(-0.5, 1.2))
	figvi = instantaneous_heatmaps(vo1s, vo2s, vo3s, :viridis, "v-vel"; fcrange=(-0.5, 0.5))
	figwi = instantaneous_heatmaps(wo1s, wo2s, wo3s, :viridis, "w-vel"; fcrange=(-4e-1, 4e-1))

	md"""
	# INSTANTANEOUS FLOW RESULTS
	
	comparing the instantaneous snapshots of buoyancy and velocity between Case 1 and Case 2,
	where the average is performed over the zonal direction and 5 years time.

	$(figbi)
	$(figui)
	$(figvi)
	$(figwi)
	"""
end

# ╔═╡ dd3df2d7-6de1-4c4f-ab19-a5154e53e953
begin
	file1 = jldopen("../channel_averages_" * string(case1) * ".jld2")
	iter1 = keys(file1["timeseries/t"])[end]
	bo1   = mean(file1["timeseries/b/" * iter1], dims = 1)[1, :, :]
	uo1   = mean(file1["timeseries/u/" * iter1], dims = 1)[1, :, :]
	vo1   = mean(file1["timeseries/v/" * iter1], dims = 1)[1, 1:end-1, :]
	
	file2 = jldopen("../channel_averages_" * string(case2) * ".jld2")
	iter2 = keys(file2["timeseries/t"])[end]
	bo2   = mean(file2["timeseries/b/" * iter2], dims = 1)[1, :, :]
	uo2   = mean(file2["timeseries/u/" * iter2], dims = 1)[1, :, :]
	vo2   = mean(file2["timeseries/v/" * iter2], dims = 1)[1, 1:end-1, :]
	
	file3 = jldopen("../channel_averages_" * string(case3) * ".jld2")
	iter3 = keys(file3["timeseries/t"])[end]
	bo3   = mean(file3["timeseries/b/" * iter3], dims = 1)[1, :, :]
	uo3   = mean(file3["timeseries/u/" * iter3], dims = 1)[1, :, :]
	vo3   = mean(file3["timeseries/v/" * iter3], dims = 1)[1, 1:end-1, :]
	
	md"""
	# MEAN FLOW RESULTS
	
	comparing the mean buoyancy and velocity between Case 1 and Case 2,
	where the average is performed over the zonal direction and 5 years time.
	
	Buoyancy 
	- Red  -> Case 1 hotter
	- Blue -> Case 2 hotter
	
	Velocity
	- Red - > Case 1 faster
	- Blue -> Case 2 faster
	"""
end

# ╔═╡ 8e4a72b1-529f-4298-85e8-655c64f8b84f
# Plotting buoyancy
plot_contours(bo1, bo2, bo3, (-1e-3, 1e-3), :thermal, "buoyancy"; levels = -0.001:0.001:0.015)	

# ╔═╡ 6ce5feaf-9a18-4514-a80c-8d71fa9d0578
plot_heatmaps(uo1, uo2, uo3, (-1e-1, 1e-1), :jet, "u-vel-mean"; fcrange = (0, 0.3))	

# ╔═╡ be76279f-f6ae-4309-8a1b-e8349b744fd2
plot_and_save(uo1, uo2, uo3, 100, "zonal velocity"; fxlim = (0.06, 0.21))

# ╔═╡ 2f70011b-9c1c-4ed0-99a2-4bbde6ed65ca
plot_heatmaps(vo1, vo2, vo3, (-1e-3, 1e-3), :jet, "v-vel-mean"; fcrange=(-0.005, 0.005)) 

# ╔═╡ e9469d80-80b9-4cdc-ae82-ae5874a9c512
begin
	Pxo1 = - mean(file1["timeseries/Pu/"   * iter1], dims = 1)[1, :, :]
	Pyo1 = - mean(file1["timeseries/Pv/"   * iter1], dims = 1)[1, :, :]
	Pzo1 = - mean(file1["timeseries/Pw/"   * iter1], dims = 1)[1, :, :]
	Bxo1 =   mean(file1["timeseries/∂xb²/" * iter1], dims = 1)[1, :, :]
	Byo1 =   mean(file1["timeseries/∂yb²/" * iter1], dims = 1)[1, :, :]
	Bzo1 =   mean(file1["timeseries/∂zb²/" * iter1], dims = 1)[1, :, :]

	Pxo2 = - mean(file2["timeseries/Pu/"   * iter2], dims = 1)[1, :, :]
	Pyo2 = - mean(file2["timeseries/Pv/"   * iter2], dims = 1)[1, :, :]
	Pzo2 = - mean(file2["timeseries/Pw/"   * iter2], dims = 1)[1, :, :]
	Bxo2 =   mean(file2["timeseries/∂xb²/" * iter2], dims = 1)[1, :, :]
	Byo2 =   mean(file2["timeseries/∂yb²/" * iter2], dims = 1)[1, :, :]
	Bzo2 =   mean(file2["timeseries/∂zb²/" * iter2], dims = 1)[1, :, :]

	Pxo3 = - mean(file3["timeseries/Pu/"   * iter3], dims = 1)[1, :, :]
	Pyo3 = - mean(file3["timeseries/Pv/"   * iter3], dims = 1)[1, :, :]
	Pzo3 = - mean(file3["timeseries/Pw/"   * iter3], dims = 1)[1, :, :]
	Bxo3 =   mean(file3["timeseries/∂xb²/" * iter3], dims = 1)[1, :, :]
	Byo3 =   mean(file3["timeseries/∂yb²/" * iter3], dims = 1)[1, :, :]
	Bzo3 =   mean(file3["timeseries/∂zb²/" * iter3], dims = 1)[1, :, :]

	md"""
	# VARIANCE DIAGNOSTIC RESULTS

	comparing the implicit dissipation between Case 1 and Case 2

	### Implicit dissipation in the three directions
	
	calculated as:
	```math
	P_d^n = \delta_d \left(b^{n+1} + b^n\right) \cdot \mathscr{A}_d - V_d \cdot \delta_d (b^n b^{n+1})
	```
	Where
	- subscript ``d`` indicates the direction (``x``, ``y``, or ``z``)
	- superscripts indicate the time step
	- ``\delta`` is a discrete difference
	- ``\mathscr{A}`` is the advective flux
	- ``V_d`` is the area transport in the ``d`` direction
	
	In MITgcm the advective flux is calculated directly, while in Oceananigans, that uses an AB2 time stepping procedure, the advective flux is calculated as an AB2 time-extrapolation of the advective fluxes at time-steps ``n`` and ``n-1``:
	```math
	\mathscr{A}_d = (1.5 + \chi) \cdot V_d^n \widetilde{b}^{n} - (0.5 + \chi) \cdot V_d^{n-1} \widetilde{b}^{n-1} 
	```
	(There are some shenanigans to account for a changing grid-size, but the impact of grid scaling is extremely small, basically negligible.)
	The results in the heatmaps are averaged over the zonal direction and 5 years time. The lines plot shows results averaged also in the meridional direction.
	
	Smaller is better
	- Blue -> Case 1 better
	- Red  -> Case 2 better
	"""
end

# ╔═╡ 536d6d46-c42f-485d-9216-a5c7a83132cb
plot_heatmaps(Pxo1, Pxo2, Pxo3, (-3e-6, 3e-6), :deep, "Px"; fcrange = (0, 3e-5), labpos = :rb)

# ╔═╡ 51a782c7-2f27-4f35-be05-9ab83e96f505
plot_heatmaps(Pyo1, Pyo2, Pyo3, (-3e-6, 3e-6), :deep, "Py"; fcrange = (0, 3e-5), labpos = :rb)

# ╔═╡ c54146a3-0613-4da0-a582-fca36cc25e4f
plot_heatmaps(Pzo1[:, 2:end], Pzo2[:, 2:end], Pzo3[:, 2:end], (-3e-6, 3e-6), :deep, "Pz"; fcrange = (-1e-7, 1e-6), fxlim = (-1e-6, 5e-6), labpos = :rb)

# ╔═╡ 51f63c7c-d11a-44fe-aae4-1925cee207eb
md"""
### Buoyancy gradient squared in the three directions

Calculated as
```math
B^d = \left( \frac{\partial b}{\partial d} \right)^2
```
averaged in time (5 years) and in the zonal direction

Larger is better? Not sure how to interpret this quantity, sure we don't want to destroy too much variance so might be that a smaller value indicates more implicit diffusivity.
So let's go with that and say that:
- Blue -> Case 1 worse
- Red  -> Case 2 worse
"""

# ╔═╡ 6dadc7c6-ca03-4945-84af-640144cdf3d6
plot_heatmaps(Bxo1, Bxo2, Bxo3, (-1e-7, 1e-7), :deep, "Bˣ"; fcrange = (0, 5e-7), labpos = :rb)

# ╔═╡ 9d5eb923-8bce-417f-9a9c-05c3bd1e950a
plot_heatmaps(Byo1, Byo2, Byo3, (-1e-7, 1e-7), :deep, "Bʸ"; fcrange = (0, 5e-7), labpos = :rb)

# ╔═╡ f0983ff0-4ab4-4ff2-88bf-305e3cb2bb13
plot_heatmaps(Bzo1[:, 2:end], Bzo2[:, 2:end], Bzo3[:, 2:end], (-3e-2, 3e-2), :deep, "Bᶻ"; fcrange = (0, 0.07), fxlim = (-0.01, 0.2), labpos = :rb)

# ╔═╡ 28f0b020-a05e-474a-a8f6-b475bb0847ee
begin 
	κxo1 = (Pxo1 ./ 2 ./ Bxo1)[3:end-3, :]
	κyo1 = (Pyo1 ./ 2 ./ Byo1)
	κzo1 = (Pzo1 ./ 2 ./ Bzo1)

	κyo1 = ((κyo1[1:end-1, :] .+ κyo1[2:end, :]) ./ 2)[3:end-3, :]
	κzo1 = ((κzo1[:, 1:end-1] .+ κzo1[:, 2:end]) ./ 2)[3:end-3, :]

	κxo2 = (Pxo2 ./ 2 ./ Bxo2)[3:end-3, :]
	κyo2 = (Pyo2 ./ 2 ./ Byo2)
	κzo2 = (Pzo2 ./ 2 ./ Bzo2)

	κyo2 = ((κyo2[1:end-1, :] .+ κyo2[2:end, :]) ./ 2)[3:end-3, :]
	κzo2 = ((κzo2[:, 1:end-1] .+ κzo2[:, 2:end]) ./ 2)[3:end-3, :]

	κxo3 = (Pxo3 ./ 2 ./ Bxo3)[3:end-3, :]
	κyo3 = (Pyo3 ./ 2 ./ Byo3)
	κzo3 = (Pzo3 ./ 2 ./ Bzo3)

	κyo3 = ((κyo3[1:end-1, :] .+ κyo3[2:end, :]) ./ 2)[3:end-3, :]
	κzo3 = ((κzo3[:, 1:end-1] .+ κzo3[:, 2:end]) ./ 2)[3:end-3, :]

	md"""
	### Implicit diffusivity in the three directions

	calculated as:
	```math
	\kappa_d = \frac{\langle P_d\rangle}{2 \langle(\partial_d b)^2\rangle}
	```
	where ``d`` is the direction (``x``, ``y`` or ``z``), and the right and left angle brackets indicate a zonal and time average (5 years)

	Smaller is definitely better
	- Blue -> Case 1 better
	- Red  -> MITgcm better
	"""
end

# ╔═╡ 4fe9de71-71b1-4162-beed-7a424e80b5de
plot_heatmaps(κxo1, κxo2, κxo3, (-5, 5), :jet, "κx"; fcrange = (0, 10), labpos = :rc)

# ╔═╡ 4e5f0a58-3c19-46f6-9904-8d3af3c46609
plot_heatmaps(κyo1, κyo2, κyo3, (-5, 5), :jet, "κy"; fcrange = (0, 10), labpos = :rc)

# ╔═╡ 31a75f83-c608-42ef-8acd-12999e1439d5
plot_heatmaps(κzo1, κzo2, κzo3, (-5e-5, 5e-5), :jet, "κz"; fcrange = (-1e-5, 1e-5), fxlim= (-1e-5, 5e-5), labpos = :rc)

# ╔═╡ 7b05777c-0ea7-48eb-a315-c2a020283ca6
begin
	Bzio1 = (Bzo1[:, 1:end-1] .+ Bzo1[:, 2:end]) / 2
	Pyio1 = (Pyo1[1:end-1, :] .+ Pyo1[2:end, :]) / 2
	κxio1 = (Pxo1  ./ 2 ./ Bzio1)[3:end-3, :]
	κyio1 = (Pyio1 ./ 2 ./ Bzio1)[3:end-3, :]
	κzio1 = κzo1

	κio1 = κxio1 + κzio1 + κyio1 
	
	Bzio2 = (Bzo2[:, 1:end-1] .+ Bzo2[:, 2:end]) / 2
	Pyio2 = (Pyo2[1:end-1, :] .+ Pyo2[2:end, :]) / 2
	κxio2 = (Pxo2  ./ 2 ./ Bzio2)[3:end-3, :]
	κyio2 = (Pyio2 ./ 2 ./ Bzio2)[3:end-3, :]
	κzio2 = κzo2

	κio2 = κxio2 + κzio2 + κyio2 
	
	Bzio3 = (Bzo3[:, 1:end-1] .+ Bzo3[:, 2:end]) / 2
	Pyio3 = (Pyo3[1:end-1, :] .+ Pyo3[2:end, :]) / 2
	κxio3 = (Pxo3  ./ 2 ./ Bzio3)[3:end-3, :]
	κyio3 = (Pyio3 ./ 2 ./ Bzio3)[3:end-3, :]
	κzio3 = κzo2

	κio3 = κxio3 + κzio3 + κyio3
	md"""
	### Another measure of implicit diffusivity
	
	Since ``Bz \gg Bx \approx By``, it makes sense to calculate the implicit diffusivity using the vertical buoyancy gradient squared instead of the horizontal one. We compute
	```math
	\kappa_{di} = \frac{\langle P_d \rangle}{2 \langle (\partial_z b)^2 \rangle}
	```
	and, finally,
	```math
	\kappa_i = \kappa_{xi} + \kappa_{yi} + \kappa_z
	```
	
	where ``i`` stands for _isotropic_. This diffusivity shows the value of an ipothetical vertical diffusivity that would result (in combination with a _perfect_ advection scheme) in the same level of diapycnal mixing obtained from the implicit diffusion of the advection scheme.

	Once again, smaller is definitely better
	- Blue -> Case 1 better
	- Red  -> Case 2 better
	"""
end

# ╔═╡ 99fb75f7-940d-4b95-8bdf-803ced0f6f1e
plot_heatmaps(κxio1, κxio2, κxio3, (-1e-4, 1e-4), :jet, "κxi"; fcrange = (0, 5e-5), fxlim=(-1e-6, 1e-4), labpos = :rc)

# ╔═╡ 61605dc0-b6a3-48d1-89b8-ab268992d1b9
plot_heatmaps(κyio1, κyio2, κyio3, (-1e-4, 1e-4), :jet, "κyi"; fcrange = (0, 5e-5), fxlim=(-1e-6, 1e-4), labpos = :rc)

# ╔═╡ 1f8d2ecf-a783-407a-a182-9ea8d7bf9550
md"""
# Total Isotropic Implicit Diffusivity!!

Everything culminates here...

The vertical dotted line in the plot on the left-hand-side shows the value of 1e-5, 
while the dashed-dotted line shows 4e-6.
"""

# ╔═╡ d3b3eca2-a80d-40d8-805d-6ce40e0de0b9
begin
	fig = plot_heatmaps(κio1, κio2, κio3, (-1e-4, 1e-4), :jet, "κi"; fcrange = (0, 2e-4), labpos = :rc, fxlim = (-5e-6, 1e-4))
	vlines!(1e-5, color = :grey, linestyle = :dot)
	# vlines!(4e-6, color = :grey, linestyle = :dash)
	fig
end

# ╔═╡ 9fa329c7-a31b-4216-92cc-86560bb11683
md""" ## Saving Figures """

# ╔═╡ 892f2f46-9875-479f-b67e-f00783bda4b6
plot_and_save(κxo1, κxo2, κxo3, 2, "κx"; fxlim = (-1, 15), labpos = :rc)

# ╔═╡ 7534194f-63b5-4a08-b026-80969301cd7f
plot_and_save(κyo1, κyo2, κyo3, 2, "κy"; fxlim = (-1, 15), labpos = :rc)

# ╔═╡ d236059b-b20e-4464-bf80-b30e854a4414
plot_and_save(κzo1, κzo2, κzo3, 5e-6, "κz"; fxlim = (-1e-5, 1.4e-5), labpos = :rc)

# ╔═╡ d02bf8c2-e899-4ec1-b273-f41e4ec51307
plot_and_save(κio1, κio2, κio3, 1e-5, "κi"; fxlim = (-1e-6, 5e-5), labpos = :rc)

# ╔═╡ Cell order:
# ╟─ec938fce-60b1-11ef-3e29-f5396e753494
# ╟─cdb48125-a977-4e3a-91ec-1861792dd9d1
# ╟─1c578cdd-cd83-499b-98a0-1292b566a0f8
# ╟─bc52c23b-5573-4931-93a5-5b598e8cc83b
# ╟─37f329a1-c664-4cc9-acab-c6c5d357c7da
# ╟─be071ec1-0067-41de-b13f-2c185b117fd1
# ╟─b40130ad-e33f-4b99-8394-2bdcaebc8567
# ╟─01b9f5ae-ec5e-4190-a795-43f6d23744cd
# ╟─dd3df2d7-6de1-4c4f-ab19-a5154e53e953
# ╟─8e4a72b1-529f-4298-85e8-655c64f8b84f
# ╠═6ce5feaf-9a18-4514-a80c-8d71fa9d0578
# ╠═be76279f-f6ae-4309-8a1b-e8349b744fd2
# ╟─2f70011b-9c1c-4ed0-99a2-4bbde6ed65ca
# ╟─e9469d80-80b9-4cdc-ae82-ae5874a9c512
# ╟─536d6d46-c42f-485d-9216-a5c7a83132cb
# ╟─51a782c7-2f27-4f35-be05-9ab83e96f505
# ╟─c54146a3-0613-4da0-a582-fca36cc25e4f
# ╟─51f63c7c-d11a-44fe-aae4-1925cee207eb
# ╟─6dadc7c6-ca03-4945-84af-640144cdf3d6
# ╟─9d5eb923-8bce-417f-9a9c-05c3bd1e950a
# ╟─f0983ff0-4ab4-4ff2-88bf-305e3cb2bb13
# ╟─28f0b020-a05e-474a-a8f6-b475bb0847ee
# ╟─4fe9de71-71b1-4162-beed-7a424e80b5de
# ╟─4e5f0a58-3c19-46f6-9904-8d3af3c46609
# ╟─31a75f83-c608-42ef-8acd-12999e1439d5
# ╟─7b05777c-0ea7-48eb-a315-c2a020283ca6
# ╟─99fb75f7-940d-4b95-8bdf-803ced0f6f1e
# ╟─61605dc0-b6a3-48d1-89b8-ab268992d1b9
# ╟─1f8d2ecf-a783-407a-a182-9ea8d7bf9550
# ╟─d3b3eca2-a80d-40d8-805d-6ce40e0de0b9
# ╟─9fa329c7-a31b-4216-92cc-86560bb11683
# ╟─892f2f46-9875-479f-b67e-f00783bda4b6
# ╟─7534194f-63b5-4a08-b026-80969301cd7f
# ╠═d236059b-b20e-4464-bf80-b30e854a4414
# ╠═d02bf8c2-e899-4ec1-b273-f41e4ec51307
