### A Pluto.jl notebook ###
# v0.19.46

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
	include("comparison.jl")
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
# a utility function for plotting:
function plot_heatmaps(vm, vo, colorrange, colormap, title; fcrange = nothing)

	L = 1:size(vm, 1)
	
	fig = Figure(size = (1200, 400), fontsize = 10)
	ax = Axis(fig[1, 1:2], xlabel = "Latitude", ylabel = "Depth")
	hm = isnothing(fcrange) ? heatmap!(ax, L, zC, vm; colormap) : heatmap!(ax, vm; colormap, colorrange = fcrange) 
	cb = Colorbar(fig[0, 1:2], hm, vertical = false, label = "MITgcm " * title)
	
	ax  = Axis(fig[1, 3:4], xlabel = "Latitude", ylabel = "Depth")
	hm = isnothing(fcrange) ? heatmap!(ax, L, zC, vo; colormap) : heatmap!(ax, vo; colormap, colorrange = fcrange) 
	cb = Colorbar(fig[0, 3:4], hm, vertical = false, label = "Oceananigans " * title)

	ax  = Axis(fig[1, 5:6], xlabel = "Latitude", ylabel = "Depth")
	hm = heatmap!(ax, L, zC, vm .- vo; colormap = :bwr, colorrange)
	cb = Colorbar(fig[0, 5:6], hm, vertical = false, label = "Difference")

	Dm = mean(vm, dims = 1)[1, :]
	Do = mean(vo, dims = 1)[1, :]
	
	ax = Axis(fig[0:1, 7]; xticklabelsvisible = false, yticklabelsvisible = false)
	lines!(ax, Dm, zC)
	lines!(ax, Do, zC, linestyle = :dash)
	return fig
end

# ╔═╡ be071ec1-0067-41de-b13f-2c185b117fd1
md"""
# Select Oceananigans' case
"""

# ╔═╡ b41bdba4-8419-4010-8238-5a53a8a87ce3
@bind mom_case NumberField(0:1, default=1)

# ╔═╡ bdc61d8c-9757-448f-94fa-e041e1d6e0ad
@bind tra_case NumberField(0:3, default=3)


# ╔═╡ b40130ad-e33f-4b99-8394-2bdcaebc8567
begin
	mom = string(mom_case)
	tra = string(tra_case)
	oceananigans_case = parse(Int, mom * tra)

	tracer_advection = if tra == "0"
		"WENO 7th order in all 3 directions"
	elseif tra == "1"
		"WENO 5th order in ``x`` and ``y`` and Centered in ``z``"
	elseif tra == "2"
		"WENO 9th order in ``x`` and ``y`` and Centered in ``z``"
	else
		"Upwind third order"
	end

	momentum = if mom == "0"
		"WENO vector invariant, 9th order for vorticity, 5th for the rest"
	else
		"Centered advection with biharmonic viscosity and ``\nu = 9e8``"
	end
	
	md"""
	# Numerical Details

	### MITgcm

	- momentum: Centered advection with biharmonic viscosity and ``\nu = 9e8``
	- closure: KKP and background viscosity of ``3e-4``
	- tracer: Upwind third order advection with flux limiter

	### Ocenanigans

	- momentum: $(momentum)
	- closure: convective adjustment with ``\kappa = \nu = 0.1`` and background viscosity of ``3e-4``
	- tracer: $(tracer_advection)
	"""
end

# ╔═╡ dd3df2d7-6de1-4c4f-ab19-a5154e53e953
md"""
# MEAN FLOW RESULTS

comparing the mean buoyancy and velocity between Oceananigans and MITgcm,
where the average is performed over the zonal direction and 5 years time.

Buoyancy 
- Red -> Oceananigans colder
- Blue -> MITgcm colder

Velocity
- Red -> Oceananigans slower
- Blue -> MITgcm slower
"""

# ╔═╡ 058b4798-e5d2-4ea6-ac43-32e63987ac28
# Start by reading mean results
begin 
	var1 = read_variable("dynDiag.0001036800.data", 9);
	bm  = mean(var1[:, :, :, 3],  dims = 1)[1, :, :] .* α   .* g
	b²m = mean(var1[:, :, :, 4],  dims = 1)[1, :, :] .* α^2 .* g^2 
	Um  = mean(var1[:, :, :, 8],  dims = 1)[1, :, :] 
	Vm  = mean(var1[:, :, :, 2],  dims = 1)[1, :, :] 

	file = jldopen("channel_averages_" * string(oceananigans_case) * ".jld2")
	iter = keys(file["timeseries/t"])[end]
	bo   = mean(file["timeseries/b/" * iter], dims = 1)[1, :, :]
	uo   = mean(file["timeseries/u/" * iter], dims = 1)[1, :, :]
	vo   = mean(file["timeseries/v/" * iter], dims = 1)[1, 1:end-1, :]

	nothing
	nothing
end

# ╔═╡ 8e4a72b1-529f-4298-85e8-655c64f8b84f
# Plotting buoyancy
plot_heatmaps(bm, bo, (-1e-3, 1e-3), :thermal, "buoyancy")	

# ╔═╡ 6ce5feaf-9a18-4514-a80c-8d71fa9d0578
plot_heatmaps(Um, uo, (-1e-1, 1e-1), :jet, "u-vel"; fcrange = (0, 0.3))	

# ╔═╡ 2f70011b-9c1c-4ed0-99a2-4bbde6ed65ca
plot_heatmaps(Vm, vo, (-1e-1, 1e-1), :jet, "v-vel"; fcrange=(-0.075, 0.075))	

# ╔═╡ f14d4d72-0715-4859-9daf-e846c2228f17
md"""
# VARIANCE DIAGNOSTIC RESULTS

comparing the implicit dissipation between Oceananigans and MITgcm
"""

# ╔═╡ dccdd6b3-89ca-4816-b143-d0c0d45ef5d7
# Start reading the MITgcm results
begin 
	var2 = read_variable("varDiag.0001036800.data", 12);
	Pxm = - mean(var2[:, :, :, 5],  dims = 1)[1, :, :] .* α^2 .* g^2
	Pym = - mean(var2[:, :, :, 6],  dims = 1)[1, :, :] .* α^2 .* g^2
	Pzm = - mean(var2[:, :, :, 7],  dims = 1)[1, :, :] .* α^2 .* g^2
	Bxm =   mean(var2[:, :, :, 9],  dims = 1)[1, :, :] .* α^2 .* g^2 ./ Δh
	Bym =   mean(var2[:, :, :, 10], dims = 1)[1, :, :] .* α^2 .* g^2 ./ Δh 
	Bzm =   mean(var2[:, :, :, 11], dims = 1)[1, :, :] .* α^2 .* g^2

	for k in 1:90
		Bzm[:, k] ./= Δz[k]
	end
	
	Pxo = - mean(file["timeseries/Pu/"   * iter], dims = 1)[1, :, :]
	Pyo = - mean(file["timeseries/Pv/"   * iter], dims = 1)[1, :, :]
	Pzo = - mean(file["timeseries/Pw/"   * iter], dims = 1)[1, :, :]
	Bxo =   mean(file["timeseries/∂xb²/" * iter], dims = 1)[1, :, :]
	Byo =   mean(file["timeseries/∂yb²/" * iter], dims = 1)[1, :, :]
	Bzo =   mean(file["timeseries/∂zb²/" * iter], dims = 1)[1, :, :]
	nothing
end

# ╔═╡ e9469d80-80b9-4cdc-ae82-ae5874a9c512
md"""
### Implicit dissipation in the three directions

calculated as:
```math
P_d^n = \delta_d \left(\frac{b^{n+1} + b^n}{2}\right) \cdot \mathscr{A}_d - V_d \cdot \delta_d (b^n b^{n+1})
```
Where
- subscript ``d`` indicates the direction
- superscripts indicate the time step
- ``\delta`` is a discrete difference
- ``\mathscr{A}`` is the advective flux
- ``V_d`` is the area transport in the ``d`` direction

In MITgcm the advective flux is calculated directly, while in Oceananigans, that uses an AB2 time stepping procedure, the advective flux is calculated as an AB2 time-extrapolation of the advective fluxes at time-steps ``n`` and ``n-1``:
```math
\mathscr{A}_d = (1.5 + \chi) \cdot V_d^n \widetilde{b}^{n} - (0.5 + \chi) \cdot V_d^{n-1} \widetilde{b}^{n-1} 
```

and averaged over the zonal direction and 5 years time.

Smaller is better
- Red  -> Oceananigans better
- Blue -> MITgcm better
"""

# ╔═╡ 536d6d46-c42f-485d-9216-a5c7a83132cb
plot_heatmaps(Pxm, Pxo, (-3e-6, 3e-6), :deep, "Px"; fcrange = (0, 3e-5))

# ╔═╡ 51a782c7-2f27-4f35-be05-9ab83e96f505
plot_heatmaps(Pym, Pyo[1:end-1, :], (-3e-6, 3e-6), :deep, "Py"; fcrange = (0, 3e-5))

# ╔═╡ c54146a3-0613-4da0-a582-fca36cc25e4f
plot_heatmaps(Pzm, Pzo[:, 1:end-1], (-3e-6, 3e-6), :deep, "Pz"; fcrange = (0, 1e-6))

# ╔═╡ 51f63c7c-d11a-44fe-aae4-1925cee207eb
md"""
### Buoyancy gradient squared in the three directions

Calculated as
```math
B^d = \left( \frac{\partial b}{\partial d} \right)^2
```
and averaged in time (5 years) and in the zonal direction

Larger is better?
- Red  -> Oceananigans worse
- Blue -> MITgcm worse
"""

# ╔═╡ 6dadc7c6-ca03-4945-84af-640144cdf3d6
plot_heatmaps(Bxm, Bxo, (-3e-7, 3e-7), :deep, "Bˣ"; fcrange = (0, 5e-7))

# ╔═╡ 9d5eb923-8bce-417f-9a9c-05c3bd1e950a
plot_heatmaps(Bym, Byo[1:end-1, :], (-3e-7, 3e-7), :deep, "Bʸ"; fcrange = (0, 5e-7))

# ╔═╡ f0983ff0-4ab4-4ff2-88bf-305e3cb2bb13
plot_heatmaps(Bzm, Bzo[:, 1:end-1], (-3e-2, 3e-2), :deep, "Bᶻ"; fcrange = (0, 0.07))

# ╔═╡ 28f0b020-a05e-474a-a8f6-b475bb0847ee
begin 
	κxo = (Pxo ./ 2 ./ Bxo)[3:end-3, :]
	κyo = (Pyo ./ 2 ./ Byo)
	κzo = (Pzo ./ 2 ./ Bzo)

	κyo = ((κyo[1:end-1, :] .+ κyo[2:end, :]) ./ 2)[3:end-3, :]
	κzo = ((κzo[:, 1:end-1] .+ κzo[:, 2:end]) ./ 2)[3:end-3, :]

	κxm = (Pxm ./ 2 ./ Bxm)[3:end-3, :]
	κym = (Pym ./ 2 ./ Bym)[3:end-3, :]
	κzm = (Pzm ./ 2 ./ Bzm)[3:end-3, :]

	md"""
	### Implicit diffusivity in the three directions

	calculated as:
	```math
	\kappa_d = \frac{\langle P_d\rangle}{2 \langle(\partial_d b)^2\rangle}
	```
	where ``d`` is the direction (``x``, ``y`` or ``z``), and the right and left angle brackets indicate a zonal and time average (5 years)

	Smaller is definitely better
	- Red  -> Oceananigans better
	- Blue -> MITgcm better
	"""
end

# ╔═╡ 4fe9de71-71b1-4162-beed-7a424e80b5de
plot_heatmaps(κxm, κxo, (-5, 5), :jet, "κx"; fcrange = (0, 30))

# ╔═╡ 4e5f0a58-3c19-46f6-9904-8d3af3c46609
plot_heatmaps(κym, κyo, (-5, 5), :jet, "κy"; fcrange = (0, 30))

# ╔═╡ 31a75f83-c608-42ef-8acd-12999e1439d5
plot_heatmaps(κzm, κzo, (-5e-5, 5e-5), :jet, "κz"; fcrange = (-1e-5, 1e-5))

# ╔═╡ 7b05777c-0ea7-48eb-a315-c2a020283ca6
begin
	Bzio = (Bzo[:, 1:end-1] .+ Bzo[:, 2:end]) / 2
	Pyio = (Pyo[1:end-1, :] .+ Pyo[2:end, :]) / 2
	κxio = (Pxo  ./ 2 ./ Bzio)[3:end-3, :]
	κyio = (Pyio ./ 2 ./ Bzio)[3:end-3, :]
	κzio = κzo

	κio = κxio + κyio + κzio
	
	κxim = (Pxm ./ 2 ./ Bzm)[3:end-3, :]
	κyim = (Pym ./ 2 ./ Bzm)[3:end-3, :]
	κzim = κzm

	κim = κxim + κyim + κzim
	
	md"""
	### Another measure of implicit diffusivity
	
	Since ``Bz \gg Bx \approx By``, it makes sense to calculate the implicit diffusivity as compared to the vertical bouyancy gradient squared. We compute
	```math
	\kappa_{di} = \frac{\langle P_d \rangle}{2 \langle (\partial_z b)^2 \rangle}
	```
	and, finally,
	```math
	\kappa_i = \kappa_{xi} + \kappa_{yi} + \kappa_z
	```
	
	where ``i`` stands for _isotropic_. This diffusivity shows the value of an ipothetical vertical diffusivity that would result (in combination with a _perfect_ advection scheme) in the same level of diapycnal mixing obtained from the implicit diffusion of the advection scheme.

	Once again, smaller is definitely better
	- Red  -> Oceananigans better
	- Blue -> MITgcm better
	"""
end

# ╔═╡ 99fb75f7-940d-4b95-8bdf-803ced0f6f1e
plot_heatmaps(κxim, κxio, (-1e-4, 1e-4), :jet, "κxi"; fcrange = (0, 2e-4))

# ╔═╡ 61605dc0-b6a3-48d1-89b8-ab268992d1b9
plot_heatmaps(κyim, κyio, (-1e-4, 1e-4), :jet, "κxi"; fcrange = (0, 2e-4))

# ╔═╡ 0afc2db0-53a1-47bf-9f80-55ea210660e9
md"""
# Meridional averages!
"""

# ╔═╡ Cell order:
# ╟─ec938fce-60b1-11ef-3e29-f5396e753494
# ╟─cdb48125-a977-4e3a-91ec-1861792dd9d1
# ╟─1c578cdd-cd83-499b-98a0-1292b566a0f8
# ╟─bc52c23b-5573-4931-93a5-5b598e8cc83b
# ╠═37f329a1-c664-4cc9-acab-c6c5d357c7da
# ╟─be071ec1-0067-41de-b13f-2c185b117fd1
# ╟─b41bdba4-8419-4010-8238-5a53a8a87ce3
# ╟─bdc61d8c-9757-448f-94fa-e041e1d6e0ad
# ╟─b40130ad-e33f-4b99-8394-2bdcaebc8567
# ╟─dd3df2d7-6de1-4c4f-ab19-a5154e53e953
# ╟─058b4798-e5d2-4ea6-ac43-32e63987ac28
# ╟─8e4a72b1-529f-4298-85e8-655c64f8b84f
# ╟─6ce5feaf-9a18-4514-a80c-8d71fa9d0578
# ╟─2f70011b-9c1c-4ed0-99a2-4bbde6ed65ca
# ╟─f14d4d72-0715-4859-9daf-e846c2228f17
# ╟─dccdd6b3-89ca-4816-b143-d0c0d45ef5d7
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
# ╠═0afc2db0-53a1-47bf-9f80-55ea210660e9
