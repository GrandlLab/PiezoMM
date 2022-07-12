### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ b060a110-ea87-11ec-0865-8d4ee9c6f216
begin
	using Plots, StatsPlots, PlutoUI, LinearAlgebra, Images, Parameters, Pkg, DataFrames, Measures, Impute, Arrow

	Pkg.develop("PiezoMM")

	using PiezoMM
	
	TableOfContents()
end

# ╔═╡ 0a55f3eb-64d8-48c4-935a-075ed82a6e2a
function simChannels2(channels, size, pInit, rStim, stimdf, τ;
	openState=2, σ₅₀=2.7, slope=0.8, D=240, seed=62)

	stim = fill(0.0, (length(channels), length(stimdf.t)))

	peakIdx = argmax(tensionMap(size, τ, rStim, maximum(stimdf.tension); D=D, clamp=false))

	@inbounds Threads.@threads for t ∈ 1:length(stimdf.t)
		if stimdf.t[t] == 0.0
			stim[:, t] = zeros(length(channels))
		else
			peak = getTension(Tuple(peakIdx) .- size/2, stimdf.t[t], rStim, D)
			
			@inbounds for c ∈ 1:length(channels)
				if sqrt(channels[c][1]^2 + channels[c][2]^2) <= rStim
					stim[c,t] = stimdf.tension[t]
				else
					temp = getTension(channels[c], stimdf.t[t], rStim, D)/peak
					stim[c, t] = temp * stimdf.tension[t]
				end
			end
		end
	end

	sim = fill(0.0, (length(channels), length(stimdf.t)))
	@inbounds Threads.@threads for c ∈ 1:length(channels)
		sim[c, :] = markovSequence(pInit, length(stimdf.t); σ=stim[c, :], σ₅₀=σ₅₀, slope=slope, τ=τ, seed=seed)
	end
		
	return sim
end

# ╔═╡ cec58c72-fc38-4342-abef-0150ff55777d
md"""
# Define Model Parameters

Model parameters are saved in the structure `mdl`. Any changes in the model should only be made here.
"""

# ╔═╡ a2f5420f-baae-4f44-baf5-d1f105dab17b
begin
	@with_kw struct mdl
        unit::Int = 10 ## pxl size in nm²
        D::Float64 = 2400.0 ## diffusion coefficient in px²/s Shi 2018
        τ::Float64 = 1e-6 ## in seconds
		rₛₜᵢₘ::Float64 = 200.0 ## radius in pxls
		σₘₐₓ::Float64 = 7.0 ## tension in mN/m
		size::Int = 600 ## size of grid in pxls
		tₛᵢₘ::Float64 = 0.2 ## simulation length in seconds
		σ₅₀::Float64 = 2.7 ## T50 in mN/m
		slope::Float64 = 0.8 ## slope factor
    end
end;

# ╔═╡ 539771cc-253d-491d-b11a-6c03b5c43c8d
md""" 
# Seed channels on membrane 

Channels are randomly placed on a grid with 10 nm² grid elements such that there can be no channel overlap. Briefly, this is done by at each point sampling a Binomial distribution with a probability of occupancy of 0.01. This probability approximately represents channel density estimates for channel overexpression from Lewis 2021.
"""

# ╔═╡ f4a7d99b-a1b5-452e-a198-2ff7d994f686
begin 
	channelLocs = seedChannels(fill(0, (mdl().size, mdl().size)), occupancy = 0.01, seed=30)
	
	gr()
	hm = heatmap(1:mdl().size, 1:mdl().size, channelLocs, aspect_ratio=:equal,
		yticks=false, yaxis=false, size=(600,600), cbar=false, xticks=false, color=:viridis, title="Channel Locations")

	savefig("C:/Users/HAL/Desktop/manuscript/figures/channels.png")
end

# ╔═╡ f8475401-b5e5-4b9f-9bbf-60597de253ab
# ╠═╡ disabled = true
#=╠═╡
begin
	using HDF5
	fid = h5open("C:/Users/HAL/Desktop/manuscript/source_data/channels.h5", "w")
	fid["channels"] = channelLocs
	close(fid)
end
  ╠═╡ =#

# ╔═╡ d7c19e9b-cdb8-4bd7-85ba-864581cd5800
md""" 
# Simulate Diffusion of Tension 

Tension diffuses following a normalized gaussian function of the form:

$P(r,t)=\frac{e^{\frac{-r^2}{4Dt}}}{4 \pi Dt}$

The diffusion constant $D$ is 0.24 μm²/s. This is 10x faster than determined in Shi 2018 but still 2 orders of magnitude below the estimated diffusion coefficient in axons from Shi 2020. The radius ($r$) determines the stimulus region.
"""

# ╔═╡ 4ecbf644-2c14-4226-b4d3-6167d270e207
begin
	membrane = fill(0.0, (mdl().size, mdl().size, Int64(mdl().tₛᵢₘ*1000)))
	for t in 1:Int64(mdl().tₛᵢₘ*1000)
		membrane[:,:,t] = tensionMap(mdl().size, t*1e-3, mdl().rₛₜᵢₘ, mdl().σₘₐₓ ; D=mdl().D)
	end

	gr()
	@gif for i ∈ 1:1:size(membrane, 3)
	   heatmap(1:mdl().size, 1:mdl().size, membrane[:,:,i],
		   aspect_ratio=:equal, axis=false, title=string(i) * "ms", ticks=false,
		   grid=false, clims=(0,mdl().σₘₐₓ))
	end
end

# ╔═╡ eb09ff0c-23df-4b38-bc54-27ad4106952c
md"""
# Discrete Markov model of channel gating

Each channel is simulated using a discrete markov model. The model is derived from a previously published 4-state kinetic model for Piezo1 ion channels (Lewis 2017). The model is shown in the figure below.
"""

# ╔═╡ 8941871b-17a8-4834-b2af-25a1421ede57
begin
	mdlImage = load("C:/Users/HAL/Desktop/manuscript/kinetic_model.png")
	mdlImage
end

# ╔═╡ d2b47f44-f632-41d6-a986-f02d8a21d501
md"""
## Determine equilibrium occupancies

The transition matrix `T` is determined using the function `generateTMatrix` and the equilibrium state occupancies are calculated using the function `equilibriumState`.

"""

# ╔═╡ 510c31ac-0a22-44c1-a468-2f5bcdf6d54a
begin
	T = generateTMatrix(0.0, mdl().σ₅₀, mdl().slope, mdl().τ)
	Eq = equilibriumState(T)
end

# ╔═╡ 7d0856e1-db55-4eff-b898-c547aef62956
begin
	function iterateList(df, sList, dList, τ, m)
		
		count=1
		for d in dList
			for seed in sList
				channels = seedChannels(fill(0, (mdl().size, mdl().size)), occupancy = 0.01, seed=seed)
		
				clocs = [[c[1] - mdl().size/2, c[2] - mdl().size/2] for c in Tuple.(findall(x -> x == 1.0, channelLocs))]
					
				s = simChannels(clocs, mdl().size, Eq', τ:τ:mdl().tₛᵢₘ,
					mdl().rₛₜᵢₘ, mdl().σₘₐₓ, τ, D=d, seed=seed, σ₅₀=mdl().σ₅₀)
					peakDelay = (Tuple(argmax(sum(s.==2.0, dims=1)'))[1]-0.05/τ) * τ
		
				df.diffusion[count] = d
				df.seed[count] = seed
				df.delay[count] = peakDelay
		
				count += 1
			end
			
			return df
		end
	end

	function saveExamplePlots(seed, tList, τ, m)
		for t in tList
			channels = seedChannels(fill(0, (mdl().size, mdl().size)), occupancy = 0.01, seed=seed)
		
			clocs = [[c[1] - mdl().size/2, c[2] - mdl().size/2] for c in Tuple.(findall(x -> x == 1.0, channelLocs))]
					
			s = simChannels(clocs, mdl().size, Eq', τ:τ:mdl().tₛᵢₘ,
					mdl().rₛₜᵢₘ, t, τ, D=mdl().D, seed=seed, σ₅₀=mdl().σ₅₀)
					peakDelay = (Tuple(argmax(sum(s.==2.0, dims=1)'))[1]-0.05/τ) * τ
	
			bl = @layout[a ; b]
			ptest = plot(
				collect(1:size(s,2)) .* τ, 
				sum(s.==1.0, dims=1)' ./ size(s, 1), 
				label="closed", 
				color=:red, 
				grid=false, 
				ylab="Channel Number", 
				xlab="Time (s)", 
				tick_direction=:out, 
				ylims=(0, 1.0)
			)
				
			plot!(
				collect(1:size(s,2)) .* τ, 
				sum(s.==2.0, dims=1)' ./ size(s, 1), 
				label="open", 
				color=:blue
			)
				
			plot!(
				collect(1:size(s,2)) .* τ, 
				sum(s.==3.0, dims=1)' ./ size(s, 1), 
				label="inactive_1", 
				color=:orange
			)
				
			plot!(
				collect(1:size(s,2)) .* τ, 
				sum(s.==4.0, dims=1)' ./ size(s, 1), 
				label="inactive_2", 
				color=:green
			)
			
			bmin = @layout[a{0.2h}; b]
				
			ptest1 = plot(
				collect(1:size(s,2)), 
				vcat(zeros(Int64(round(0.05/τ))), repeat([5.0], length(collect(τ:τ:mdl().tₛᵢₘ))), zeros(Int64(round(0.05/τ)))), 
				color=:maroon, 
				xaxis=false, 
				xticks=false, 
				ylab="", 
				yticks=false,
				yaxis=false,
				label=""
			)
		
			
			plot!([0.01, 0.01], [3.0, 5.0], color=:black, label="", linewidth=2)
			annotate!()
					
			ptest2 = plot(
				collect(1:size(s,2)) .* τ, 
				-2.2 .* sum(s.==2.0, dims=1)', 
				color=:black, 
				xlab="", 
				ylab="", 
				label="",
				xticks=false,
				yticks=false,
				xaxis=false,
				yaxis=false,
				title="delay= " * string(round(peakDelay * 1000, digits=2)) * " ms"
			)
			
			vline!([0.05], linestyle=:dash, color=:red, legend=:none)
			plot!([0.01, 0.01], [-500.0, -400.0], color=:black, label="", linewidth=2)
			plot!([0.01, 0.02], [-500.0, -500.0], color=:black, label="", linewidth=2)
			
			pt = plot(ptest1, ptest2, layout=bmin, grid=false, tick_direction=:out)
			
			plot(ptest, pt, layout=bl, size=(600, 800), left_margin=5mm, bottom_margin=5mm)
			
			savefig("C:/Users/HAL/Desktop/example_" * string(t) * ".png")
		end
	end
end

# ╔═╡ afd8ae76-9907-4b6a-baf8-a2d04d6d6fcf
md"""
## Simulate Channels

Simulate the current generated given the markov model of channel activity and tension diffusion as previously defined.
"""

# ╔═╡ 6e1565a1-2803-4780-af05-45dc159c24c5
sim = simChannels(findall(x -> x == 1, channelLocs), 
	mdl().size, Eq', mdl().τ:mdl().τ:mdl().tₛᵢₘ, 
	mdl().rₛₜᵢₘ, mdl().σₘₐₓ, mdl().τ; openState=2, σ₅₀=mdl().σ₅₀, slope=mdl().slope, D=mdl().D, seed=62)

# ╔═╡ 02df5141-e965-4c14-ad00-c433b5f12e30
begin
	lR = @layout [a{0.3h}; b]
	l = @layout [a b]

	p1 = plot(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1),
		sum(x ->x == 1.0, sim, dims=1)'/size(sim, 1), 
		color=:blue, label="closed", ylims=(0,1.0), grid=false, tick_direction=:out, ylab="State occupancy (%)", xlab="Time (s)")
	plot!(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1),
		sum(x ->x == 2.0, sim, dims=1)'/size(sim, 1), 
		color=:black, label="open")
	plot!(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1),
		sum(x ->x == 3.0, sim, dims=1)'/size(sim, 1), 
		color=:red, label="inactive₁")
	plot!(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1),
		sum(x ->x == 4.0, sim, dims=1)'/size(sim, 1), 
		color=:purple, label="inactive₂")

	p2 = plot(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1), vcat(zeros(Int64(0.05/mdl().τ)),
		repeat([mdl().σₘₐₓ], Int64(mdl().tₛᵢₘ/mdl().τ)), zeros(Int64(0.05/mdl().τ))), color=:maroon, xaxis=false, xticks=false, label="")

	p3 = plot(mdl().τ:mdl().τ:(mdl().tₛᵢₘ + 0.1),
		-2.2 .* sum(x ->x == 2.0, sim, dims=1)', color=:black, label="", ylab="Current (pA)", xlab="Time (s)", xlims=(0,0.3), ylims=(-600, 0))

	pR = plot(p2, p3, layout=lR, grid=false)

	pF = plot(p1, pR, layout=l, size=(700,300))
end

# ╔═╡ c1ae762d-353a-4143-9116-4c9b7b388766
# ╠═╡ disabled = true
#=╠═╡
begin
	ex = DataFrame(Arrow.Table("E:/Research/Thesis/bychannel/mp1/20191213 (done)/c1/20191213_hek293t_mp1_c1_scan-80_preprocessed.feather"))
	ex = ex[Int64(4*size(ex, 1)/5) + 1:end, :]
	ex = filter(x -> x.ti > 450 && x.ti <= 950, ex)
	oversample = Int64(round(((ex.ti[2] - ex.ti[1])/1000)/1e-5))
	stim = ex.force ./ (maximum(ex.force)/5.0)

	dftemp = DataFrame(:t => 1e-5:1e-5:0.5, :tension => [[missing, 0.0][Int64(i%2) + 1] for i in 1:length(1e-5:1e-5:0.5)])

	dftemp.tension .= missing
	@views dftemp.tension[4:4:50000] = stim 

	test = Impute.interp(dftemp) |> Impute.locf() |> Impute.nocb()

	membrane2 = fill(0.0, (mdl().size, mdl().size, length(100:100:length(test.t))))
	for t in 100:100:length(test.t)
		membrane2[:,:,div(t, 100)] = tensionMap(mdl().size, test.t[t], mdl().rₛₜᵢₘ, test.tension[t] ; D=mdl().D)
	end
	gr()
	@gif for i ∈ 1:size(membrane2, 3)
	   heatmap(1:mdl().size, 1:mdl().size, membrane2[:,:,i],
		   aspect_ratio=:equal, axis=false, title=string(i) * "ms", ticks=false,
		   grid=false, clims=(0,mdl().σₘₐₓ))
	end
end
  ╠═╡ =#

# ╔═╡ c50322ae-0dd7-4ad2-b308-0242a6cca528
# ╠═╡ disabled = true
#=╠═╡
begin 

	channels2 = seedChannels(fill(0, (mdl().size, mdl().size)), occupancy = 0.01, seed=10)
		
	clocs2 = [[c[1] - mdl().size/2, c[2] - mdl().size/2] for c in Tuple.(findall(x -> x == 1.0, channels2))]
					
	s2 = simChannels2(clocs2, mdl().size, Eq',
		mdl().rₛₜᵢₘ, test, mdl().τ, D=mdl().D, seed=10, σ₅₀=mdl().σ₅₀)
	
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═b060a110-ea87-11ec-0865-8d4ee9c6f216
# ╟─0a55f3eb-64d8-48c4-935a-075ed82a6e2a
# ╟─cec58c72-fc38-4342-abef-0150ff55777d
# ╠═a2f5420f-baae-4f44-baf5-d1f105dab17b
# ╟─539771cc-253d-491d-b11a-6c03b5c43c8d
# ╠═f4a7d99b-a1b5-452e-a198-2ff7d994f686
# ╟─d7c19e9b-cdb8-4bd7-85ba-864581cd5800
# ╠═4ecbf644-2c14-4226-b4d3-6167d270e207
# ╟─f8475401-b5e5-4b9f-9bbf-60597de253ab
# ╟─eb09ff0c-23df-4b38-bc54-27ad4106952c
# ╟─8941871b-17a8-4834-b2af-25a1421ede57
# ╟─d2b47f44-f632-41d6-a986-f02d8a21d501
# ╠═510c31ac-0a22-44c1-a468-2f5bcdf6d54a
# ╟─7d0856e1-db55-4eff-b898-c547aef62956
# ╟─afd8ae76-9907-4b6a-baf8-a2d04d6d6fcf
# ╠═6e1565a1-2803-4780-af05-45dc159c24c5
# ╠═02df5141-e965-4c14-ad00-c433b5f12e30
# ╠═c1ae762d-353a-4143-9116-4c9b7b388766
# ╠═c50322ae-0dd7-4ad2-b308-0242a6cca528
