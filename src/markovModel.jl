function generateTMatrix(σ, σ₅₀, slope, τ)
	## Publication rates in s⁻¹
	k₁ = 5.1 * exp(σ/slope)
	k₋₁ = 5 * exp(σ₅₀/slope)
	k₂ = 8.0 
	k₋₂ = 0.4 
	k₃ = 34.6 * exp(-σ/slope)
	k₋₃ = (k₁*k₂*k₃)/(k₋₁*k₋₂)
	k₄ = 4.0 
	k₋₄ = 0.6 

	## Convert to probabilities
	P₀₁ =  k₁ * τ
	P₀₂ = k₋₃ * τ
	P₀₀ = 1 - (P₀₁ + P₀₂)
	
	P₁₀ = k₋₁ * τ
	P₁₂ = k₂ * τ
	P₁₃ = k₄ * τ
	P₁₁ = 1 - (P₁₀ + P₁₂ + P₁₃)
	
	P₂₀ = k₃ * τ 
	P₂₁ = k₋₂ * τ
	P₂₂ = 1 - (P₂₀ + P₂₁)
	 
		
	P₃₁ = k₋₄ * τ
	P₃₃ = 1 - P₃₁
		
	pTransition = [
		P₀₀ P₀₁ P₀₂ 0;
		P₁₀ P₁₁ P₁₂ P₁₃;
		P₂₀ P₂₁ P₂₂ 0;
		0 P₃₁ 0 P₃₃
	]

	return pTransition
end

## Markov simulation
function markovSequence(pInit, nSamples; σ=:none, σ₅₀=2.7, slope=0.8, τ=1e-5, seed=62)

	σ == :none && (σ = zeros(nSamples))

	pTransition = generateTMatrix(σ[1], σ₅₀, slope, τ)
		
	dInit = Multinomial(1, pInit[1,:])
	states = zeros(Int64, nSamples)
	states[1] = findfirst(x -> x .== 1, rand(dInit, 1))[1][1]
			
	@inbounds for i in 2:nSamples
		σ[i] != σ[i-1] && (pTransition = generateTMatrix(σ[i], σ₅₀, slope, τ))
			
		check = sum(pTransition, dims=2)
		@assert sum(check .≈ 1.0) == 4 "$check Each row must sum to 1.0 for a probability matrix"

		pT = pTransition[states[i-1], :]
		d = Multinomial(1, pT)
		states[i] = findfirst(x -> x .== 1, rand(d, 1))[1][1]
	end
	return states
end

	
## Simulate multiple channels
	
function simChannels(channels, size, pInit, tArr, rStim, mStim, τ;
	openState=2, σ₅₀=2.7, slope=0.8, D=240, seed=62)

	tsteps = vcat(zeros(Int64(round(0.05/τ))), collect(tArr), zeros(Int64(round(0.05/τ))))
	stim = fill(0.0, (length(channels), length(tsteps)))

	peakIdx = argmax(tensionMap(size, τ, rStim, mStim; D=D, clamp=false))

	@inbounds Threads.@threads for t ∈ 1:length(tsteps)
		if tsteps[t] == 0.0
			stim[:, t] = zeros(length(channels))
		else
			peak = getTension(Tuple(peakIdx) .- size/2, tsteps[t], rStim, D)
			
			@inbounds for c ∈ 1:length(channels)
				if sqrt(channels[c][1]^2 + channels[c][2]^2) <= rStim
					stim[c,t] = mStim
				else
					temp = getTension(channels[c], tsteps[t], rStim, D)/peak
					stim[c,t] = temp * mStim
				end
			end
		end
	end

	sim = fill(0.0, (length(channels), length(tsteps)))
	@inbounds Threads.@threads for c ∈ 1:length(channels)
		sim[c, :] = markovSequence(pInit, length(tsteps); σ=stim[c, :], σ₅₀=σ₅₀, slope=slope, τ=τ, seed=seed)
	end
		
	return sim
end

function equilibriumState(T)
	nStates = size(T, 1)
	a = [T' - I ; ones(nStates)']
	b = [zeros(4, 1); 1]
	return (a' * a)\(a' * b)
end
