begin
	function getTension(position, t, rStim, D)
		if D != 0.0
			rChannel = sqrt(position[1]^2 + position[2]^2) - rStim
			return exp(-rChannel^2/(4*D*t))/(4π * D * t)
		else 
			return 0.0
		end
	end

	function tensionMap(size, t, rStim, mStim ; D=24000.0, clamp=true)
	
		membrane = fill(0.0, (size, size))

		center= [size/2, size/2]
		@inbounds for I ∈ CartesianIndices(membrane)
			membrane[I] = getTension(Tuple(I) .- center, t, rStim, D)
		end

        if clamp == true
			mask = BitArray(sqrt(x^2+y^2) <= rStim for x = 1-0.5*size:0.5*size, y = 1-0.5*size:0.5*size) ##stimulus mask
			
            membrane .*= mStim/maximum(membrane)
            membrane[mask] .= mStim
        end

        return membrane
	end

	function seedChannels(arr; occupancy=0.01, seed=62)
		Random.seed!(seed)
		d = Binomial(1, occupancy)

		@inbounds for I in CartesianIndices(arr)
			arr[I] = rand(d, 1)[1]
		end
		
		return arr
	end
end