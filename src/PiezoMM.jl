module PiezoMM

using LinearAlgebra, Distributions, Random

include("markovModel.jl")
include("tensionDiffusion.jl")

export generateTMatrix, markovSequence, simChannels, getTension, simulateTension, equilibriumState, tensionMap, seedChannels

end
