"""
    AbstractPotentials

Module for defining and implementing potentials.
"""
module AbstractPotentials

using Unitful

export AbstractPotential
export LJParameters, lj_energy
export CompositeLJParameters
export LennardJonesParameterSets
export ASECalculator, ASELennardJones

abstract type AbstractPotential end

include("lennardjones.jl")
include("ase_calculators.jl")

end # module Potentials