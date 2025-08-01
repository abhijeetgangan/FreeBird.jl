"""
    EnergyEval
    
Module for evaluating energy-related quantities for a system.
"""
module EnergyEval

using AtomsBase
using Unitful
using StaticArrays
using PythonCall
using ..AbstractWalkers
using ..AbstractPotentials
using ..AbstractHamiltonians

export pbc_dist
export interacting_energy, frozen_energy
export single_site_energy


include("atomistic_energies.jl")

include("lattice_energies.jl")

end # module EnergyEval