using PythonCall
using AtomsBase, StaticArrays, Unitful

export ASECalculator, ASELennardJones

"""
    struct ASECalculator

Light-weight wrapper around an `ase` calculator (accessed via PythonCall) so it
can be passed around in Julia code just like the built-in `LJParameters`.

The wrapped calculator *must* be able to return energies in *electron-volts* via
`atoms.get_potential_energy()`.
"""
struct ASECalculator <: LennardJonesParameterSets
    calc::Py
end

"""
    ASELennardJones(;epsilon=1.0, sigma=1.0, cutoff=Inf)

Convenience constructor that creates an `ase.calculators.lj.LennardJones`
calculator with the supplied parameters and wraps it in an `ASECalculator`.
All energy/length units follow the ASE convention (eV / Å).
"""
function ASELennardJones(; epsilon = 1.0, sigma = 1.0, cutoff = Inf)
    lj_mod = pyimport("ase.calculators.lj")
    # ASE uses the keyword `rc` for the cut-off radius
    rc = isfinite(cutoff) ? cutoff * sigma : cutoff
    py_calc = lj_mod.LennardJones(; epsilon = epsilon, sigma = sigma, rc = rc)
    return ASECalculator(py_calc)
end

# Convert an AtomsBase `AbstractSystem` to an `ase.Atoms` object.
function _to_ase(system::AtomsBase.AbstractSystem)
    ase = pyimport("ase")
    n = length(system)
    symbols = [string(atomic_symbol(system, i)) for i in 1:n]
    # collect positions (convert to Å and strip units)
    pos = zeros(Float64, n, 3)
    for i in 1:n
        p = position(system, i)
        pos[i, 1] = ustrip(uconvert(u"Å", p[1]))
        pos[i, 2] = ustrip(uconvert(u"Å", p[2]))
        pos[i, 3] = ustrip(uconvert(u"Å", p[3]))
    end

    py_atoms = ase.Atoms(symbols, positions = pos)

    # set cell and pbc if available
    try
        cell_vecs = cell_vectors(system)
        if !isnothing(cell_vecs)
            cell_mat = zeros(Float64, 3, 3)
            for i in 1:3
                for j in 1:3
                    cell_mat[i, j] = ustrip(uconvert(u"Å", cell_vecs[i][j]))
                end
            end
            py_atoms.cell = cell_mat
        end
    catch err
        # Fallback: ignore if not orthorhombic / not available
    end
    try
        py_atoms.pbc = tuple(periodicity(system)...)
    catch
        # ignore if periodicity not defined
    end

    return py_atoms
end

# Dummy frozen_energy for compatibility with LJAtomWalkers constructor when no
# atoms are actually frozen (returns zero).
function frozen_energy(system::AtomsBase.AbstractSystem,
                       calc::ASECalculator,
                       list_num_par::Vector{Int},
                       frozen::Vector{Bool})
    # For now we simply return 0 because ASE calculators do not currently
    # support separating frozen–frozen interactions.
    return 0.0u"eV"
end 