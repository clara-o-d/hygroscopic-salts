from ase.io import read, write
from ase.build import make_supercell

atoms = read("MgCl2.cif")
atoms = atoms.repeat((4,4,2))   # small crystal slab
write("MgCl2_crystal.pdb", atoms)