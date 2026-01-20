from ase.io import read, write

# Read the XYZ file
atoms = read('start_Mg_Cl.data', format='lammps-data')

# Write to a PDB file
write('start_Mg_Cl.pdb', atoms, format='proteindatabank')