import math
import sys

# --- CONFIGURATION ---
INPUT_FILE = "start_Mg_Cl.pdb"    # <--- REPLACE THIS
OUTPUT_FILE = "start_Mg_Cl_fixed.pdb"
# ---------------------

def get_pbc_dist_sq(atom1, atom2, box_dims):
    """Calculates squared distance considering Periodic Boundary Conditions (PBC)"""
    dx = atom1['x'] - atom2['x']
    dy = atom1['y'] - atom2['y']
    dz = atom1['z'] - atom2['z']
    
    # Apply Minimum Image Convention (wrap deltas to be within -L/2 to L/2)
    if box_dims:
        Lx, Ly, Lz = box_dims
        dx -= Lx * round(dx / Lx)
        dy -= Ly * round(dy / Ly)
        dz -= Lz * round(dz / Lz)
        
    return dx*dx + dy*dy + dz*dz

print(f"Reading {INPUT_FILE}...")

atoms = {'O': [], 'H': [], 'Ion': []}
box_dims = None # [Lx, Ly, Lz]

# 1. READ AND PARSE
with open(INPUT_FILE, 'r') as f:
    for line in f:
        if line.startswith("CRYST1"):
            # Parse box dimensions: CRYST1   20.004   20.004 ...
            parts = line.split()
            box_dims = [float(parts[1]), float(parts[2]), float(parts[3])]
            print(f"Detected Box Size: {box_dims}")

        if line.startswith("ATOM") or line.startswith("HETATM"):
            elem = line[76:78].strip().upper()
            name = line[12:16].strip().upper()
            
            # Extract coordinates
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            atom_data = {
                'line_orig': line,
                'x': x, 'y': y, 'z': z,
                'id': line[6:11].strip()
            }
            
            if 'MG' in name or 'CL' in name or 'NA' in name:
                atoms['Ion'].append(atom_data)
            elif elem == 'O' or 'O' in name:
                atoms['O'].append(atom_data)
            elif elem == 'H' or 'H' in name:
                atoms['H'].append(atom_data)

print(f"Found {len(atoms['O'])} Oxygens, {len(atoms['H'])} Hydrogens, {len(atoms['Ion'])} Ions.")

if box_dims is None:
    print("WARNING: No CRYST1 line found. Assuming non-periodic (this may fail for your file).")

# 2. RECONNECT WATER MOLECULES
with open(OUTPUT_FILE, 'w') as out:
    # Write Header
    if box_dims:
        out.write(f"CRYST1{box_dims[0]:9.3f}{box_dims[1]:9.3f}{box_dims[2]:9.3f}  90.00  90.00  90.00 P 1\n")
    
    residue_count = 1
    used_hydrogens = set()
    
    # Loop through every Oxygen and find its 2 closest Hydrogens
    for oxy in atoms['O']:
        # Sort ALL hydrogens by distance to this Oxygen (computationally heavy but safe)
        # We only check available (unused) hydrogens
        candidates = []
        for i, hyd in enumerate(atoms['H']):
            if i in used_hydrogens: continue
            d2 = get_pbc_dist_sq(oxy, hyd, box_dims)
            if d2 < 9.0: # Optimization: Only care if within 3 Angstroms (3^2=9)
                candidates.append((d2, i, hyd))
        
        # Pick top 2 closest
        candidates.sort(key=lambda x: x[0])
        best_2 = candidates[:2]
        
        if len(best_2) < 2:
            print(f"Warning: Oxygen {oxy['id']} only found {len(best_2)} hydrogens! (Is the file truncated?)")
        
        # Form the Residue
        # 1. Write Oxygen
        out.write(f"ATOM  {1:5d}  O   HOH {residue_count:4d}    {oxy['x']:8.3f}{oxy['y']:8.3f}{oxy['z']:8.3f}  1.00  0.00           O  \n")
        
        # 2. Write Hydrogens
        h_counter = 1
        for _, idx, hyd in best_2:
            used_hydrogens.add(idx)
            # Wrap Hydrogen relative to Oxygen to fix visual bonds
            # (Optional: If H is 200 and O is 5, move H to 5 +/- bond_length)
            # For strict PDB compliance we just write coords, OpenMM handles wrapping if topology is good.
            out.write(f"ATOM  {1:5d}  H{h_counter}  HOH {residue_count:4d}    {hyd['x']:8.3f}{hyd['y']:8.3f}{hyd['z']:8.3f}  1.00  0.00           H  \n")
            h_counter += 1
            
        residue_count += 1

    # 3. WRITE IONS
    for ion in atoms['Ion']:
        # Parse original name/element to keep it correct (MG vs CL)
        orig_name = ion['line_orig'][12:16].strip()
        orig_elem = ion['line_orig'][76:78].strip()
        # Ensure residue name is unique/standard (e.g. MG, CL)
        res_name = orig_name
        
        out.write(f"ATOM  {1:5d} {orig_name:>4} {res_name:>3} {residue_count:4d}    {ion['x']:8.3f}{ion['y']:8.3f}{ion['z']:8.3f}  1.00  0.00          {orig_elem:>2}  \n")
        residue_count += 1

print(f"Done! Saved to {OUTPUT_FILE}")