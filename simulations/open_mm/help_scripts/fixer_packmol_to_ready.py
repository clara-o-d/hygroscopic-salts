import sys
import re

# --- CONFIGURATION ---
INPUT_FILE = 'mixture_raw.pdb'   # Your Packmol output file
OUTPUT_FILE = 'mixture_ready.pdb'
BOX_SIZE = 40.0                      # Simulation box size in Angstroms
# ---------------------

def fix_packmol_pdb(infile, outfile, box_size):
    print(f"Reading {infile}...")
    
    fixed_lines = []
    
    # 1. Add CRYST1 Header
    cryst_line = f"CRYST1{box_size:9.3f}{box_size:9.3f}{box_size:9.3f}  90.00  90.00  90.00 P 1           1\n"
    fixed_lines.append(cryst_line)

    atom_serial = 0
    water_h_counts = {} # To track H1 vs H2

    try:
        with open(infile, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find {infile}")
        return

    for line in lines:
        parts = line.split()
        
        # Filter for atom lines only
        if not parts or parts[0] not in ['ATOM', 'HETATM']:
            continue
            
        try:
            # --- ROBUST COORDINATE DETECTION ---
            # Find the first column that contains a '.' AND looks like a number.
            # This prevents integers (like serial numbers) from confusing the parser.
            coord_start_idx = -1
            for i, part in enumerate(parts):
                # Check for decimal point to distinguish coords from IDs
                if '.' in part: 
                    # Verify it is actually a number (handles negative signs too)
                    clean_part = part.replace('.', '', 1).lstrip('-')
                    if clean_part.isdigit():
                        # Ensure we have X, Y, Z available
                        if i + 2 < len(parts):
                            coord_start_idx = i
                            break
            
            if coord_start_idx == -1:
                print(f"Skipping line (no coords found): {line.strip()}")
                continue

            # Parse Coords
            x = float(parts[coord_start_idx])
            y = float(parts[coord_start_idx+1])
            z = float(parts[coord_start_idx+2])
            
            # --- IDENTIFY METADATA ---
            # Packmol columns are usually: [Record, Serial, Name, ResName, (Chain), ResNum, X...]
            # We work backwards from X to find ResNum, ResName, Name.
            
            # 1. Residue Number (Usually immediately before X)
            raw_res_num = parts[coord_start_idx - 1]
            # Extract only digits (handle cases like 'A74' or '1')
            res_num_digits = "".join(filter(str.isdigit, raw_res_num))
            
            if res_num_digits:
                res_num = int(res_num_digits)
            else:
                # If immediate neighbor isn't the ID, try one further back (skip Chain)
                raw_res_num = parts[coord_start_idx - 2]
                res_num_digits = "".join(filter(str.isdigit, raw_res_num))
                res_num = int(res_num_digits) if res_num_digits else 1

            # 2. Residue Name & Atom Name
            # Typically: Name is index 2, ResName is index 3
            original_atom_name = parts[2]
            res_name = parts[3]
            
            # --- FIX ATOM NAMING ---
            
            # Fix Water (H -> H1/H2)
            final_atom_name = original_atom_name
            if res_name in ["HOH", "SOL", "WAT"]:
                if original_atom_name.startswith("H"):
                    count = water_h_counts.get(res_num, 0) + 1
                    water_h_counts[res_num] = count
                    final_atom_name = f"H{count}"
            
            # Fix Chlorine (CLA -> CL)
            if "CL" in res_name.upper():
                res_name = "Cl"
                final_atom_name = "Cl"

            # Fix Magnesium
            if "MG" in res_name.upper():
                res_name = "Mg"
                final_atom_name = "Mg"

            # Determine Element (Strip numbers)
            element = "".join([i for i in final_atom_name if not i.isdigit()])[:2]
            
            atom_serial += 1
            
            # --- WRITE STANDARD PDB LINE ---
            # Handles correct spacing for atom names (4 chars vs <4 chars)
            if len(final_atom_name) < 4:
                name_fmt = f" {final_atom_name:<3s}" # " O  "
            else:
                name_fmt = f"{final_atom_name:<4s}"  # "H1  "

            new_line = (
                f"ATOM  {atom_serial:5d} {name_fmt} {res_name:<3s} A{res_num:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
            )
            fixed_lines.append(new_line)

        except Exception as e:
            print(f"Skipping weird line: {line.strip()} | Error: {e}")
            continue

    with open(outfile, 'w') as f:
        f.writelines(fixed_lines)
        f.write("END\n")
    
    print(f"Success! Processed {atom_serial} atoms.")
    print(f"Saved to: {outfile}")

if __name__ == "__main__":
    fix_packmol_pdb(INPUT_FILE, OUTPUT_FILE, BOX_SIZE)