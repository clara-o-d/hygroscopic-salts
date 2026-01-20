import sys

# --- CONFIGURATION ---
input_filename = "input-liquid.pdb"  # <--- REPLACE THIS with your actual filename
output_filename = "input-liquid-fixed.pdb"
# ---------------------

print(f"Reading {input_filename}...")

with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
    residue_counter = 0
    atom_counter = 0
    
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Identify the element (Column 77-78 usually, or parsed from name)
            # In your file, the name is in cols 12-16.
            atom_name = line[12:16].strip()
            
            # --- LOGIC: Group every 3 atoms into a new residue ---
            # Your file order is H, H, O (3 atoms per water)
            if atom_counter % 3 == 0:
                residue_counter += 1
            
            atom_counter += 1

            # --- RENAME ATOMS (Standardize for AMBER/CHARMM) ---
            # Your file has H, H, O. 
            # Standard templates often prefer specific names like H1, H2, O.
            if atom_name == 'O':
                new_name = " O  " # Standard padding
            elif atom_name == 'H':
                # Determine if it's the first or second H in the trio
                if atom_counter % 3 == 1:
                    new_name = " H1 "
                else:
                    new_name = " H2 "
            else:
                new_name = line[12:16] # Keep original if unknown

            # --- REWRITE THE LINE ---
            # Columns (0-indexed logic):
            # 17-20: Residue Name -> Change "MOL " to "HOH "
            # 22-26: Residue Seq  -> Change "1" to unique number
            
            # Slice the original line parts
            part1 = line[:12]   # ATOM + Serial
            part2 = line[16:17] # Alternate loc
            part3 = line[21:22] # Chain ID
            coords = line[26:]  # XYZ and the rest

            # Format new columns
            new_res_name = "HOH "
            new_res_seq = f"{residue_counter:4d}"
            
            # Assemble
            new_line = f"{part1}{new_name}{part2}{new_res_name}{part3}{new_res_seq}{coords}"
            outfile.write(new_line)
            
        else:
            # Copy header/footer lines (CRYST1, END) exactly as is
            outfile.write(line)

print(f"Success! Fixed file written to: {output_filename}")
print(f"Found {residue_counter} water molecules.")