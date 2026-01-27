import sys

def fix_pdb_split(input_path, output_path, crystal_resname="MOL"):
    print(f"Processing {input_path}...")

    # OpenMM standard residue names for ions
    # If your Force Field uses different names (e.g., CLA for Chloride), change them here.
    ION_MAP = {
        'MG': ('MG', 'MG'), # (Residue Name, Atom Name)
        'CL': ('CL', 'CL'), 
        'NA': ('NA', 'NA'),
        'K':  ('K',  'K')
    }

    current_res_seq = 0
    
    # Solvent trackers
    prev_res_num_infile = None
    prev_chain_infile = None
    solvent_atom_counts = {'H': 0, 'O': 0}

    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                
                # --- PARSE ---
                line = line.rstrip('\n') # Remove newline to handle manually
                
                old_res_name = line[17:20].strip()
                res_seq_infile = line[22:26]
                chain_id = line[21]
                
                # Determine Element
                # 1. Try Col 76-78
                element = line[76:78].strip().upper()
                # 2. Fallback to parsing Atom Name if Element is missing
                if not element:
                    old_atom_name = line[12:16]
                    element = "".join([c for c in old_atom_name if c.isalpha()]).upper()

                # --- LOGIC: CRYSTAL (MOL) -> SPLIT INTO IONS ---
                if old_res_name == crystal_resname:
                    
                    # Every ION gets a NEW Residue Number
                    current_res_seq += 1
                    
                    # Determine standard naming
                    # Default to Element name if not in map
                    if element in ION_MAP:
                        new_res_name, new_atom_name = ION_MAP[element]
                    else:
                        new_res_name = element.ljust(3)
                        new_atom_name = element.ljust(4)

                # --- LOGIC: SOLVENT (HOH) ---
                else:
                    # Increment residue only when the input file indicates a new molecule
                    if res_seq_infile != prev_res_num_infile or chain_id != prev_chain_infile:
                        current_res_seq += 1
                        solvent_atom_counts = {'H': 0, 'O': 0}
                    
                    prev_res_num_infile = res_seq_infile
                    prev_chain_infile = chain_id

                    new_res_name = "HOH" # Standard Water Name
                    
                    # Standard Water Atom Names
                    if element == 'H':
                        solvent_atom_counts['H'] += 1
                        new_atom_name = f"H{solvent_atom_counts['H']}"
                    elif element == 'O':
                        new_atom_name = "O"
                    else:
                        new_atom_name = element

                # --- FORMAT OUTPUT ---
                
                # 1. Atom Name (Cols 13-16)
                if len(new_atom_name) < 4:
                    new_atom_field = f" {new_atom_name:<3}" 
                else:
                    new_atom_field = f"{new_atom_name:<4}"

                # 2. Residue Name (Cols 18-20)
                new_res_name_field = f"{new_res_name:<3}"

                # 3. Residue Sequence (Cols 23-26)
                # Wrap around 10000 if needed (PDB format limit), or let it expand 
                # OpenMM is okay with >9999, but strict PDB requires mod 10000.
                # We will output full integer; OpenMM usually handles it.
                new_res_seq_field = f"{current_res_seq:>4}"
                if len(new_res_seq_field) > 4:
                     # Hex is complex, let's just use modulo if strict PDB is needed, 
                     # but for OpenMM usually 10000+ is fine if columns shift.
                     # Let's stick to 4 chars modulo to keep columns aligned for safety.
                     new_res_seq_field = f"{current_res_seq % 10000:>4}"

                # 4. Construct Line
                # [0:12]  Serial
                # [12:16] Atom Name
                # [16]    AltLoc
                # [17:20] Res Name
                # [20:22] Chain
                # [22:26] Res Seq
                # [26:76] Coords etc.
                
                # Take coords from original line (safe slice)
                # Ensure we don't grab old element column
                line_coords = line[26:76]
                
                new_line_start = (
                    f"{line[:12]}"          # Serial
                    f"{new_atom_field}"     # Atom Name
                    f"{line[16]}"           # AltLoc
                    f"{new_res_name_field}" # Res Name
                    f" {line[21]}"          # Chain (Space at 20, Chain at 21)
                    f"{new_res_seq_field}"  # Res Seq
                    f"{line_coords}"        # Coords
                )

                # 5. Force Element Column (77-78)
                padding = " " * (76 - len(new_line_start))
                fmt_element = f"{element:>2}"
                
                final_line = f"{new_line_start}{padding}{fmt_element}\n"
                f_out.write(final_line)

            else:
                f_out.write(line)

    print(f"Done! Saved to {output_path}")
    print("NOTE: Load this into OpenMM. If you get 'missing parameters', ensure your forcefield.xml includes 'MG' and 'CL'.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fix_pdb_split.py <input.pdb> <output.pdb>")
    else:
        fix_pdb_split(sys.argv[1], sys.argv[2])



# import sys

# def fix_pdb_complete(input_path, output_path, crystal_resname="MOL"):
#     print(f"Processing {input_path}...")

#     # --- TRACKERS ---
#     crystal_counts = {}          # e.g., {'MG': 5, 'CL': 10}
#     solvent_res_seq = 1          # Start at 1, first water becomes 2
    
#     prev_res_num_infile = None
#     prev_chain_infile = None
#     solvent_atom_counts = {'H': 0, 'O': 0}

#     with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
#         for line in f_in:
#             if line.startswith("ATOM") or line.startswith("HETATM"):
                
#                 # --- 1. PARSE INPUT ---
#                 # We strip newlines immediately to handle line length manually later
#                 line = line.rstrip('\n')
                
#                 old_atom_name = line[12:16]
#                 res_name = line[17:20].strip()
#                 chain_id = line[21]
#                 res_seq_infile = line[22:26]
                
#                 # specific logic to find the element. 
#                 # If column 76-78 is empty, we guess from the atom name.
#                 element_in_col = line[76:78].strip() if len(line) >= 78 else ""
                
#                 if element_in_col:
#                     element = element_in_col.upper()
#                 else:
#                     # Guess element from name (e.g. " H1 " -> "H", "Mg" -> "MG")
#                     element = "".join([c for c in old_atom_name if c.isalpha()]).upper()

#                 # --- 2. LOGIC: CRYSTAL (MOL) ---
#                 if res_name == crystal_resname:
#                     new_res_seq = 1
                    
#                     # Track counts for Mg1, Mg2, Cl1...
#                     if element not in crystal_counts:
#                         crystal_counts[element] = 0
#                     crystal_counts[element] += 1
                    
#                     # Format Name
#                     fmt_el = element.capitalize() if len(element) == 2 else element
#                     new_atom_name = f"{fmt_el}{crystal_counts[element]}"

#                 # --- 3. LOGIC: SOLVENT (HOH) ---
#                 else:
#                     # Check for new residue
#                     if res_seq_infile != prev_res_num_infile or chain_id != prev_chain_infile:
#                         solvent_res_seq += 1
#                         solvent_atom_counts = {'H': 0, 'O': 0} # Reset for new water
                        
#                     prev_res_num_infile = res_seq_infile
#                     prev_chain_infile = chain_id
#                     new_res_seq = solvent_res_seq

#                     # Determine Name (H1, H2, O)
#                     if element == 'H':
#                         solvent_atom_counts['H'] += 1
#                         new_atom_name = f"H{solvent_atom_counts['H']}"
#                     elif element == 'O':
#                         new_atom_name = "O"
#                     else:
#                         new_atom_name = element # Fallback

#                 # --- 4. FORMATTING OUTPUT ---
                
#                 # A. Atom Name (Cols 13-16)
#                 if len(new_atom_name) < 4:
#                     new_atom_field = f" {new_atom_name:<3}" 
#                 else:
#                     new_atom_field = f"{new_atom_name:<4}"

#                 # B. Residue Seq (Cols 23-26)
#                 new_res_field = f"{new_res_seq:>4}"

#                 # C. Reconstruct the Start of the Line
#                 # [0:12] Record Type + Serial
#                 # [12:16] New Atom Name
#                 # [16:22] Alt + ResName + Chain
#                 # [22:26] New Res Seq
#                 # [26:54] Coords + Occupancy (standard PDB coords end around 54)
#                 # We just take everything from 26 up to wherever the line ended previously
#                 # BUT we must ensure we don't include the old element column if it existed.
                
#                 # Safely grab coords part (usually up to col 66 for Temp Factor, but let's take up to 76 max)
#                 # PDB width is usually 80. Coords/Occ/Temp end at 66. Segment ID 72-76. Element 77-78.
#                 # We will just take up to column 76 from the original line to preserve data.
#                 line_body = line[26:76] 
                
#                 # Construct new beginning
#                 new_line_start = line[:12] + new_atom_field + line[16:22] + new_res_field + line_body

#                 # D. Force Element Column (Cols 77-78)
#                 # Pad with spaces to reach column 76 (0-indexed 76 is the 77th char)
#                 padding = " " * (76 - len(new_line_start))
                
#                 # Format element right-justified in 2 chars
#                 fmt_element_col = f"{element:>2}"
                
#                 final_line = f"{new_line_start}{padding}{fmt_element_col}\n"
#                 f_out.write(final_line)

#             else:
#                 f_out.write(line)

#     print(f"Done! Saved to {output_path}")

# if __name__ == "__main__":
#     if len(sys.argv) < 3:
#         print("Usage: python fix_pdb_complete.py <input.pdb> <output.pdb>")
#     else:
#         fix_pdb_complete(sys.argv[1], sys.argv[2])

# import sys

# def fix_pdb_enumeration(input_path, output_path, crystal_resname="MOL"):
#     print(f"Processing {input_path}...")

#     # Dictionary to keep track of counts for each element in the crystal
#     # e.g., {'MG': 2, 'CL': 4} -> next Mg will be Mg3
#     element_counts = {}

#     # Start solvent residue numbering at 1, so the first water becomes 2
#     solvent_res_seq = 1
    
#     # Trackers for input solvent residues to detect when to increment
#     prev_res_num_infile = None
#     prev_chain_infile = None

#     with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
#         for line in f_in:
#             if line.startswith("ATOM") or line.startswith("HETATM"):
#                 # PDB Column Parsing
#                 atom_name_old = line[12:16]
#                 res_name = line[17:20].strip()
#                 chain_id = line[21]
#                 res_seq_infile = line[22:26]
                
#                 # Get Element (Columns 77-78), or fallback to parsing Atom Name
#                 element = line[76:78].strip()
#                 if not element:
#                     # Fallback: extract alpha characters from atom name (e.g. "Mg" from "Mg  ")
#                     element = "".join([c for c in atom_name_old if c.isalpha()])

#                 # --- LOGIC FOR CRYSTAL (MOL) ---
#                 if res_name == crystal_resname:
#                     # 1. Force Residue Number to 1
#                     new_res_seq_str = f"{1:>4}"

#                     # 2. Increment Element Counter
#                     # Normalize element to uppercase for dictionary key, but keep casing for output if needed
#                     el_key = element.upper()
#                     if el_key not in element_counts:
#                         element_counts[el_key] = 0
#                     element_counts[el_key] += 1
                    
#                     # 3. Create New Atom Name (e.g., "Mg" + "1" = "Mg1")
#                     # Note: We capitalize the first letter and lower the second for aesthetics (Mg, Cl)
#                     # if the element is 2 letters.
#                     fmt_element = element.capitalize() if len(element) == 2 else element
#                     new_atom_name = f"{fmt_element}{element_counts[el_key]}"

#                     # 4. Format Atom Name Field (Cols 13-16)
#                     # Standard PDB: 4 chars. Align left-ish.
#                     # If name is short (e.g. Mg1 -> 3 chars), standard is " Mg1"
#                     # If name is long (e.g. Mg100 -> 5 chars), we risk overflow but OpenMM usually handles it.
#                     if len(new_atom_name) < 4:
#                         # Pad with leading space for alignment (standard PDB style)
#                         new_atom_name_field = f" {new_atom_name:<3}" 
#                     elif len(new_atom_name) == 4:
#                         new_atom_name_field = f"{new_atom_name:<4}"
#                     else:
#                         # If it exceeds 4 chars (Mg1000), we just print it. 
#                         # This breaks strict PDB column width but works in most parsers.
#                         new_atom_name_field = new_atom_name[:4] # Clip or handle overflow? 
#                         # Better to overflow slightly than clip unique numbers:
#                         new_atom_name_field = new_atom_name
                        
#                     # Reconstruct line. Note: if new_atom_name_field > 4 chars, 
#                     # we must be careful not to eat into Residue Name.
#                     # We assume standard parsing separates by space if columns merge.
                    
#                     # Safe reconstruction for standard width:
#                     if len(new_atom_name_field) <= 4:
#                          line = line[:12] + new_atom_name_field + line[16:22] + new_res_seq_str + line[26:]
#                     else:
#                         # Handle Overflow: Manually construct spaced string
#                         # This effectively shifts the rest of the line right, which OpenMM usually accepts
#                         part1 = line[:12]
#                         part2 = new_atom_name_field # The long name
#                         part3 = f" {res_name} {chain_id}{new_res_seq_str}" # Reconstruction of middle
#                         part4 = line[26:] # Coords
#                         line = f"{part1}{part2}{part3}{part4}"


#                 # --- LOGIC FOR SOLVENT (HOH) ---
#                 else:
#                     # Detect if this is a new residue in the input file
#                     if res_seq_infile != prev_res_num_infile or chain_id != prev_chain_infile:
#                         solvent_res_seq += 1
                    
#                     # Update trackers
#                     prev_res_num_infile = res_seq_infile
#                     prev_chain_infile = chain_id
                    
#                     # Create Residue String
#                     new_res_seq_str = f"{solvent_res_seq:>4}"
                    
#                     # Keep original atom name (H1, H2, O), just update residue number
#                     line = line[:22] + new_res_seq_str + line[26:]

#                 f_out.write(line)

#             else:
#                 # Header/Connect/Footer
#                 f_out.write(line)

#     print(f"Done! Fixed file saved to: {output_path}")

# if __name__ == "__main__":
#     if len(sys.argv) < 3:
#         print("Usage: python fix_pdb_enum.py <input.pdb> <output.pdb>")
#     else:
#         fix_pdb_enumeration(sys.argv[1], sys.argv[2])




# import sys

# def fix_pdb_crystal(input_path, output_path, crystal_resname="MOL"):
#     """
#     Splits a 'clumped' crystal residue into individual residues for each ion
#     and renumbers subsequent solvent residues sequentially.
#     """
    
#     # We start residue numbering at 0 so the first one becomes 1
#     current_res_seq = 0
    
#     # Trackers to know when we hit a new solvent residue
#     prev_res_num_infile = None
#     prev_chain_infile = None

#     print(f"Processing {input_path}...")
    
#     with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
#         for line in f_in:
#             if line.startswith("ATOM") or line.startswith("HETATM"):
#                 # PDB Column definitions:
#                 # 17-20: Residue Name (e.g., "MOL", "HOH")
#                 # 21: Chain ID
#                 # 22-26: Residue Sequence Number
                
#                 res_name = line[17:20].strip()
#                 chain_id = line[21]
#                 res_seq_infile = line[22:26]

#                 # LOGIC:
#                 # 1. If it is the Crystal (MOL), EVERY atom is a new residue.
#                 if res_name == crystal_resname:
#                     current_res_seq += 1
                
#                 # 2. If it is solvent (HOH), we only increment residue ID
#                 # if the input file's residue ID (or chain) changes.
#                 else:
#                     if res_seq_infile != prev_res_num_infile or chain_id != prev_chain_infile:
#                         current_res_seq += 1

#                 # Update trackers
#                 prev_res_num_infile = res_seq_infile
#                 prev_chain_infile = chain_id

#                 # FORMATTING:
#                 # Create the new 4-character residue string (right justified)
#                 # We cap at 9999 to adhere to standard PDB format, though OpenMM 
#                 # can often handle overflows or hybrid36 if strictly necessary.
#                 new_res_str = f"{current_res_seq:>4}"
#                 if len(new_res_str) > 4:
#                     # Fallback for >9999 residues: hex or mod 10000 
#                     # (Simple mod fix to keep column width valid)
#                     new_res_str = f"{current_res_seq % 10000:>4}"

#                 # Reconstruct the line preserving fixed width
#                 # 0-22: Start of line up to residue number
#                 # 22-26: New Residue Number
#                 # 26+:  Rest of the line (coords, temp factor, etc.)
#                 new_line = line[:22] + new_res_str + line[26:]
#                 f_out.write(new_line)

#             else:
#                 # Copy header/footer/CONECT lines exactly as is
#                 f_out.write(line)

#     print(f"Done! Fixed file saved to: {output_path}")

# if __name__ == "__main__":
#     # Usage: python fix_pdb.py input.pdb fixed.pdb
#     if len(sys.argv) < 3:
#         print("Usage: python fix_pdb.py <input_pdb> <output_pdb>")
#     else:
#         fix_pdb_crystal(sys.argv[1], sys.argv[2])











# import sys

# # --- CONFIGURATION ---
# input_filename = "input-liquid.pdb"  # <--- REPLACE THIS with your actual filename
# output_filename = "input-liquid-fixed.pdb"
# # ---------------------

# print(f"Reading {input_filename}...")

# with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
#     residue_counter = 0
#     atom_counter = 0
    
#     for line in infile:
#         if line.startswith("ATOM") or line.startswith("HETATM"):
#             # Identify the element (Column 77-78 usually, or parsed from name)
#             # In your file, the name is in cols 12-16.
#             atom_name = line[12:16].strip()
            
#             # --- LOGIC: Group every 3 atoms into a new residue ---
#             # Your file order is H, H, O (3 atoms per water)
#             if atom_counter % 3 == 0:
#                 residue_counter += 1
            
#             atom_counter += 1

#             # --- RENAME ATOMS (Standardize for AMBER/CHARMM) ---
#             # Your file has H, H, O. 
#             # Standard templates often prefer specific names like H1, H2, O.
#             if atom_name == 'O':
#                 new_name = " O  " # Standard padding
#             elif atom_name == 'H':
#                 # Determine if it's the first or second H in the trio
#                 if atom_counter % 3 == 1:
#                     new_name = " H1 "
#                 else:
#                     new_name = " H2 "
#             else:
#                 new_name = line[12:16] # Keep original if unknown

#             # --- REWRITE THE LINE ---
#             # Columns (0-indexed logic):
#             # 17-20: Residue Name -> Change "MOL " to "HOH "
#             # 22-26: Residue Seq  -> Change "1" to unique number
            
#             # Slice the original line parts
#             part1 = line[:12]   # ATOM + Serial
#             part2 = line[16:17] # Alternate loc
#             part3 = line[21:22] # Chain ID
#             coords = line[26:]  # XYZ and the rest

#             # Format new columns
#             new_res_name = "HOH "
#             new_res_seq = f"{residue_counter:4d}"
            
#             # Assemble
#             new_line = f"{part1}{new_name}{part2}{new_res_name}{part3}{new_res_seq}{coords}"
#             outfile.write(new_line)
            
#         else:
#             # Copy header/footer lines (CRYST1, END) exactly as is
#             outfile.write(line)

# print(f"Success! Fixed file written to: {output_filename}")
# print(f"Found {residue_counter} water molecules.")