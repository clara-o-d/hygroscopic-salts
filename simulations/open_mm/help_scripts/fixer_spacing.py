import re

# Update these filenames as needed
input_file = 'mixture_fixed.pdb'
output_file = 'mixture_ready.pdb'

def brute_force_fix(infile, outfile):
    print(f"Brute-force reformatting {infile}...")
    
    try:
        with open(infile, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: {infile} not found.")
        return

    fixed_lines = []
    atom_count = 0
    
    # Standard PDB format string for ATOM records
    # Forces exact column alignment:
    # ATOM (1-4) | Serial (7-11) | Name (13-16) | Res (18-20) | ResNum (23-26) | X,Y,Z (31-54)
    pdb_format = "ATOM  {:5d} {:^4s} {:3s} A{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}\n"

    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('CRYST1'):
            fixed_lines.append(line + "\n")
            continue
            
        parts = line.split()
        # Check if line starts with ATOM or HETATM and has enough data
        if len(parts) >= 7 and parts[0] in ['ATOM', 'HETATM']:
            atom_count += 1
            try:
                # Basic mapping based on your snippet: 
                # [0]ATOM [1]1 [2]H1 [3]HOH [4]1 [5]X [6]Y [7]Z
                name = parts[2]
                res_name = parts[3]
                res_num = int(parts[4])
                x = float(parts[5])
                y = float(parts[6])
                z = float(parts[7])
                # Guess element from name
                element = "".join([c for c in name if c.isalpha()])[:2]

                new_line = pdb_format.format(
                    atom_count, name[:4], res_name[:3], res_num, x, y, z, element
                )
                fixed_lines.append(new_line)
            except (ValueError, IndexError):
                continue

    with open(outfile, 'w') as f:
        f.writelines(fixed_lines)
        f.write("END\n")

    print(f"Done! Created {outfile} with {atom_count} atoms.")

if __name__ == "__main__":
    brute_force_fix(input_file, output_file)