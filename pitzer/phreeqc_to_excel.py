import re
import pandas as pd
import os

# --- CONFIGURATION ---
file_path = 'pitzerCopy.dat'
output_filename = 'pitzer_parameters_charges.xlsx'

def parse_pitzer_dat(file_content):
    """
    Parses a PHREEQC Pitzer database file to extract species charges and Pitzer coefficients.
    Robustly handles block switching and species regex.
    """
    
    BLOCK_HEADERS = {
        'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES', 'PHASES', 'PITZER', 
        'EXCHANGE_MASTER_SPECIES', 'EXCHANGE_SPECIES', 
        'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES', 
        'MEAN_GAMMAS', 'RATES', 'KINETICS', 'END'
    }

    # Regex for Pitzer sub-tags (e.g., -B0, -THETA)
    re_sub_tag = re.compile(r'^\s*-(B0|B1|B2|C0|THETA|LAMBDA|ZETA|PSI|MU|ETA|ALPHAS)', re.IGNORECASE)
    
    # Regex to identify species tokens in a reaction line
    # Looks for alphanumeric strings ending in + or - (optionally with a number)
    # Examples: Na+, Cl-, SO4-2, Ca+2, H+
    re_charged_species = re.compile(r'([A-Za-z0-9\(\)]+[\+\-]\d?)')

    data_pitzer = []
    data_species = {}
    
    current_block = None
    current_tag = None
    
    lines = file_content.split('\n')
    
    for line in lines:
        clean_line = line.split('#')[0].strip()
        if not clean_line:
            continue
            
        # --- BLOCK SWITCHING ---
        first_word = clean_line.split()[0].upper()
        if first_word in BLOCK_HEADERS:
            current_block = first_word
            current_tag = None
            continue

        # --- PARSE SPECIES CHARGE (SOLUTION_SPECIES) ---
        if current_block == 'SOLUTION_SPECIES':
            # Skip lines that are just properties (start with -) or keywords
            if clean_line.startswith('-') or '=' not in clean_line:
                continue

            # Find all tokens that look like charged species
            matches = re_charged_species.findall(clean_line)
            for token in matches:
                # Determine charge
                charge = 0
                if token.endswith('+'):
                    charge = 1
                elif token.endswith('-'):
                    charge = -1
                elif re.search(r'\+\d+$', token): # e.g., +2
                    charge = int(token.split('+')[-1])
                elif re.search(r'-\d+$', token): # e.g., -2
                    charge = -int(token.split('-')[-1])
                
                # Sanity check: Ensure it looks like a chemical formula (not just a number)
                if charge != 0 and token[0].isalpha():
                    data_species[token] = charge

        # --- PARSE PITZER COEFFICIENTS ---
        if current_block == 'PITZER':
            tag_match = re_sub_tag.match(clean_line)
            if tag_match:
                current_tag = tag_match.group(1).upper()
                continue
            
            if current_tag:
                parts = clean_line.split()
                # Data lines start with a Species name (Letter), not number or parenthesis
                if not parts[0][0].isalpha() and not parts[0].startswith('('):
                    continue

                num_species = 3 if current_tag in ['ZETA', 'PSI', 'MU'] else 2
                
                specs = parts[:num_species]
                coeffs_raw = parts[num_species:]
                
                coeffs = []
                for c in coeffs_raw:
                    try:
                        coeffs.append(float(c))
                    except ValueError:
                        break 
                
                while len(coeffs) < 6: coeffs.append(0.0)
                
                entry = {
                    'Parameter': current_tag,
                    'Species 1': specs[0] if len(specs) > 0 else '',
                    'Species 2': specs[1] if len(specs) > 1 else '',
                    'Species 3': specs[2] if len(specs) > 2 else '',
                    'A0 (25C)': coeffs[0],
                    'A1 (T)': coeffs[1],
                    'A2 (lnT)': coeffs[2],
                    'A3 (T-Tr)': coeffs[3],
                    'A4 (T^2)': coeffs[4],
                    'A5 (1/T)': coeffs[5]
                }
                data_pitzer.append(entry)

    return data_pitzer, data_species

# --- EXECUTION ---
if not os.path.exists(file_path):
    print(f"Error: {file_path} not found.")
else:
    print("Parsing file...")
    pitzer_data, species_map = parse_pitzer_dat(open(file_path).read())

    df_pitzer = pd.DataFrame(pitzer_data)
    df_species = pd.DataFrame(list(species_map.items()), columns=['Species', 'Charge'])
    
    # Ensure correct columns exist
    cols = ['Parameter', 'Species 1', 'Species 2', 'Species 3', 
            'A0 (25C)', 'A1 (T)', 'A2 (lnT)', 'A3 (T-Tr)', 'A4 (T^2)', 'A5 (1/T)']
    for c in cols:
        if c not in df_pitzer.columns: df_pitzer[c] = 0.0
    
    with pd.ExcelWriter(output_filename) as writer:
        df_pitzer[cols].to_excel(writer, sheet_name='Coefficients', index=False)
        df_species.to_excel(writer, sheet_name='Species_Charges', index=False)

    print(f"Success! Captured {len(df_species)} species charges (check for Na+: {'Na+' in species_map}).")
    print(f"File saved: {output_filename}")