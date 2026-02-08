#!/usr/bin/env python3
"""
Parse thermo_ymp.R2.tdat thermodynamic database and export to CSV files.

This script extracts:
- Elements
- Basis species
- Redox couples
- Aqueous species
- Minerals
- Gases
- Pitzer parameters (binary, theta, lambda, psi)
"""

import csv
import re
from pathlib import Path


def parse_elements(lines, start_idx):
    """Parse elements section."""
    elements = []
    i = start_idx
    
    # Find number of elements
    while i < len(lines):
        match = re.match(r'\s*(\d+)\s+elements', lines[i])
        if match:
            num_elements = int(match.group(1))
            i += 2  # Skip blank line
            break
        i += 1
    
    # Parse element entries
    count = 0
    while i < len(lines) and count < num_elements:
        line = lines[i].strip()
        if line and not line.startswith('*') and not line.startswith('-end-'):
            # Format: "Element (Symbol) mole wt.= value g"
            match = re.match(r'(\w+)\s+\((\w+)\)\s+mole wt\.=\s+([\d.]+)\s+g', line)
            if match:
                elements.append({
                    'name': match.group(1),
                    'symbol': match.group(2),
                    'mole_weight': float(match.group(3))
                })
                count += 1
        i += 1
    
    return elements


def parse_species_section(lines, start_idx, section_name):
    """Parse basis species, aqueous species, minerals, or gases."""
    species_list = []
    i = start_idx
    
    # Find number of species
    while i < len(lines):
        match = re.match(rf'\s*(\d+)\s+{section_name}', lines[i])
        if match:
            num_species = int(match.group(1))
            i += 2
            break
        i += 1
    
    # Parse species entries
    current_species = None
    while i < len(lines):
        line = lines[i]
        
        if line.strip() == '-end-':
            if current_species:
                species_list.append(current_species)
            break
        
        # New species name (not indented, not blank, not comment)
        if line and not line.startswith(' ') and not line.startswith('*') and line.strip():
            if current_species:
                species_list.append(current_species)
            
            current_species = {
                'name': line.strip(),
                'charge': None,
                'ion_size': None,
                'mole_weight': None,
                'formula': None,
                'mole_vol': None
            }
        
        # Parse charge, ion size, mole weight
        elif current_species and 'charge=' in line:
            charge_match = re.search(r'charge=\s*([-\d]+)', line)
            size_match = re.search(r'ion size=\s*([\d.]+)', line)
            mw_match = re.search(r'mole wt\.=\s*([\d.]+)', line)
            
            if charge_match:
                current_species['charge'] = int(charge_match.group(1))
            if size_match:
                current_species['ion_size'] = float(size_match.group(1))
            if mw_match:
                current_species['mole_weight'] = float(mw_match.group(1))
        
        # Parse formula for minerals
        elif current_species and 'formula=' in line:
            formula_match = re.search(r'formula=\s*(\S+)', line)
            if formula_match:
                current_species['formula'] = formula_match.group(1)
        
        # Parse mole volume for minerals
        elif current_species and 'mole vol' in line:
            vol_match = re.search(r'mole vol\.=\s*([\d.]+)', line)
            if vol_match:
                current_species['mole_vol'] = float(vol_match.group(1))
        
        i += 1
    
    return species_list


def parse_pitzer_binary_parameters(lines, start_idx):
    """Parse Pitzer binary (ca) parameters."""
    params = []
    i = start_idx
    
    # Find start of ca combinations
    while i < len(lines):
        if 'ca combinations' in lines[i].lower() or re.match(r'^[A-Za-z0-9\(\)\+\-]+\s+[A-Za-z0-9\(\)\+\-]+\s*$', lines[i].strip()):
            break
        i += 1
    
    current_param = None
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for end of binary parameters section
        if '-end-' in line and 'lambda' in lines[i-1:i+5] or 'theta' in str(lines[i-1:i+5]):
            if current_param:
                params.append(current_param)
            break
        
        # New parameter pair (cation-anion or species pair)
        if line and not line.startswith('*') and not line.startswith('beta') and not line.startswith('cphi') and not line.startswith('alpha'):
            match = re.match(r'^([A-Za-z0-9\(\)\+\-]+)\s+([A-Za-z0-9\(\)\+\-]+)\s*$', line)
            if match:
                if current_param:
                    params.append(current_param)
                
                current_param = {
                    'species1': match.group(1),
                    'species2': match.group(2),
                    'beta0': None,
                    'beta1': None,
                    'beta2': None,
                    'cphi': None,
                    'alpha1': None,
                    'alpha2': None,
                    'source': None
                }
        
        # Parse parameters
        elif current_param:
            if line.startswith('beta0'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['beta0'] = ' '.join(vals)
            elif line.startswith('beta1'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['beta1'] = ' '.join(vals)
            elif line.startswith('beta2'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['beta2'] = ' '.join(vals)
            elif line.startswith('cphi'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['cphi'] = ' '.join(vals)
            elif line.startswith('alpha1'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['alpha1'] = vals[0]
            elif line.startswith('alpha2'):
                vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
                if vals:
                    current_param['alpha2'] = vals[0]
            elif line.startswith('* Source:'):
                current_param['source'] = line.replace('* Source:', '').strip()
        
        i += 1
    
    return params


def parse_theta_parameters(lines, start_idx):
    """Parse theta parameters (like-charge interactions)."""
    params = []
    i = start_idx
    
    # Find theta section
    while i < len(lines):
        if 'begin with theta' in lines[i] or (i > 0 and '-end-' in lines[i-1] and 'theta' in lines[i+1]):
            i += 1
            break
        i += 1
    
    current_param = None
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for end of theta section
        if '-end-' in line and 'lambda' in str(lines[i:i+5]):
            if current_param:
                params.append(current_param)
            break
        
        # New theta pair
        if line and not line.startswith('*') and not line.startswith('theta'):
            match = re.match(r'^([A-Za-z0-9\(\)\+\-]+)\s+([A-Za-z0-9\(\)\+\-]+)\s*$', line)
            if match:
                if current_param:
                    params.append(current_param)
                
                current_param = {
                    'species1': match.group(1),
                    'species2': match.group(2),
                    'theta': None,
                    'source': None
                }
        
        # Parse theta value
        elif current_param and line.startswith('theta'):
            vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
            if vals:
                current_param['theta'] = ' '.join(vals)
        
        # Parse source
        elif current_param and line.startswith('* Source:'):
            current_param['source'] = line.replace('* Source:', '').strip()
        
        i += 1
    
    return params


def parse_lambda_parameters(lines, start_idx):
    """Parse lambda parameters (neutral species interactions)."""
    params = []
    i = start_idx
    
    # Find lambda section
    while i < len(lines):
        if 'begin with lambda' in lines[i]:
            i += 1
            break
        i += 1
    
    current_param = None
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for end of lambda section
        if '-end-' in line and 'psi' in str(lines[i:i+5]):
            if current_param:
                params.append(current_param)
            break
        
        # New lambda pair
        if line and not line.startswith('*') and not line.startswith('lambda'):
            match = re.match(r'^([A-Za-z0-9\(\)\+\-]+)\s+([A-Za-z0-9\(\)\+\-]+)\s*$', line)
            if match:
                if current_param:
                    params.append(current_param)
                
                current_param = {
                    'species1': match.group(1),
                    'species2': match.group(2),
                    'lambda': None,
                    'source': None
                }
        
        # Parse lambda value
        elif current_param and line.startswith('lambda'):
            vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
            if vals:
                current_param['lambda'] = ' '.join(vals)
        
        # Parse source
        elif current_param and line.startswith('* Source:'):
            current_param['source'] = line.replace('* Source:', '').strip()
        
        i += 1
    
    return params


def parse_psi_parameters(lines, start_idx):
    """Parse psi parameters (ternary interactions)."""
    params = []
    i = start_idx
    
    # Find psi section
    while i < len(lines):
        if 'begin with psi' in lines[i]:
            i += 1
            break
        i += 1
    
    current_param = None
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for end of psi section
        if '-end-' in line:
            if current_param:
                params.append(current_param)
            break
        
        # New psi triplet
        if line and not line.startswith('*') and not line.startswith('psi'):
            match = re.match(r'^([A-Za-z0-9\(\)\+\-]+)\s+([A-Za-z0-9\(\)\+\-]+)\s+([A-Za-z0-9\(\)\+\-]+)\s*$', line)
            if match:
                if current_param:
                    params.append(current_param)
                
                current_param = {
                    'species1': match.group(1),
                    'species2': match.group(2),
                    'species3': match.group(3),
                    'psi': None,
                    'source': None
                }
        
        # Parse psi value
        elif current_param and line.startswith('psi'):
            vals = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', line)
            if vals:
                current_param['psi'] = ' '.join(vals)
        
        # Parse source
        elif current_param and line.startswith('* Source:'):
            current_param['source'] = line.replace('* Source:', '').strip()
        
        i += 1
    
    return params


def main():
    # Read the file (assumes script is run from data/ directory)
    db_file = Path('thermo_ymp.R2.tdat')
    print(f"Reading {db_file}...")
    
    if not db_file.exists():
        raise FileNotFoundError(f"Could not find {db_file}. Please run this script from the data/ directory.")
    
    # Try different encodings
    for encoding in ['utf-8', 'latin-1', 'cp1252']:
        try:
            with open(db_file, 'r', encoding=encoding) as f:
                lines = f.readlines()
            print(f"Successfully read file with {encoding} encoding")
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError("Could not decode file with any standard encoding")
    
    print(f"Total lines: {len(lines)}")
    
    # Create output directory within data/
    output_dir = Path('parsed_thermodb')
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Parse elements
    print("\nParsing elements...")
    elements = parse_elements(lines, 0)
    print(f"Found {len(elements)} elements")
    
    with open(output_dir / 'elements.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['name', 'symbol', 'mole_weight'])
        writer.writeheader()
        writer.writerows(elements)
    
    # Parse basis species
    print("Parsing basis species...")
    basis_species = parse_species_section(lines, 0, 'basis species')
    print(f"Found {len(basis_species)} basis species")
    
    with open(output_dir / 'basis_species.csv', 'w', newline='') as f:
        fieldnames = ['name', 'charge', 'ion_size', 'mole_weight']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(basis_species)
    
    # Parse redox couples
    print("Parsing redox couples...")
    redox_couples = parse_species_section(lines, 0, 'redox couples')
    print(f"Found {len(redox_couples)} redox couples")
    
    with open(output_dir / 'redox_couples.csv', 'w', newline='') as f:
        fieldnames = ['name', 'charge', 'ion_size', 'mole_weight']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(redox_couples)
    
    # Parse aqueous species
    print("Parsing aqueous species...")
    aqueous_species = parse_species_section(lines, 0, 'aqueous species')
    print(f"Found {len(aqueous_species)} aqueous species")
    
    with open(output_dir / 'aqueous_species.csv', 'w', newline='') as f:
        fieldnames = ['name', 'charge', 'ion_size', 'mole_weight']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(aqueous_species)
    
    # Parse minerals
    print("Parsing minerals...")
    minerals = parse_species_section(lines, 0, 'minerals')
    print(f"Found {len(minerals)} minerals")
    
    with open(output_dir / 'minerals.csv', 'w', newline='') as f:
        fieldnames = ['name', 'formula', 'mole_vol', 'mole_weight']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(minerals)
    
    # Parse gases
    print("Parsing gases...")
    gases = parse_species_section(lines, 0, 'gases')
    print(f"Found {len(gases)} gases")
    
    with open(output_dir / 'gases.csv', 'w', newline='') as f:
        fieldnames = ['name', 'mole_weight']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(gases)
    
    # Parse Pitzer parameters
    print("\nParsing Pitzer binary parameters...")
    binary_params = parse_pitzer_binary_parameters(lines, 7000)
    print(f"Found {len(binary_params)} binary parameter sets")
    
    with open(output_dir / 'pitzer_binary.csv', 'w', newline='') as f:
        fieldnames = ['species1', 'species2', 'beta0', 'beta1', 'beta2', 
                      'cphi', 'alpha1', 'alpha2', 'source']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(binary_params)
    
    print("Parsing Pitzer theta parameters...")
    theta_params = parse_theta_parameters(lines, 9000)
    print(f"Found {len(theta_params)} theta parameter sets")
    
    with open(output_dir / 'pitzer_theta.csv', 'w', newline='') as f:
        fieldnames = ['species1', 'species2', 'theta', 'source']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(theta_params)
    
    print("Parsing Pitzer lambda parameters...")
    lambda_params = parse_lambda_parameters(lines, 9400)
    print(f"Found {len(lambda_params)} lambda parameter sets")
    
    with open(output_dir / 'pitzer_lambda.csv', 'w', newline='') as f:
        fieldnames = ['species1', 'species2', 'lambda', 'source']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(lambda_params)
    
    print("Parsing Pitzer psi parameters...")
    psi_params = parse_psi_parameters(lines, 9700)
    print(f"Found {len(psi_params)} psi parameter sets")
    
    with open(output_dir / 'pitzer_psi.csv', 'w', newline='') as f:
        fieldnames = ['species1', 'species2', 'species3', 'psi', 'source']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(psi_params)
    
    print(f"\n✓ All data exported to {output_dir}/")
    print("\nGenerated CSV files:")
    for csv_file in sorted(output_dir.glob('*.csv')):
        size = csv_file.stat().st_size / 1024
        print(f"  - {csv_file.name} ({size:.1f} KB)")


if __name__ == '__main__':
    main()
