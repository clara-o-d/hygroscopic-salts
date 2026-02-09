#!/usr/bin/env python3
"""
Process Pitzer Parameter Data for Machine Learning

This script combines Pitzer parameter data (psi, lambda, theta) with ion properties
from baseline_numeric_only.csv to create ML-ready datasets for predicting
Pitzer parameters at 25°C.

For each parameter type:
- psi: 3-species interaction (cation-anion-cation or cation-anion-anion)
- lambda: 2-species interaction (neutral-ion or ion-neutral)
- theta: 2-species interaction (cation-cation or anion-anion)

Output: Processed CSV files in respective directories (lambda/, psi/, theta/)
"""

import pandas as pd
import numpy as np
import os
import re
from collections import defaultdict


def parse_ion_name(ion_string):
    """Parse ion name to get element and charge."""
    # Remove (aq) suffix if present
    ion_string = ion_string.replace('(aq)', '')
    
    # Try to find charge pattern
    # Patterns: Na+, Ca++, Cl-, SO4--, etc.
    charge_match = re.search(r'(\++|--|-|\+)$', ion_string)
    
    if charge_match:
        charge_str = charge_match.group(1)
        element = ion_string[:charge_match.start()]
        
        # Convert charge string to integer
        if '+' in charge_str:
            charge = charge_str.count('+')
        else:
            charge = -charge_str.count('-')
    else:
        # Neutral species
        element = ion_string
        charge = 0
    
    return element, charge


def extract_ion_properties_from_baseline(baseline_df):
    """
    Extract individual ion properties from baseline_numeric_only.csv.
    
    The baseline file has electrolyte pairs with separate cation and anion properties.
    We'll create a dictionary mapping ion names to their properties.
    """
    
    ion_properties = {}
    
    # Ion property columns in baseline
    cation_props = [col for col in baseline_df.columns if col.startswith('cation_1_')]
    anion_props = [col for col in baseline_df.columns if col.startswith('anion_1_')]
    
    # Property names without prefix
    prop_names_cation = [col.replace('cation_1_', '') for col in cation_props]
    prop_names_anion = [col.replace('anion_1_', '') for col in anion_props]
    
    # Process each electrolyte
    for idx, row in baseline_df.iterrows():
        electrolyte = row['electrolyte']
        
        # Parse electrolyte name to get cation and anion
        # Common patterns: NaCl, CaCl2, Na2SO4, etc.
        cation_name, anion_name = parse_electrolyte_name(electrolyte)
        
        if cation_name and cation_name not in ion_properties:
            ion_properties[cation_name] = {}
            for prop, col in zip(prop_names_cation, cation_props):
                ion_properties[cation_name][prop] = row[col]
            ion_properties[cation_name]['charge'] = determine_charge_from_name(cation_name)
        
        if anion_name and anion_name not in ion_properties:
            ion_properties[anion_name] = {}
            for prop, col in zip(prop_names_anion, anion_props):
                ion_properties[anion_name][prop] = row[col]
            ion_properties[anion_name]['charge'] = determine_charge_from_name(anion_name)
    
    return ion_properties


def parse_electrolyte_name(electrolyte):
    """
    Parse electrolyte name to extract cation and anion.
    
    Examples:
    - NaCl -> Na+, Cl-
    - CaCl2 -> Ca++, Cl-
    - Na2SO4 -> Na+, SO4--
    - MgSO4 -> Mg++, SO4--
    """
    
    # Common electrolyte patterns
    patterns = {
        # Format: (cation, anion, cation_charge, anion_charge)
        # 1-1 electrolytes
        r'^([A-Z][a-z]?)Cl$': (r'\1', 'Cl', 1, -1),
        r'^([A-Z][a-z]?)Br$': (r'\1', 'Br', 1, -1),
        r'^([A-Z][a-z]?)I$': (r'\1', 'I', 1, -1),
        r'^([A-Z][a-z]?)F$': (r'\1', 'F', 1, -1),
        r'^([A-Z][a-z]?)NO3$': (r'\1', 'NO3', 1, -1),
        r'^([A-Z][a-z]?)NO2$': (r'\1', 'NO2', 1, -1),
        r'^([A-Z][a-z]?)OH$': (r'\1', 'OH', 1, -1),
        r'^([A-Z][a-z]?)ClO4$': (r'\1', 'ClO4', 1, -1),
        r'^([A-Z][a-z]?)ClO3$': (r'\1', 'ClO3', 1, -1),
        r'^([A-Z][a-z]?)BrO3$': (r'\1', 'BrO3', 1, -1),
        r'^([A-Z][a-z]?)SCN$': (r'\1', 'SCN', 1, -1),
        r'^([A-Z][a-z]?)H2PO4$': (r'\1', 'H2PO4', 1, -1),
        r'^([A-Z][a-z]?)HCO3$': (r'\1', 'HCO3', 1, -1),
        
        # 1-2 electrolytes (M2X)
        r'^([A-Z][a-z]?)2SO4$': (r'\1', 'SO4', 1, -2),
        r'^([A-Z][a-z]?)2CO3$': (r'\1', 'CO3', 1, -2),
        r'^([A-Z][a-z]?)2S$': (r'\1', 'S', 1, -2),
        
        # 2-1 electrolytes (MX2)
        r'^([A-Z][a-z]?)Cl2$': (r'\1', 'Cl', 2, -1),
        r'^([A-Z][a-z]?)Br2$': (r'\1', 'Br', 2, -1),
        r'^([A-Z][a-z]?)I2$': (r'\1', 'I', 2, -1),
        r'^([A-Z][a-z]?)\(NO3\)2$': (r'\1', 'NO3', 2, -1),
        r'^([A-Z][a-z]?)\(ClO4\)2$': (r'\1', 'ClO4', 2, -1),
        
        # 2-2 electrolytes
        r'^([A-Z][a-z]?)SO4$': (r'\1', 'SO4', 2, -2),
        r'^([A-Z][a-z]?)CO3$': (r'\1', 'CO3', 2, -2),
        
        # 3-1 electrolytes
        r'^([A-Z][a-z]?)Cl3$': (r'\1', 'Cl', 3, -1),
        
        # 4-1 electrolytes  
        r'^([A-Z][a-z]?)\(NO3\)4$': (r'\1', 'NO3', 4, -1),
        
        # Hydrogen compounds
        r'^HCl$': ('H', 'Cl', 1, -1),
        r'^HBr$': ('H', 'Br', 1, -1),
        r'^HI$': ('H', 'I', 1, -1),
        r'^HNO3$': ('H', 'NO3', 1, -1),
        r'^HClO4$': ('H', 'ClO4', 1, -1),
        
        # Ammonium compounds
        r'^NH4Cl$': ('NH4', 'Cl', 1, -1),
        r'^NH4NO3$': ('NH4', 'NO3', 1, -1),
        r'^NH4ClO4$': ('NH4', 'ClO4', 1, -1),
        r'^\(NH4\)2SO4$': ('NH4', 'SO4', 1, -2),
    }
    
    for pattern, (cat, an, cat_charge, an_charge) in patterns.items():
        match = re.match(pattern, electrolyte)
        if match:
            cation = match.expand(cat)
            anion = an
            # Add charge symbols
            cation_name = cation + '+' * cat_charge
            anion_name = anion + '-' * abs(an_charge)
            return cation_name, anion_name
    
    # If no pattern matches, return None
    return None, None


def determine_charge_from_name(ion_name):
    """Determine charge from ion name."""
    if '+' in ion_name:
        return ion_name.count('+')
    elif '-' in ion_name:
        return -ion_name.count('-')
    else:
        return 0


def get_ion_properties(ion_name, ion_properties_dict):
    """Get properties for a specific ion, return dict of properties or None."""
    # Try exact match first
    if ion_name in ion_properties_dict:
        return ion_properties_dict[ion_name]
    
    # Try without (aq) suffix
    clean_name = ion_name.replace('(aq)', '')
    if clean_name in ion_properties_dict:
        return ion_properties_dict[clean_name]
    
    # If not found, return None
    return None


def process_lambda_data(lambda_df, ion_properties):
    """
    Process lambda parameters (neutral-ion interactions).
    
    Lambda involves 2 species: typically one neutral and one charged.
    """
    print("\n" + "="*80)
    print("Processing Lambda Parameters")
    print("="*80)
    
    processed_rows = []
    
    for idx, row in lambda_df.iterrows():
        species1 = row['species1']
        species2 = row['species2']
        
        # Get properties for both species
        props1 = get_ion_properties(species1, ion_properties)
        props2 = get_ion_properties(species2, ion_properties)
        
        # Create feature row
        feature_row = {
            'species1': species1,
            'species2': species2,
            'lambda_a1': row['lambda_a1'],
            'lambda_a2': row.get('lambda_a2', np.nan),
            'lambda_a3': row.get('lambda_a3', np.nan),
            'lambda_a4': row.get('lambda_a4', np.nan),
            'source': row['source']
        }
        
        # Add species1 properties
        if props1:
            for prop_name, prop_value in props1.items():
                feature_row[f'species1_{prop_name}'] = prop_value
        
        # Add species2 properties
        if props2:
            for prop_name, prop_value in props2.items():
                feature_row[f'species2_{prop_name}'] = prop_value
        
        processed_rows.append(feature_row)
    
    df = pd.DataFrame(processed_rows)
    print(f"✓ Processed {len(df)} lambda parameter entries")
    print(f"  Features: {len(df.columns)} columns")
    
    return df


def process_psi_data(psi_df, ion_properties):
    """
    Process psi parameters (triple ion interactions).
    
    Psi involves 3 species: typically cation-anion-cation or cation-anion-anion.
    """
    print("\n" + "="*80)
    print("Processing Psi Parameters")
    print("="*80)
    
    processed_rows = []
    
    for idx, row in psi_df.iterrows():
        species1 = row['species1']
        species2 = row['species2']
        species3 = row['species3']
        
        # Get properties for all three species
        props1 = get_ion_properties(species1, ion_properties)
        props2 = get_ion_properties(species2, ion_properties)
        props3 = get_ion_properties(species3, ion_properties)
        
        # Create feature row
        feature_row = {
            'species1': species1,
            'species2': species2,
            'species3': species3,
            'psi_a1': row['psi_a1'],
            'psi_a2': row.get('psi_a2', np.nan),
            'psi_a3': row.get('psi_a3', np.nan),
            'psi_a4': row.get('psi_a4', np.nan),
            'source': row['source']
        }
        
        # Add species1 properties
        if props1:
            for prop_name, prop_value in props1.items():
                feature_row[f'species1_{prop_name}'] = prop_value
        
        # Add species2 properties
        if props2:
            for prop_name, prop_value in props2.items():
                feature_row[f'species2_{prop_name}'] = prop_value
        
        # Add species3 properties
        if props3:
            for prop_name, prop_value in props3.items():
                feature_row[f'species3_{prop_name}'] = prop_value
        
        processed_rows.append(feature_row)
    
    df = pd.DataFrame(processed_rows)
    print(f"✓ Processed {len(df)} psi parameter entries")
    print(f"  Features: {len(df.columns)} columns")
    
    return df


def process_theta_data(theta_df, ion_properties):
    """
    Process theta parameters (same-charge ion interactions).
    
    Theta involves 2 species of the same charge sign.
    """
    print("\n" + "="*80)
    print("Processing Theta Parameters")
    print("="*80)
    
    processed_rows = []
    
    for idx, row in theta_df.iterrows():
        species1 = row['species1']
        species2 = row['species2']
        
        # Get properties for both species
        props1 = get_ion_properties(species1, ion_properties)
        props2 = get_ion_properties(species2, ion_properties)
        
        # Create feature row
        feature_row = {
            'species1': species1,
            'species2': species2,
            'theta_a1': row['theta_a1'],
            'theta_a2': row.get('theta_a2', np.nan),
            'theta_a3': row.get('theta_a3', np.nan),
            'theta_a4': row.get('theta_a4', np.nan),
            'source': row['source']
        }
        
        # Add species1 properties
        if props1:
            for prop_name, prop_value in props1.items():
                feature_row[f'species1_{prop_name}'] = prop_value
        
        # Add species2 properties
        if props2:
            for prop_name, prop_value in props2.items():
                feature_row[f'species2_{prop_name}'] = prop_value
        
        processed_rows.append(feature_row)
    
    df = pd.DataFrame(processed_rows)
    print(f"✓ Processed {len(df)} theta parameter entries")
    print(f"  Features: {len(df.columns)} columns")
    
    return df


def main():
    """Main execution function."""
    print("\n" + "="*80)
    print("PITZER PARAMETER DATA PROCESSING")
    print("Combining Parameter Data with Ion Properties")
    print("="*80)
    
    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, '../data')
    parsed_dir = os.path.join(data_dir, 'parsed_thermodb')
    
    # Load baseline data
    print("\n" + "="*80)
    print("Loading Data")
    print("="*80)
    
    baseline_path = os.path.join(data_dir, 'baseline_numeric_only.csv')
    baseline_df = pd.read_csv(baseline_path)
    print(f"✓ Loaded baseline data: {len(baseline_df)} electrolytes")
    
    # Extract ion properties
    print("\nExtracting ion properties from baseline data...")
    ion_properties = extract_ion_properties_from_baseline(baseline_df)
    print(f"✓ Extracted properties for {len(ion_properties)} unique ions")
    print(f"  Sample ions: {', '.join(list(ion_properties.keys())[:10])}...")
    
    # Load parameter data
    lambda_path = os.path.join(parsed_dir, 'pitzer_lambda.csv')
    psi_path = os.path.join(parsed_dir, 'pitzer_psi.csv')
    theta_path = os.path.join(parsed_dir, 'pitzer_theta.csv')
    
    lambda_df = pd.read_csv(lambda_path)
    psi_df = pd.read_csv(psi_path)
    theta_df = pd.read_csv(theta_path)
    
    print(f"\n✓ Loaded lambda data: {len(lambda_df)} entries")
    print(f"✓ Loaded psi data: {len(psi_df)} entries")
    print(f"✓ Loaded theta data: {len(theta_df)} entries")
    
    # Process each parameter type
    lambda_processed = process_lambda_data(lambda_df, ion_properties)
    psi_processed = process_psi_data(psi_df, ion_properties)
    theta_processed = process_theta_data(theta_df, ion_properties)
    
    # Save processed data
    print("\n" + "="*80)
    print("Saving Processed Data")
    print("="*80)
    
    lambda_out = os.path.join(script_dir, 'lambda', 'lambda_processed.csv')
    psi_out = os.path.join(script_dir, 'psi', 'psi_processed.csv')
    theta_out = os.path.join(script_dir, 'theta', 'theta_processed.csv')
    
    lambda_processed.to_csv(lambda_out, index=False)
    print(f"✓ Saved lambda data to: {lambda_out}")
    print(f"  Shape: {lambda_processed.shape}")
    
    psi_processed.to_csv(psi_out, index=False)
    print(f"✓ Saved psi data to: {psi_out}")
    print(f"  Shape: {psi_processed.shape}")
    
    theta_processed.to_csv(theta_out, index=False)
    print(f"✓ Saved theta data to: {theta_out}")
    print(f"  Shape: {theta_processed.shape}")
    
    # Summary
    print("\n" + "="*80)
    print("PROCESSING COMPLETE!")
    print("="*80)
    print(f"\nGenerated Files:")
    print(f"  - lambda/lambda_processed.csv ({len(lambda_processed)} samples)")
    print(f"  - psi/psi_processed.csv ({len(psi_processed)} samples)")
    print(f"  - theta/theta_processed.csv ({len(theta_processed)} samples)")
    print("\n" + "="*80)


if __name__ == "__main__":
    main()
