"""
Data Merger for Master Equation Symbolic Regression
---------------------------------------------------
1. Loads full activity curves (Activity Data).
2. Loads ion properties (Properties Data).
3. Normalizes salt names (removes '()', handles 'ZnCl' vs 'ZnCl2').
4. Downsamples activity curves (100 points -> 10 points).
5. Merges datasets into a single Master CSV.

Author: Generated for AWH ML Discussion
Date: 2026-01-28
"""

import pandas as pd
import numpy as np
import re

# --- CONFIGURATION ---
# UPDATE THESE FILENAMES TO MATCH YOUR ACTUAL FILES
ACTIVITY_FILE = '../../../data/water_activity_all_salts_combined.csv'  # The file with 100 rows per salt
PROPERTIES_FILE = '../../../data/baseline_numeric_only.csv' # The file with ion properties
OUTPUT_FILE = '../../../data/master_long_format.csv'

POINTS_PER_SALT = 10  # Downsample target (e.g., 10 points per curve)

def normalize_name(name):
    """
    Basic cleanup: remove parentheses and whitespace.
    """
    if not isinstance(name, str):
        return str(name)
    # Remove () and whitespace
    clean = name.replace('(', '').replace(')', '').strip()
    return clean

def get_name_mapping(activity_salts, prop_salts):
    """
    Creates a dictionary mapping: {Activity_File_Name -> Property_File_Name}
    Logic:
      1. Exact match (after normalization).
      2. "Missing 2" match (e.g., ZnCl matches ZnCl2).
    """
    mapping = {}
    
    # Pre-normalize property salts for faster lookup
    prop_norm_map = {normalize_name(s): s for s in prop_salts}
    prop_norm_set = set(prop_norm_map.keys())
    
    print(f"Finding matches for {len(activity_salts)} salts...")
    
    for salt in activity_salts:
        norm = normalize_name(salt)
        
        match_found = None
        
        # 1. Direct Match
        if norm in prop_norm_set:
            match_found = prop_norm_map[norm]
            
        # 2. "Missing 2" Match (Check if adding or removing a '2' helps)
        elif (norm + "2") in prop_norm_set:
            match_found = prop_norm_map[norm + "2"]
            # print(f"  Matched {salt} -> {match_found} (Added '2')")
            
        elif norm.endswith("2") and (norm[:-1] in prop_norm_set):
            match_found = prop_norm_map[norm[:-1]]
            # print(f"  Matched {salt} -> {match_found} (Removed '2')")
            
        if match_found:
            mapping[salt] = match_found
            
    print(f"  -> Found {len(mapping)} valid matches between datasets.")
    return mapping

def main():
    # 1. Load Data
    print("Loading datasets...")
    try:
        df_act = pd.read_csv(ACTIVITY_FILE)
        df_prop = pd.read_csv(PROPERTIES_FILE)
    except FileNotFoundError as e:
        print(f"\nERROR: Could not find file. {e}")
        print("Please ensure the filenames in the 'CONFIGURATION' section match your files.")
        return

    # 2. Downsample Activity Data (100 rows -> 10 rows evenly spaced)
    print(f"Downsampling activity curves to ~{POINTS_PER_SALT} points per salt...")
    
    # We group by Salt and pick indices
    def get_evenly_spaced(group):
        if len(group) <= POINTS_PER_SALT:
            return group
        # Select evenly spaced indices
        indices = np.round(np.linspace(0, len(group) - 1, POINTS_PER_SALT)).astype(int)
        return group.iloc[indices]

    # Assuming column name is 'Salt' in activity file (based on your image)
    # Adjust if different (e.g., 'electrolyte')
    salt_col_act = 'Salt' 
    if salt_col_act not in df_act.columns:
        # Fallback: try to find a likely column
        possible = [c for c in df_act.columns if 'salt' in c.lower() or 'electrolyte' in c.lower()]
        if possible:
            salt_col_act = possible[0]
            print(f"  Assuming salt column in Activity file is: '{salt_col_act}'")
        else:
            raise ValueError("Could not identify Salt column in Activity file.")

    df_act_down = df_act.groupby(salt_col_act, group_keys=False).apply(get_evenly_spaced)
    
    # 3. Match Names
    act_salts = df_act_down[salt_col_act].unique()
    prop_salts = df_prop['electrolyte'].unique() # Assuming 'electrolyte' is col name in props file
    
    name_map = get_name_mapping(act_salts, prop_salts)
    
    # 4. Filter and Merge
    # Create a 'merge_key' column in activity data
    df_act_down['merge_key'] = df_act_down[salt_col_act].map(name_map)
    
    # Drop rows that didn't match any property (salt not in second database)
    df_act_final = df_act_down.dropna(subset=['merge_key'])
    
    # Merge!
    # Inner join ensures we only keep salts present in BOTH
    print("Merging datasets...")
    df_master = pd.merge(
        df_act_final,
        df_prop,
        left_on='merge_key',
        right_on='electrolyte',
        how='inner'
    )
    
    # 5. Clean up
    # Drop the temporary merge key if desired
    df_master = df_master.drop(columns=['merge_key'])
    
    # 6. Save
    print(f"Saving Master Dataset to {OUTPUT_FILE}...")
    df_master.to_csv(OUTPUT_FILE, index=False)
    
    print("Done!")
    print(f"Final shape: {df_master.shape} (Rows, Cols)")
    print(f"Unique Salts included: {df_master['electrolyte'].nunique()}")

if __name__ == "__main__":
    main()