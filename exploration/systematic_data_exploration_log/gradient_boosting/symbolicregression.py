"""
Symbolic Regression for Hygroscopicity (PySR)
---------------------------------------------
Goal: Discover an interpretable physical equation that predicts
      the activity coefficient (ln_gamma) at 90% RH based on ion properties.

Target: ln_gamma_at_90RH
        (Lower values = More Hygroscopic/Stable in solution)

Feature Engineering:
    - Lattice Energy Proxy (Kapustinskii)
    - Enthalpy of Solution Proxy
    - Water Affinity Mismatch
    - Coulomb Interaction

Author: Generated for AWH ML Discussion
Date: 2026-01-27
"""

import pandas as pd
import numpy as np
import os
from pysr import PySRRegressor

# --- CONFIGURATION ---
DATA_FILE = '../../../data/baseline_numeric_only.csv'
OUTPUT_DIR = '../../../figures/pysr_hygroscopicity/'
TARGET_COL = 'ln_gamma_at_90RH'

# PySR Settings
N_ITERATIONS = 40
POPULATIONS = 20
BINARY_OPERATORS = ["+", "-", "*", "/"]
UNARY_OPERATORS = ["exp", "abs", "sqrt", "inv(x) = 1/x", "square(x) = x^2"]

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_and_engineer_features(filepath):
    print(f"Loading data from {filepath}...")
    
    # Read the file
    # Assuming the file structure matches previous descriptions (Header on line 0)
    df = pd.read_csv(filepath)
    
    # Filter for valid target rows
    df = df.dropna(subset=[TARGET_COL])
    
    print(f"Original columns: {len(df.columns)}")
    
    # --- FEATURE ENGINEERING (The 'Secret Sauce') ---
    
    # 1. Coulomb Interaction Term (z+ z- / r+ + r-)
    # We approximate z based on electrolyte type or assume 1 if not explicitly available.
    # Here we use the raw radii provided: r_M_angstrom (Cation), r_X_angstrom (Anion)
    # Note: 'r_M' and 'r_X' are typically Crystal Radii.
    
    df['inv_sum_radii'] = 1.0 / (df['r_M_angstrom'] + df['r_X_angstrom'])
    
    # 2. Lattice Energy Proxy (Kapustinskii-like)
    # U_lat ~ 1 / (r_c + r_a)
    df['feat_lattice_proxy'] = df['inv_sum_radii']
    
    # 3. Hydration Energy Sum (Total interactions with water)
    df['feat_hydration_sum'] = df['cation_1_delta_G_hydration'] + df['anion_1_delta_G_hydration']
    
    # 4. Enthalpy of Solution Proxy (Competition: Lattice vs Hydration)
    # This is a critical term for solubility/hygroscopicity.
    # Note: Using '+' because hydration G is typically negative.
    # If |Hydration| > |Lattice|, salt dissolves well.
    df['feat_enthalpy_sol_proxy'] = df['feat_lattice_proxy'] + (df['feat_hydration_sum'] / 1000.0) 
    
    # 5. Law of Matching Water Affinities (Mismatch)
    # | DeltaG_cat - DeltaG_an |
    df['feat_G_mismatch'] = (df['cation_1_delta_G_hydration'] - df['anion_1_delta_G_hydration']).abs()
    
    # 6. Ionic Volume / Size mismatch
    df['feat_radius_ratio'] = df['r_M_angstrom'] / df['r_X_angstrom']
    
    # --- SELECT FEATURES FOR MODEL ---
    # We select the engineered features + key raw physical properties
    feature_cols = [
        # Engineered High-Value Features
        'feat_lattice_proxy',
        'feat_hydration_sum',
        'feat_enthalpy_sol_proxy',
        'feat_G_mismatch',
        'feat_radius_ratio',
        
        # Raw Physical Properties (giving the model flexibility)
        'r_M_angstrom',
        'r_X_angstrom',
        'cation_1_delta_G_hydration',
        'anion_1_delta_G_hydration',
        'cation_1_viscosity_jones_dole',
        'anion_1_viscosity_jones_dole'
    ]
    
    X = df[feature_cols]
    y = df[TARGET_COL]
    
    print(f"Engineered {len(feature_cols)} physical features.")
    return X, y, feature_cols

def run_symbolic_regression():
    # 1. Prepare Data
    X, y, feature_names = load_and_engineer_features(DATA_FILE)
    
    # 2. Initialize Model
    print("\nInitializing PySR Regressor...")
    print("Searching for equations... (This may take a few minutes)")
    
    model = PySRRegressor(
        niterations=N_ITERATIONS,
        populations=POPULATIONS,
        binary_operators=BINARY_OPERATORS,
        unary_operators=UNARY_OPERATORS,
        model_selection="best",  # Selects best mix of accuracy vs complexity
        temp_equation_file=os.path.join(OUTPUT_DIR, "hall_of_fame.csv"),
        feature_names=feature_names,
        equation_file=os.path.join(OUTPUT_DIR, "pysr_equations.csv")
    )
    
    # 3. Fit
    model.fit(X, y)
    
    # 4. Results
    print("\n" + "="*60)
    print("DISCOVERED EQUATIONS (Ranked by Score)")
    print("="*60)
    print(model)
    
    print("\n" + "="*60)
    print("BEST EQUATION (Trade-off between simplicity and accuracy):")
    print("="*60)
    print(model.sympy())
    
    # 5. Save Analysis
    print(f"\nSaving full equation list to: {OUTPUT_DIR}pysr_equations.csv")

if __name__ == "__main__":
    run_symbolic_regression()