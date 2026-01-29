"""
Symbolic Regression for Hygroscopicity (Robust / High-Complexity)
-----------------------------------------------------------------
Goal: Discover an interpretable physical equation for ln_gamma at 90% RH.
      * ROBUST MODE: Uses L1 Loss to ignore outliers.
      * COMPLEX MODE: Allows longer, more intricate equations.

Target: ln_gamma_at_90RH

Author: Generated for AWH ML Discussion
Date: 2026-01-28
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pysr import PySRRegressor

# --- CONFIGURATION ---
DATA_FILE = '../../../data/baseline_numeric_only.csv'
OUTPUT_DIR = '../../../figures/symbolic_regression/'
TARGET_COL = 'ln_gamma_at_90RH'

# --- IMPROVED SETTINGS ---
N_ITERATIONS = 200          # Increased from 40 (search longer)
POPULATIONS = 30            # Increased diversity
MAX_EQUATION_LENGTH = 30    # Allow more complicated formulas
LOSS_FUNCTION = "L1DistLoss()" # ROBUST LOSS: Ignores outliers (Mean Absolute Error)

BINARY_OPERATORS = ["+", "-", "*", "/", "^"] # Added power (^)
UNARY_OPERATORS = ["exp", "abs", "sqrt", "log", "inv(x) = 1/x"]

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_and_engineer_features(filepath):
    print(f"Loading data from {filepath}...")
    df = pd.read_csv(filepath)
    df = df.dropna(subset=[TARGET_COL])
    
    # --- FEATURE ENGINEERING ---
    # 1. Coulomb / Lattice Terms
    df['inv_sum_radii'] = 1.0 / (df['r_M_angstrom'] + df['r_X_angstrom'])
    df['feat_lattice_proxy'] = df['inv_sum_radii']
    
    # 2. Hydration Energy Terms
    df['feat_hydration_sum'] = df['cation_1_delta_G_hydration'] + df['anion_1_delta_G_hydration']
    
    # 3. Competition Terms (Enthalpy of Solution Proxy)
    # This roughly correlates with solubility
    df['feat_enthalpy_sol_proxy'] = df['feat_lattice_proxy'] + (df['feat_hydration_sum'] / 1000.0) 
    
    # 4. Mismatch Terms (Matching Water Affinities)
    df['feat_G_mismatch'] = (df['cation_1_delta_G_hydration'] - df['anion_1_delta_G_hydration']).abs()
    
    # 5. Size Ratios
    df['feat_radius_ratio'] = df['r_M_angstrom'] / df['r_X_angstrom']
    
    # --- SELECT FEATURES ---
    feature_cols = [
        'feat_lattice_proxy',
        'feat_hydration_sum',
        'feat_enthalpy_sol_proxy',
        'feat_G_mismatch',
        'feat_radius_ratio',
        'r_M_angstrom',
        'r_X_angstrom',
        'cation_1_delta_G_hydration',
        'anion_1_delta_G_hydration'
    ]
    
    X = df[feature_cols]
    y = df[TARGET_COL]
    
    print(f"Engineered {len(feature_cols)} physical features.")
    print(f"Target Range: {y.min()} to {y.max()}")
    return X, y

def run_symbolic_regression():
    X, y = load_and_engineer_features(DATA_FILE)
    
    print("\nInitializing PySR Regressor (Robust Mode)...")
    print(f"Loss Function: {LOSS_FUNCTION} (minimizes absolute error to handle outliers)")
    
    model = PySRRegressor(
        niterations=N_ITERATIONS,
        populations=POPULATIONS,
        binary_operators=BINARY_OPERATORS,
        unary_operators=UNARY_OPERATORS,
        maxsize=MAX_EQUATION_LENGTH,      # Allow longer equations
        elementwise_loss=LOSS_FUNCTION,   # KEY FIX for outliers
        model_selection="best",
        temp_equation_file=os.path.join(OUTPUT_DIR, "hall_of_fame.csv"),
        extra_sympy_mappings={
            "inv": lambda x: 1/x
        }
    )
    
    model.fit(X, y)
    
    # --- RESULTS & SAVING ---
    final_output_path = os.path.join(OUTPUT_DIR, "pysr_equations_robust.csv")
    model.equations_.to_csv(final_output_path)
    print(f"\nEquations saved to: {final_output_path}")
    
    best_eqn = model.sympy()
    print(f"\nBEST EQUATION: {best_eqn}")

    # --- VISUALIZATION ---
    print("\nGenerating visualization...")
    y_pred = model.predict(X)
    
    # Identify dominant feature for plotting
    correlations = X.corrwith(y).abs()
    best_feat_name = correlations.idxmax()
    print(f"Plotting against dominant feature: {best_feat_name}")
    
    plt.figure(figsize=(10, 6))
    
    # Plot Observed
    plt.scatter(X[best_feat_name], y, color='blue', alpha=0.4, s=60, label='Observed')
    
    # Plot Predicted
    plt.scatter(X[best_feat_name], y_pred, color='red', alpha=0.6, marker='x', s=60, label='Predicted')
    
    plt.xlabel(best_feat_name)
    plt.ylabel(TARGET_COL)
    
    # Safe Title String (No Latex)
    eqn_str = str(best_eqn)
    if len(eqn_str) > 70:
        eqn_str = eqn_str[:35] + " ... " + eqn_str[-35:]
        
    plt.title(f"Robust Fit (L1 Loss)\nEq: {eqn_str}", fontsize=10)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plot_path = os.path.join(OUTPUT_DIR, "robust_fit_plot.png")
    plt.savefig(plot_path, dpi=150)
    print(f"Plot saved to: {plot_path}")

if __name__ == "__main__":
    run_symbolic_regression()