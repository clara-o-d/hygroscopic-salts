"""
Master Equation Symbolic Regression (High-Volume Visualization)
---------------------------------------------------------------
1. PHYSICS: Uses 'feat_DH_shape' to enforce correct limiting behavior.
2. ROBUST: Auto-detects ln(gamma) vs gamma.
3. VIZ: Generates a 16-salt grid AND a global parity plot.

Author: Generated for AWH ML Discussion
Date: 2026-01-28
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pysr import PySRRegressor
from sklearn.metrics import r2_score, mean_squared_error

# --- CONFIGURATION ---
DATA_FILE = './master_long_format_for_sr.csv'
OUTPUT_DIR = '../../../figures/symbolic_regression/'

# Search Settings
N_ITERATIONS = 100
POPULATIONS = 40
MAX_EQUATION_LENGTH = 50
LOSS_FUNCTION = "L1DistLoss()" # Robust MAE loss

BINARY_OPERATORS = ["+", "-", "*", "/", "^"]
UNARY_OPERATORS = ["exp", "abs", "sqrt", "log", "inv(x) = 1/x"]

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_and_engineer_features(filepath):
    print(f"Loading Master Dataset from {filepath}...")
    df = pd.read_csv(filepath)
    
    # --- 1. TARGET HANDLING ---
    # Try to find the target column
    possible_targets = ['ln_gamma', 'ln_gammaw', 'Activity_Coefficient', 'gamma']
    target_col = next((c for c in possible_targets if c in df.columns), None)
    
    if target_col is None:
        raise ValueError(f"Target not found. Columns: {df.columns.tolist()}")

    # Intelligent transformation
    sample_mean = df[target_col].mean()
    if sample_mean < 0:
        print(f"  -> Target is negative (avg {sample_mean:.2f}). Assuming already ln(gamma).")
        df['y_final'] = df[target_col]
    else:
        print(f"  -> Target is positive (avg {sample_mean:.2f}). Applying log transform.")
        df = df[df[target_col] > 1e-9] # Safety clip
        df['y_final'] = np.log(df[target_col])

    # Clean extreme outliers
    df = df[(df['y_final'] > -6.0) & (df['y_final'] < 2.0)]
    
    # --- 2. DYNAMIC FEATURES ---
    if 'Ionic_Strength' not in df.columns:
        df['Ionic_Strength'] = df['Molality_m'] if 'Molality_m' in df.columns else np.nan
        if df['Ionic_Strength'].isna().all():
            raise ValueError("Need 'Ionic_Strength' or 'Molality_m'")

    df['sqrt_I'] = np.sqrt(df['Ionic_Strength'])
    
    # Physics Helper: Debye-Huckel Shape Term
    # This helps the model find the curvature at low concentrations
    df['feat_DH_shape'] = df['sqrt_I'] / (1.0 + 1.5 * df['sqrt_I'])
    
    # --- 3. STATIC FEATURES ---
    df['inv_sum_radii'] = 1.0 / (df['r_M_angstrom'] + df['r_X_angstrom'])
    df['feat_lattice_proxy'] = df['inv_sum_radii']
    df['feat_hydration_sum'] = df['cation_1_delta_G_hydration'] + df['anion_1_delta_G_hydration']
    df['feat_enthalpy_sol_proxy'] = df['feat_lattice_proxy'] + (df['feat_hydration_sum'] / 1000.0) 
    df['feat_G_mismatch'] = (df['cation_1_delta_G_hydration'] - df['anion_1_delta_G_hydration']).abs()
    
    # --- SELECT FEATURES ---
    feature_cols = [
        'feat_DH_shape', 'sqrt_I',          # Dynamic
        'feat_lattice_proxy', 'feat_hydration_sum', 
        'feat_enthalpy_sol_proxy', 'feat_G_mismatch',
        'r_M_angstrom', 'r_X_angstrom',
        'cation_1_delta_G_hydration', 'anion_1_delta_G_hydration'
    ]
    
    X = df[feature_cols]
    y = df['y_final']
    
    print(f"Training on {len(df)} points. Target Range: {y.min():.2f} to {y.max():.2f}")
    return X, y, df

def run_symbolic_regression():
    X, y, df_full = load_and_engineer_features(DATA_FILE)
    
    print("\nInitializing PySR...")
    model = PySRRegressor(
        niterations=N_ITERATIONS,
        populations=POPULATIONS,
        binary_operators=BINARY_OPERATORS,
        unary_operators=UNARY_OPERATORS,
        maxsize=MAX_EQUATION_LENGTH,
        elementwise_loss=LOSS_FUNCTION, 
        model_selection="best",
        parsimony=0.0001,
        temp_equation_file=os.path.join(OUTPUT_DIR, "hall_of_fame.csv"),
        extra_sympy_mappings={"inv": lambda x: 1/x}
    )
    
    model.fit(X, y)
    
    # Save Equations
    final_output_path = os.path.join(OUTPUT_DIR, "pysr_equations_master.csv")
    model.equations_.to_csv(final_output_path)
    print(f"\nEquations saved to: {final_output_path}")
    
    best_eqn = model.sympy()
    print(f"\nBEST EQUATION: {best_eqn}")

    # --- VISUALIZATION SUITE ---
    print("\nGenerating visualizations...")
    y_pred = model.predict(X)
    df_full['y_pred'] = y_pred
    
    # 1. GLOBAL PARITY PLOT (The "Truth" Plot)
    plt.figure(figsize=(8, 8))
    plt.scatter(y, y_pred, alpha=0.3, s=15, color='blue', label='Data Points')
    
    # Reference Line
    min_val, max_val = min(y.min(), y_pred.min()), max(y.max(), y_pred.max())
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2, label='Perfect Fit')
    
    r2 = r2_score(y, y_pred)
    mse = mean_squared_error(y, y_pred)
    
    plt.title(f"Global Parity Plot\n$R^2$ = {r2:.3f} | MSE = {mse:.4f}")
    plt.xlabel("Observed ln($\\gamma$)")
    plt.ylabel("Predicted ln($\\gamma$)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(OUTPUT_DIR, "global_parity_plot.png"), dpi=150)
    plt.close()

    # 2. 16-SALT GRID (Curve Check)
    unique_salts = df_full['electrolyte'].unique()
    # Pick 16 salts (or all if <16)
    n_plot = min(16, len(unique_salts))
    plot_salts = np.random.choice(unique_salts, size=n_plot, replace=False)
    
    fig, axes = plt.subplots(4, 4, figsize=(16, 12))
    axes = axes.flatten()
    
    for i, salt in enumerate(plot_salts):
        ax = axes[i]
        subset = df_full[df_full['electrolyte'] == salt].sort_values('sqrt_I')
        
        ax.scatter(subset['sqrt_I'], subset['y_final'], color='black', s=20, alpha=0.6)
        ax.plot(subset['sqrt_I'], subset['y_pred'], color='red', lw=2)
        
        ax.set_title(salt, fontsize=10)
        if i >= 12: ax.set_xlabel("sqrt(I)")
        if i % 4 == 0: ax.set_ylabel("ln($\\gamma$)")
        ax.grid(True, alpha=0.2)
        
    plt.suptitle(f"Master Equation Fits (Random 16 Salts)", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "salt_grid_plot.png"), dpi=150)
    plt.close()
    
    print(f"All plots saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    run_symbolic_regression()