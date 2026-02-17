"""
Lattice Occupation Model for LiCl Water Activity Coefficient
=============================================================
Uses a simple lattice/regular-solution heuristic for the energetics of
aqueous LiCl at different temperatures. Compares the heuristic activity
coefficient of water (gamma_w) against the actual gamma_w from experimental
data, at various water mole fractions (x_w).

Ion properties (radius, valence) are loaded from project data files.
Reference: Shannon (1976) ionic radii; baseline_numeric_only.csv.

Author: Generated for hygroscopic-salts project
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

# --- Paths ---
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
FIG_DIR = PROJECT_ROOT / "figures" / "simulations"

# --- Constants ---
R_GAS = 8.314462618  # J/(mol·K)
MW_WATER = 18.015  # g/mol
MW_LICL = 42.4  # g/mol


def load_ion_data():
    """Load ion radius and valence from project data files."""
    baseline_path = DATA_DIR / "baseline_numeric_only.csv"
    if baseline_path.exists():
        df = pd.read_csv(baseline_path)
        licl = df[df["electrolyte"] == "LiCl"].iloc[0]
        r_cat_A = float(licl["r_M_angstrom"])  # Li+ radius in Angstroms
        r_an_A = float(licl["r_X_angstrom"])   # Cl- radius in Angstroms
    else:
        # Fallback: Shannon (1976) 6-coordinate radii in Angstroms
        r_cat_A = 0.76  # Li+
        r_an_A = 1.81   # Cl-
    z_cat = 1
    z_an = 1
    return {"r_cat": r_cat_A, "r_an": r_an_A, "z_cat": z_cat, "z_an": z_an}


def lattice_omega(ion_data, T_K, alpha=-5000.0):
    """
    Interaction parameter omega (J/mol) for regular-solution lattice model.
    omega = alpha * (1 + charge_density_factor): heuristic scaling by ion properties.
    Charge density = z/r^3; higher charge density -> stronger hydration -> gamma_w < 1.
    omega < 0 gives ln(gamma_w) < 0 (water activity depressed below ideal).
    """
    r_cat = ion_data["r_cat"]  # Angstroms
    r_an = ion_data["r_an"]
    z_cat = ion_data["z_cat"]
    z_an = ion_data["z_an"]
    rho_cat = z_cat / (r_cat ** 3)
    rho_an = z_an / (r_an ** 3)
    mean_rho = 0.5 * (rho_cat + rho_an)
    # Normalize: LiCl has mean_rho ~ 1.6; use as scaling factor
    charge_factor = mean_rho / 1.6  # ~1 for LiCl
    omega = alpha * charge_factor
    return omega


def gamma_w_lattice(x_w, ion_data, T_K, alpha=-5000.0):
    """
    Activity coefficient of water from Bragg-Williams / regular solution model.
    ln(gamma_w) = (omega / RT) * (1 - x_w)^2
    omega < 0 -> gamma_w < 1 (water activity depressed by ions).
    """
    omega = lattice_omega(ion_data, T_K, alpha=alpha)
    x_s = 1.0 - x_w
    ln_gamma = (omega / (R_GAS * T_K)) * (x_s ** 2)
    return np.exp(ln_gamma)


def calculate_mf_licl(RH, T_C):
    """
    Mass fraction of LiCl at given RH and T (C).
    From: https://doi.org/10.1016/j.ijthermalsci.2003.09.003
    Solves for xi (mass fraction salt) such that RH = f(xi, T).
    """
    if RH <= 0 or RH >= 1:
        raise ValueError("RH must be in (0, 1)")
    if T_C > 100:
        raise ValueError("T should be in Celsius, < 100")

    p_0, p_1, p_2 = 0.28, 4.3, 0.60
    p_3, p_4, p_5 = 0.21, 5.10, 0.49
    p_6, p_7, p_8, p_9 = 0.362, -4.75, -0.40, 0.03

    theta = (T_C + 273.15) / 647.0  # reduced temperature

    def residual(xi):
        term1 = (1 + (xi / p_6) ** p_7) ** p_8
        term2 = p_9 * np.exp(-((xi - 0.1) ** 2) / 0.005)
        bracket1 = 1 - term1 - term2
        term3 = (1 + (xi / p_0) ** p_1) ** p_2
        term4 = ((1 + (xi / p_3) ** p_4) ** p_5 - 1) * theta
        bracket2 = 2 - term3 + term4
        return RH - bracket1 * bracket2

    try:
        xi = brentq(residual, 0.01, 0.75)
    except ValueError:
        xi = np.nan
    return xi


def mass_to_mole_fraction_water(mf_salt):
    """Convert mass fraction salt to mole fraction water."""
    mf_water = 1.0 - mf_salt
    n_water = mf_water / MW_WATER
    n_salt = mf_salt / MW_LICL
    x_w = n_water / (n_water + n_salt)
    return x_w


def get_actual_licl_data():
    """Load actual LiCl activity coefficient data from water_activity CSV."""
    csv_path = DATA_DIR / "water_activity_all_salts_combined.csv"
    df = pd.read_csv(csv_path)
    licl = df[df["Salt"] == "LiCl"].copy()
    licl = licl.rename(columns={
        "Mole_Fraction_Water": "x_w",
        "Activity_Coefficient_Water": "gamma_w_actual",
        "Temperature_C": "T_C",
    })
    return licl[["T_C", "x_w", "gamma_w_actual"]]


def generate_licl_at_temperatures(T_C_list, n_rh=80):
    """
    Generate (x_w, gamma_w) at specified temperatures using calculate_mf_LiCl.
    gamma_w = a_w / x_w = RH / x_w
    """
    rh_vec = np.linspace(0.15, 0.95, n_rh)
    rows = []
    for T_C in T_C_list:
        for RH in rh_vec:
            try:
                mf_salt = calculate_mf_licl(RH, T_C)
                x_w = mass_to_mole_fraction_water(mf_salt)
                gamma_w = RH / x_w if x_w > 0 else np.nan
                rows.append({"T_C": T_C, "x_w": x_w, "gamma_w_actual": gamma_w, "RH": RH})
            except (ValueError, Exception):
                continue
    return pd.DataFrame(rows)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ion_data = load_ion_data()
    print("Ion data (LiCl):", ion_data)

    # --- Load actual data ---
    actual_df = get_actual_licl_data()
    T_actual = actual_df["T_C"].unique()
    print(f"Actual LiCl data: {len(actual_df)} points at T = {T_actual} °C")

    # --- Generate data at multiple temperatures ---
    T_list = [10, 25, 40]
    gen_df = generate_licl_at_temperatures(T_list)
    print(f"Generated LiCl data: {len(gen_df)} points at T = {T_list} °C")

    # Combine: use actual for 25°C (denser), generated for 10 and 40°C
    actual_25 = actual_df[actual_df["T_C"] == 25].copy()
    actual_25["RH"] = np.nan
    gen_other = gen_df[gen_df["T_C"] != 25]
    combined = pd.concat([actual_25[["T_C", "x_w", "gamma_w_actual", "RH"]], gen_other], ignore_index=True)

    # --- Lattice model: fit alpha to 25°C data ---
    ref = combined[combined["T_C"] == 25].dropna()
    if len(ref) < 5:
        ref = actual_df.dropna()
    x_ref = ref["x_w"].values
    gamma_ref = ref["gamma_w_actual"].values
    T_ref = 298.15  # K

    # Simple least-squares fit for alpha (negative for gamma_w < 1)
    def err(alpha):
        pred = np.array([gamma_w_lattice(x, ion_data, T_ref, alpha=alpha) for x in x_ref])
        return np.mean((np.log(pred + 1e-10) - np.log(gamma_ref + 1e-10)) ** 2)

    from scipy.optimize import minimize_scalar
    res = minimize_scalar(err, bounds=(-5e4, -5e2), method="bounded")
    alpha_fit = res.x
    print(f"Fitted alpha = {alpha_fit:.3e}")

    # --- Plot: Lattice vs Actual at different x_w and T ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: gamma_w vs x_w at different temperatures
    ax1 = axes[0]
    colors = {10: "#2166ac", 25: "#4393c3", 40: "#92c5de"}
    for T_C in sorted(combined["T_C"].unique()):
        sub = combined[combined["T_C"] == T_C].sort_values("x_w")
        x_w = sub["x_w"].values
        gamma_act = sub["gamma_w_actual"].values
        T_K = T_C + 273.15
        gamma_lat = np.array([gamma_w_lattice(x, ion_data, T_K, alpha=alpha_fit) for x in x_w])
        ax1.plot(x_w, gamma_act, "o", color=colors.get(T_C, "gray"), ms=4, label=f"Actual, T={T_C}°C")
        ax1.plot(x_w, gamma_lat, "-", color=colors.get(T_C, "gray"), lw=2, label=f"Lattice, T={T_C}°C")

    ax1.set_xlabel(r"Water mole fraction $x_w$", fontsize=12)
    ax1.set_ylabel(r"Activity coefficient of water $\gamma_w$", fontsize=12)
    ax1.set_title("LiCl: Lattice Occupation Model vs Actual", fontsize=14)
    ax1.legend(loc="best", fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0.5, 1.0)

    # Right: Parity plot (lattice vs actual)
    ax2 = axes[1]
    for T_C in sorted(combined["T_C"].unique()):
        sub = combined[combined["T_C"] == T_C].dropna()
        x_w = sub["x_w"].values
        gamma_act = sub["gamma_w_actual"].values
        T_K = T_C + 273.15
        gamma_lat = np.array([gamma_w_lattice(x, ion_data, T_K, alpha=alpha_fit) for x in x_w])
        ax2.scatter(gamma_act, gamma_lat, c=colors.get(T_C, "gray"), s=30, label=f"T={T_C}°C", alpha=0.7)

    lims = [0.1, 1.0]
    ax2.plot(lims, lims, "k--", lw=2, label="1:1")
    ax2.set_xlabel(r"Actual $\gamma_w$", fontsize=12)
    ax2.set_ylabel(r"Lattice model $\gamma_w$", fontsize=12)
    ax2.set_title("Parity: Lattice vs Actual", fontsize=14)
    ax2.legend(loc="best", fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_aspect("equal")

    plt.tight_layout()
    out_path = FIG_DIR / "lattice_occupation_licl_activity.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")

    # --- Summary statistics ---
    all_act = []
    all_lat = []
    for _, row in combined.dropna().iterrows():
        T_K = row["T_C"] + 273.15
        g_act = row["gamma_w_actual"]
        g_lat = gamma_w_lattice(row["x_w"], ion_data, T_K, alpha=alpha_fit)
        all_act.append(g_act)
        all_lat.append(g_lat)
    all_act = np.array(all_act)
    all_lat = np.array(all_lat)
    rmse = np.sqrt(np.mean((all_lat - all_act) ** 2))
    mae = np.mean(np.abs(all_lat - all_act))
    print(f"RMSE (gamma_w): {rmse:.4f}")
    print(f"MAE (gamma_w):  {mae:.4f}")


if __name__ == "__main__":
    main()
