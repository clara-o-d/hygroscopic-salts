"""
Lattice Occupation Model for Water Activity Coefficient
========================================================
Regular-solution heuristic vs actual gamma_w. Uses MATLAB calculate_mf for full
RH range (low x_w regime). Ion data from baseline_numeric_only.csv.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
FIG_DIR = PROJECT_ROOT / "figures" / "simulations"
R_GAS, MW_WATER = 8.314462618, 18.015

# 12-15 salts: name, MW, RH_min, RH_max, T_C (from load_salt_data)
SALT_CONFIG = [
    ("LiCl", 42.4, 0.13, 0.96, 25), ("NaCl", 58.443, 0.77, 0.99, 25), ("KCl", 74.551, 0.86, 0.99, 25),
    ("MgCl2", 95.2, 0.34, 0.96, 25), ("CaCl2", 111, 0.32, 0.96, 25), ("LiBr", 86.85, 0.08, 0.96, 25),
    ("ZnCl2", 136.3, 0.08, 0.96, 25), ("LiI", 133.85, 0.19, 0.96, 25), ("NaBr", 102.89, 0.62, 0.93, 25),
    ("NaI", 149.89, 0.59, 0.96, 25), ("HCl", 36.5, 0.18, 0.96, 25), ("LiOH", 24, 0.86, 0.96, 25),
    ("NH4Cl", 53.491, 0.82, 0.99, 25), ("CsCl", 168.363, 0.83, 0.99, 25), ("Na2SO4", 142.04, 0.90, 0.99, 25),
]
N_RH = 120  # Dense sampling for low x_w


def load_ion_data(salt_name):
    df = pd.read_csv(DATA_DIR / "baseline_numeric_only.csv")
    row = df[df["electrolyte"] == salt_name]
    if row.empty:
        row = df[df["electrolyte_stripped"] == salt_name]
    if row.empty:
        raise ValueError(f"{salt_name} not in baseline")
    r = row.iloc[0]
    z_cat, z_an = (2, 1) if salt_name in ("MgCl2", "CaCl2", "ZnCl2") else (1, 2) if "SO4" in salt_name else (1, 1)
    return {"r_cat": float(r["r_M_angstrom"]), "r_an": float(r["r_X_angstrom"]), "z_cat": z_cat, "z_an": z_an}


def gamma_w_lattice(x_w, ion_data, T_K, alpha=-5000.0):
    r_cat, r_an = ion_data["r_cat"], ion_data["r_an"]
    rho = 0.5 * (ion_data["z_cat"] / r_cat**3 + ion_data["z_an"] / r_an**3)
    omega = alpha * (rho / 1.6)
    return np.exp((omega / (R_GAS * T_K)) * (1 - x_w)**2)


def call_matlab_batch(all_rows):
    """Call MATLAB once with all (salt, RH, T) rows. One MATLAB startup for all salts."""
    in_df = pd.DataFrame(all_rows)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        in_df.to_csv(f.name, index=False)
        in_path = f.name
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        out_path = f.name
    try:
        sim_dir = str(SCRIPT_DIR).replace("\\", "/")
        r = subprocess.run(["matlab", "-batch", f"cd('{sim_dir}'); run_calculate_mf_batch('{in_path}','{out_path}')"],
                           capture_output=True, text=True, timeout=120, cwd=str(PROJECT_ROOT))
        if r.returncode != 0:
            raise RuntimeError(r.stderr)
        return pd.read_csv(out_path)
    finally:
        Path(in_path).unlink(missing_ok=True)
        Path(out_path).unlink(missing_ok=True)


def mf_to_xw(mf, mw):
    nw, ns = (1 - mf) / MW_WATER, mf / mw
    return nw / (nw + ns)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    # Build one batch: all (salt, RH, T) for all salts
    all_rows = []
    mw_by_salt = {}
    for name, mw, rh_lo, rh_hi, t in SALT_CONFIG:
        try:
            load_ion_data(name)
        except (ValueError, FileNotFoundError):
            print(f"Skip {name}: no ion data")
            continue
        mw_by_salt[name] = mw
        rh = np.linspace(rh_lo, rh_hi, N_RH)
        for rv in rh:
            all_rows.append({"salt": name, "RH": rv, "T": t})
    all_data = []
    try:
        out = call_matlab_batch(all_rows)
    except (RuntimeError, FileNotFoundError) as e:
        print(f"MATLAB failed: {e}. Fallback to CSV.")
        out = None
    if out is not None:
        for name in mw_by_salt:
            sub = out[out["salt"] == name]
            mw = mw_by_salt[name]
            rows = []
            for _, r in sub.iterrows():
                mf = r["mf"]
                if np.isnan(mf) or mf <= 0 or mf >= 1:
                    continue
                xw = mf_to_xw(mf, mw)
                rows.append({"salt": name, "T_C": r["T"], "x_w": xw, "gamma_w_actual": r["RH"] / xw})
            if rows:
                all_data.append(pd.DataFrame(rows))
    if not all_data:
        print("No data. Fallback: use water_activity CSV.")
        df = pd.read_csv(DATA_DIR / "water_activity_all_salts_combined.csv")
        cfg = {c[0]: c for c in SALT_CONFIG}
        for name in cfg:
            sub = df[df["Salt"] == name]
            if sub.empty:
                continue
            try:
                load_ion_data(name)
            except (ValueError, FileNotFoundError):
                continue
            mw = cfg[name][1]
            sub = sub.rename(columns={"Mole_Fraction_Water": "x_w", "Activity_Coefficient_Water": "gamma_w_actual", "Temperature_C": "T_C"})
            sub["salt"] = name
            all_data.append(sub[["salt", "T_C", "x_w", "gamma_w_actual"]])
    combined = pd.concat(all_data, ignore_index=True)

    # Fit alpha per salt
    alpha_by_salt = {}
    for name, _, _, _, t in SALT_CONFIG:
        try:
            ion_data = load_ion_data(name)
        except (ValueError, FileNotFoundError):
            continue
        ref = combined[(combined["salt"] == name) & (combined["T_C"] == t)].dropna()
        if len(ref) < 5:
            ref = combined[combined["salt"] == name].dropna()
        if len(ref) < 5:
            alpha_by_salt[name] = -15000.0
            continue
        x_ref, g_ref = ref["x_w"].values, ref["gamma_w_actual"].values
        from scipy.optimize import minimize_scalar
        err = lambda a: np.mean((np.log(np.array([gamma_w_lattice(x, ion_data, t + 273.15, a) for x in x_ref]) + 1e-10) - np.log(g_ref + 1e-10))**2)
        alpha_by_salt[name] = minimize_scalar(err, bounds=(-5e4, -5e2), method="bounded").x

    # Separate subplot per salt (4x4 grid)
    n_salts = len(SALT_CONFIG)
    n_cols = 4
    n_rows = (n_salts + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 3.5 * n_rows))
    axes = axes.flatten()
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    for idx, (name, _, _, _, t) in enumerate(SALT_CONFIG):
        ax = axes[idx]
        try:
            ion_data = load_ion_data(name)
        except (ValueError, FileNotFoundError):
            ax.text(0.5, 0.5, f"{name}\n(no data)", ha="center", va="center"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            continue
        sub = combined[combined["salt"] == name].sort_values("x_w")
        if sub.empty:
            ax.text(0.5, 0.5, f"{name}\n(no data)", ha="center", va="center"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            continue
        xw, g_act = sub["x_w"].values, sub["gamma_w_actual"].values
        alpha = alpha_by_salt.get(name, -15000.0)
        g_lat = np.array([gamma_w_lattice(x, ion_data, t + 273.15, alpha) for x in xw])
        c = colors[idx % 10]
        ax.plot(xw, g_act, "o", color=c, ms=3, label="Actual"); ax.plot(xw, g_lat, "-", color=c, lw=2, label="Lattice")
        ax.set_xlabel(r"$x_w$"); ax.set_ylabel(r"$\gamma_w$"); ax.set_title(f"{name} (T={int(t)}°C)")
        ax.legend(loc="best", fontsize=8); ax.grid(True, alpha=0.3); ax.set_xlim(0.15, 1.0)
    for j in range(n_salts, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle("Regular Solution Model vs Actual (full calculate_mf range, low $x_w$)", fontsize=12)
    plt.tight_layout(); plt.savefig(FIG_DIR / "lattice_occupation_multi_salt_activity.png", dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved: {FIG_DIR / 'lattice_occupation_multi_salt_activity.png'}")

    # Parity plot
    fig2, ax2 = plt.subplots(1, 1, figsize=(7, 6))
    for idx, (name, _, _, _, t) in enumerate(SALT_CONFIG):
        try:
            ion_data = load_ion_data(name)
        except (ValueError, FileNotFoundError):
            continue
        sub = combined[combined["salt"] == name].dropna()
        if sub.empty:
            continue
        alpha = alpha_by_salt.get(name, -15000.0)
        g_act = sub["gamma_w_actual"].values
        g_lat = np.array([gamma_w_lattice(r["x_w"], ion_data, r["T_C"] + 273.15, alpha) for _, r in sub.iterrows()])
        ax2.scatter(g_act, g_lat, c=[colors[idx % 10]], s=15, label=name, alpha=0.7)
    ax2.plot([0.05, 1.05], [0.05, 1.05], "k--", lw=2)
    ax2.set_xlabel(r"Actual $\gamma_w$"); ax2.set_ylabel(r"Lattice $\gamma_w$"); ax2.set_title("Parity: Regular Solution vs Actual")
    ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=1, fontsize=8)
    ax2.grid(True, alpha=0.3); ax2.set_xlim(0.05, 1.05); ax2.set_ylim(0.05, 1.05); ax2.set_aspect("equal")
    plt.tight_layout(); plt.savefig(FIG_DIR / "lattice_occupation_multi_salt_parity.png", dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved: {FIG_DIR / 'lattice_occupation_multi_salt_parity.png'}")

    for name, _, _, _, t in SALT_CONFIG:
        sub = combined[combined["salt"] == name].dropna()
        if sub.empty:
            continue
        try:
            ion_data = load_ion_data(name)
        except (ValueError, FileNotFoundError):
            continue
        alpha = alpha_by_salt.get(name, -15000.0)
        act = sub["gamma_w_actual"].values
        lat = np.array([gamma_w_lattice(r["x_w"], ion_data, r["T_C"] + 273.15, alpha) for _, r in sub.iterrows()])
        print(f"  {name} RMSE: {np.sqrt(np.mean((lat - act)**2)):.4f}")


if __name__ == "__main__":
    main()
