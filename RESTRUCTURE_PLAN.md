# Repository Restructuring Plan

This document outlines suggested changes to make the hygroscopic-salts repository more concise, maintainable, and uniformly named.

**Status: Implemented (Feb 2025)** — Consolidation complete. Exploration folder left unchanged per user request. Baseline CSV not renamed (exploration depends on it).

---

## 1. Consolidate `isotherms_screening` into One File

### Current State
- **`Isotherms_screening_endothermic.m`**: Uses `load_salt_data()` to iterate over all salts; produces uptake curves, enthalpy vs DRH, solubility validation, and G_excess plots. Outputs to `figures/uptake/`.
- **`Isotherms_screening_exothermic.m`**: Manually hardcodes ~12 exothermic salts (LiCl, CaCl2, MgCl2, etc.) with individual loops; different output names (`Uptake_screening_exothermic_gg`, etc.).

### Recommendation
**Merge into a single script** `isotherms_screening_all_salts.m` that:
- Uses `load_salt_data()` as the single source of truth (already contains both endothermic and exothermic salts)
- Handles `func_args` (temperature) for salts like LiCl, CaCl2 that need it
- Produces uptake curves for all salts in one pass
- Optionally supports filtering by salt type (e.g., `'endothermic'`, `'exothermic'`, `'sulfates'`) via a parameter at the top for subset plots

The exothermic script’s manual approach is redundant because `load_salt_data()` already includes those salts with correct RH ranges and function names. The main fix is to call `calculate_mf_*` with temperature when `func_args == 1` (as the endothermic script does via `calc_func(rh_vec(j))`—but some functions need `T`; check `calculate_mf_LiCl(RH, T)`).

### Implementation Notes
- In the processing loop, use `r{6}` (func_args) to decide: `calc_func(rh)` vs `calc_func(rh, T)`. **Note:** The current endothermic script does *not* pass T—it always calls `calc_func(rh_vec(j))`. This is a bug for LiCl, CaCl2, etc. Use the pattern from `save_all_data_to_csv.m` and `plot_water_activity_all_salts.m` which correctly check `func_args` and pass T from `r{9}` when needed.
- Remove the duplicate `MgNO32` entry in `load_salt_data` if it's redundant with `MgNO3`
- Output figures with a consistent naming scheme (see Section 3)

---

## 2. Separate Solubility–Humidity Comparison Script

### Current State
The solubility validation (literature vs calculated from DRH) lives in `Isotherms_screening_endothermic.m` (lines 232–280), producing `Solubility_vs_Est_Solubility.png`.

### Recommendation
Create **`scripts/validate_solubility_vs_rh.m`** (or `analysis/validate_solubility_vs_rh.m`) that:
- Loads `load_salt_data()` and `baseline_numeric_only.csv`
- Computes estimated solubility from DRH for each salt
- Plots literature vs calculated solubility (1:1 comparison)
- Saves to `figures/uptake/` or `figures/validation/`

This script can be run standalone or called from the main isotherms script. Keeping it separate makes the validation analysis reusable and easier to extend (e.g., adding more reference data sources).

---

## 3. Uniform Naming Scheme

### Current Inconsistencies
| Category | Examples |
|----------|----------|
| MATLAB scripts | `Isotherms_screening_endothermic` (PascalCase) vs `plot_water_activity_all_salts` (snake_case) |
| calculate_mf | `calculate_mf_CaCl` (CaCl vs CaCl2) vs `calculate_mf_BaNO32` (inconsistent subscripts) |
| Data files | `baseline_numeric_only.csv` vs `water_activity_all_salts_combined.csv` |
| Figure outputs | `RH_vs_Uptake_gg.png` vs `Uptake_screening_exothermic_gg.png` |

### Proposed Convention

#### MATLAB Files
- **Scripts**: `snake_case` with descriptive action + domain  
  - `isotherms_screening_all_salts.m`  
  - `validate_solubility_vs_rh.m`  
  - `plot_water_activity_all_salts.m` (already good)
- **Functions**: `snake_case` or `camelCase` (MATLAB convention often uses camelCase for functions)  
  - `load_salt_data.m` (already good)  
  - `calculate_mf_NaCl.m` (keep for compatibility; these are called by name)

#### Data Files
- **Pattern**: `{domain}_{descriptor}.csv`  
  - `salt_baseline_properties.csv` (instead of `baseline_numeric_only.csv`)  
  - `water_activity_all_salts.csv` (instead of `water_activity_all_salts_combined.csv`)

#### Figure Outputs
- **Pattern**: `{category}_{quantity}_vs_{xaxis}.png`  
  - `uptake_rh_vs_gg.png`  
  - `uptake_rh_vs_molmol.png`  
  - `validation_solubility_lit_vs_calc.png`  
  - `uptake_enthalpy_vs_drh.png`  
  - `uptake_gmix_excess_vs_mole_fraction.png`

#### calculate_mf Functions
- Keep existing names to avoid breaking `load_salt_data` and other callers
- Document the naming rule: `calculate_mf_{SaltFormula}` where formula uses minimal subscripts (e.g., `CaCl` for CaCl2, `MgNO3` for Mg(NO3)2)

---

## 4. Suggested Directory Structure

```
hygroscopic-salts/
├── data/                    # All input data
│   ├── salt_baseline_properties.csv
│   ├── load_salt_data.m
│   └── ...
├── scripts/                 # Main analysis scripts (NEW or RENAME)
│   ├── isotherms_screening_all_salts.m
│   ├── validate_solubility_vs_rh.m
│   └── ...
├── calculate_mf/            # Salt-specific calculation functions (unchanged)
├── plot_water_activity/     # Water activity plotting
├── pitzer/                  # Pitzer model work
├── exploration/             # Exploratory analysis
├── figures/                 # Output figures by category
│   ├── uptake/
│   ├── validation/
│   └── ...
└── util/
```

Alternatively, keep `isotherms_screening/` as a folder but with a single file: `isotherms_screening/run_all_salts.m`.

---

## 5. Migration Checklist

- [x] Create `isotherms_screening_all_salts.m` (merge endothermic + exothermic logic)
- [x] Add temperature argument handling for `func_args == 1` salts
- [x] Create `validate_solubility_vs_rh.m` and remove that block from the main isotherms script
- [x] Update figure output names to new convention (uptake_*, validation_*)
- [ ] Rename `baseline_numeric_only.csv` — skipped (exploration folder depends on it)
- [ ] Gradually rename other MATLAB scripts to snake_case (low priority)
- [x] Remove `Isotherms_screening_endothermic.m` and `Isotherms_screening_exothermic.m`
- [x] Consolidate `plot_water_activity/subsets/` into `plot_water_activity_subset.m` (parameterized)

---

## 6. Additional Concision Ideas

- **plot_water_activity/subsets/**: Consider a single parameterized script `plot_water_activity_subset.m` that accepts a salt filter (e.g., `'sulfates'`, `'halides'`) instead of six separate files.
- **exploration/**: The `systematic_data_exploration_log/` subfolder has many small scripts; consider consolidating related analyses (e.g., PLS + gradient boosting + mixed effects) into a single `run_exploration_pipeline.m` with sections.
- **calculate_mf/**: Keep as-is; the 55+ files are necessary for the function-name dispatch pattern. No consolidation needed.
