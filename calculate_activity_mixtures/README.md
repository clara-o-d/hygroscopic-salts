# Calculate Activity Mixtures

This directory contains MATLAB functions to calculate water activity for binary salt mixtures based on experimental data.

## Overview

Each function calculates water activity (`aw`) as a function of the mass fractions of two salt components. The functions use polynomial fits (either univariate or bivariate) trained on experimental data from `mixture_sources_list.xlsx`.

## Available Mixtures

1. **NaCl + LiCl** - Sodium chloride and lithium chloride
2. **NaNO3 + LiNO3** - Sodium nitrate and lithium nitrate
3. **NaC2H3O2 + LiC2H3O2** - Sodium acetate and lithium acetate
4. **LiCl + BaCl2** - Lithium chloride and barium chloride
5. **LiCl + MgCl2** - Lithium chloride and magnesium chloride
6. **LiCl + CaCl2** - Lithium chloride and calcium chloride
7. **HClO4 + LiClO4** - Perchloric acid and lithium perchlorate
8. **HClO4 + NaClO4** - Perchloric acid and sodium perchlorate
9. **LiClO4 + NaClO4** - Lithium perchlorate and sodium perchlorate
10. **HCl + NaCl** - Hydrochloric acid and sodium chloride
11. **HCl + BaCl2** - Hydrochloric acid and barium chloride
12. **HCl + NaClO4** - Hydrochloric acid and sodium perchlorate
13. **HCl + Ba(ClO4)2** - Hydrochloric acid and barium perchlorate
14. **HCl + Na2SO4** - Hydrochloric acid and sodium sulfate
15. **LiCl + KCl** - Lithium chloride and potassium chloride
16. **LiCl + CsCl** - Lithium chloride and cesium chloride
17. **NaCl + CsCl** - Sodium chloride and cesium chloride
18. **MgCl2 + NaCl** - Magnesium chloride and sodium chloride
19. **MgCl2 + CaCl2** - Magnesium chloride and calcium chloride
20. **NH4Cl + LiCl** - Ammonium chloride and lithium chloride

## Usage

Each function takes two inputs: the mass fractions of the two components.

```matlab
% Example: Calculate water activity for NaCl + LiCl mixture
mf_NaCl = 0.10;  % 10% NaCl by mass
mf_LiCl = 0.05;  % 5% LiCl by mass

aw = calculate_activity_NaCl_LiCl(mf_NaCl, mf_LiCl);
fprintf('Water activity: %.4f\n', aw);
```

### Converting Molality to Mass Fraction

If you have molality values (mol/kg H₂O), convert to mass fractions:

```matlab
% Given molalities
m_NaCl = 2.0;  % mol/kg H2O
m_LiCl = 1.5;  % mol/kg H2O

% Molecular weights
MW_NaCl = 58.443;  % g/mol
MW_LiCl = 42.394;  % g/mol

% Calculate masses
mass_water = 1000;  % g (by definition of molality)
mass_NaCl = m_NaCl * MW_NaCl;
mass_LiCl = m_LiCl * MW_LiCl;
total_mass = mass_water + mass_NaCl + mass_LiCl;

% Calculate mass fractions
mf_NaCl = mass_NaCl / total_mass;
mf_LiCl = mass_LiCl / total_mass;

% Calculate water activity
aw = calculate_activity_NaCl_LiCl(mf_NaCl, mf_LiCl);
```

## Fit Quality

All functions use bivariate polynomial fits (degree 3) that capture the interaction between both components. The root-mean-square errors (RMSE) for each mixture are:

| Mixture | RMSE | Data Points |
|---------|------|-------------|
| NaCl + LiCl | 0.000135 | 36 |
| NaNO3 + LiNO3 | 0.004313 | 26 |
| NaC2H3O2 + LiC2H3O2 | 0.003989 | 26 |
| LiCl + BaCl2 | 0.003033 | 28 |
| LiCl + MgCl2 | 0.002975 | 34 |
| LiCl + CaCl2 | 0.006177 | 43 |
| HClO4 + LiClO4 | 0.009408 | 15 |
| HClO4 + NaClO4 | 0.004765 | 42 |
| LiClO4 + NaClO4 | 0.002810 | 33 |
| HCl + NaCl | 0.003699 | 10 |
| HCl + BaCl2 | 0.002157 | 12 |
| HCl + NaClO4 | 0.000293 | 7 |
| HCl + Ba(ClO4)2 | 0.003412 | 12 |
| HCl + Na2SO4 | 0.000218 | 13 |
| LiCl + KCl | 0.002768 | 45 |
| LiCl + CsCl | 0.009326 | 43 |
| NaCl + CsCl | 0.006851 | 40 |
| MgCl2 + NaCl | 0.002177 | 32 |
| MgCl2 + CaCl2 | 0.015789 | 36 |
| NH4Cl + LiCl | 0.002840 | 21 |

## Calibrated Ranges

Each function includes warnings when inputs fall outside the calibrated range. The calibrated ranges are documented in each function's header. For example, `calculate_activity_NaCl_LiCl.m` is calibrated for:
- NaCl: 0.0057 to 0.2366 mass fraction
- LiCl: 0.0047 to 0.1549 mass fraction
- Water activity: 0.7604 to 0.9826

## Generation

These functions were automatically generated using the Python script `../data/generate_mixture_activity_files.py`, which:
1. Reads experimental data from `mixture_sources_list.xlsx`
2. Fits bivariate polynomials to the water activity data
3. Generates MATLAB functions with appropriate documentation and validation

## Visualization

For 3D visualization of mixture water activity, see:
- `../multi_salt/visualize_nacl_licl_activity_3d.m` - 3D surface plots and cross-sections for NaCl + LiCl

## References

The experimental data sources are documented in `mixture_sources_list.xlsx`.
