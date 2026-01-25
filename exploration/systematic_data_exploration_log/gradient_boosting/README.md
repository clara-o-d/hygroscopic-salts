# Gradient Boosting Model for ln(γ) at 90% RH

This directory contains a complete implementation of a gradient boosting regression model to predict the natural logarithm of the activity coefficient (ln(γ)) at 90% relative humidity for electrolyte solutions.

## Overview

**Objective:** Predict ln(γ) at 90% RH using ion properties and Pitzer parameters

**Algorithm:** Custom gradient boosting implementation with regression trees

**Dataset:** 108 electrolyte samples from `baseline_numeric_only.csv`

**Date:** January 25, 2026

---

## Files

### Scripts
- **`gradient_boosting.m`** - Main Octave/MATLAB script for model training
  - Loads and preprocesses data
  - Implements custom gradient boosting algorithm
  - Performs train/test split and evaluation
  - Computes feature importance
  - Saves results

- **`create_visualizations.py`** - Python script for generating plots
  - Predicted vs Observed plots
  - Residual analysis
  - Feature importance visualization
  - Error analysis by electrolyte

### Output Files (in `../../../figures/gradient_boosting/`)

1. **`gradient_boosting_results.txt`** - Detailed numerical results
   - Model hyperparameters
   - Performance metrics (RMSE, MAE, R²)
   - Feature importances
   - Individual predictions for test set

2. **`predictions.csv`** - All predictions with residuals
   - Columns: Electrolyte, Observed, Predicted, Residual, Dataset
   - 86 training samples + 22 test samples

3. **`RESULTS_SUMMARY.md`** - Comprehensive analysis summary
   - Executive summary of findings
   - Feature importance interpretation
   - Model performance analysis
   - Recommendations for improvement

4. **Visualizations:**
   - `predicted_vs_observed.png` - Scatter plots for train/test sets
   - `residuals_analysis.png` - Residual plots and distributions
   - `feature_importance.png` - Top 15 most important features
   - `error_analysis.png` - Worst predictions and error patterns

---

## Model Architecture

### Hyperparameters

```
Number of trees:     100
Learning rate:       0.1
Max tree depth:      4
Min leaf size:       5
Subsample ratio:     0.8
Train/test split:    80/20
Random seed:         42
```

### Algorithm Details

The implementation uses **sequential gradient boosting**:

1. **Initialize:** F₀(x) = mean(y)
2. **For each iteration m = 1 to M:**
   - Compute residuals: rᵢ = yᵢ - F_{m-1}(xᵢ)
   - Subsample 80% of training data
   - Fit regression tree to residuals
   - Update predictions: F_m(x) = F_{m-1}(x) + η * tree_m(x)
3. **Final prediction:** F_M(x)

**Loss function:** L2 (squared error)

**Split criterion:** Variance reduction

---

## Results Summary

### Performance Metrics

| Dataset  | RMSE    | MAE     | R²      |
|----------|---------|---------|---------|
| Training | 0.01182 | 0.00602 | 0.399   |
| Test     | 0.01460 | 0.00957 | -0.517  |

### Key Findings

1. **Overfitting:** The model shows severe overfitting with negative test R², indicating poor generalization.

2. **Feature Importance:** 
   - **Cation molecular weight** (26%) and **B_MX_0** Pitzer parameter (24%) are dominant
   - Cation properties generally more important than anion properties
   - Several features have zero importance (redundant)

3. **Error Patterns:**
   - Largest errors for zinc salts (ZnCl₂, ZnI₂)
   - Model struggles with high ln(γ) values
   - Bias toward predicting near-zero values

4. **Learning Curves:**
   - Training error decreases throughout
   - Test error increases after iteration 1-10
   - **Early stopping recommended** at ~10 iterations

---

## Top 5 Most Important Features

1. **cation_1_molecular_weight** (25.9%)
2. **B_MX_0_original** (23.7%) - Pitzer parameter
3. **r_M_angstrom** (8.1%) - Cation radius
4. **max_molality_original** (7.7%)
5. **C_MX_phi_original** (7.2%) - Pitzer parameter

---

## Worst Predictions (Test Set)

| Electrolyte | Observed | Predicted | Error    |
|-------------|----------|-----------|----------|
| ZnCl₂       | 0.0439   | -0.0093   | 0.0532   |
| LiClO₄      | 0.0136   | -0.0108   | 0.0244   |
| ZnI₂        | -0.0269  | -0.0082   | -0.0187  |
| CsBr        | 0.0195   | 0.0042    | 0.0153   |
| NH₄Cl       | 0.0019   | -0.0094   | 0.0112   |

---

## How to Run

### Train the Model (Octave/MATLAB)

```bash
cd exploration/systematic_data_exploration_log/gradient_boosting/
octave --no-gui gradient_boosting.m
```

**Requirements:**
- Octave 6.0+ or MATLAB R2018b+
- Statistics package: `pkg load statistics` (Octave)

**Output:** Results and predictions saved to `figures/gradient_boosting/`

### Create Visualizations (Python)

```bash
cd exploration/systematic_data_exploration_log/gradient_boosting/
python3 create_visualizations.py
```

**Requirements:**
- Python 3.7+
- pandas, numpy, matplotlib

**Output:** PNG images saved to `figures/gradient_boosting/`

---

## Recommendations for Improvement

### 1. Reduce Overfitting
- **Early stopping:** Stop at iteration 5-10 when test error minimizes
- **Increase regularization:** Lower learning rate (0.01-0.05)
- **Simplify model:** Reduce tree depth, increase min_leaf_size
- **Cross-validation:** Use k-fold CV for hyperparameter tuning

### 2. Feature Engineering
- Create interaction features (cation × anion properties)
- Add charge-to-radius ratios (charge density)
- Include polynomial terms for key features
- Remove zero-importance features

### 3. Data Strategy
- Collect more samples (especially divalent cations, transition metals)
- Use data augmentation or synthetic samples
- Consider multi-task learning (predict at multiple RH values)

### 4. Alternative Models
- **Random Forest:** May generalize better with less overfitting
- **Regularized Linear Models:** Ridge, Lasso, Elastic Net
- **Ensemble Methods:** Combine multiple model types
- **Neural Networks:** Physics-informed architectures

### 5. Physics-Based Validation
- Check thermodynamic consistency
- Compare with Pitzer model predictions
- Investigate chemical reasons for outliers (e.g., Zn²⁺ hydrolysis)

---

## Dataset Details

**Source:** `data/baseline_numeric_only.csv`

**Features (22 total):**

**Ion Properties:**
- Radii: `r_M_angstrom`, `r_X_angstrom`
- Hydration: `*_delta_G_hydration`, `*_radius_hydrated`
- Molecular: `*_molecular_weight`, `*_n_atoms`
- Transport: `*_diffusion_coefficient`, `*_viscosity_jones_dole`

**Pitzer Parameters:**
- `B_MX_0_original`, `B_MX_1_original`, `B_MX_2_original`
- `C_MX_phi_original`

**Other:**
- `max_molality_original`
- `electrolyte_type_numeric`
- `cation_type_numeric`, `anion_type_numeric`

**Target:** `ln_gamma_at_90RH`

**Samples:** 108 (after filtering for non-missing targets)

---

## Interpretation of Results

### Why the Model Overfits

1. **Small dataset:** 108 samples with 22 features (ratio ~5:1)
2. **High complexity:** 100 trees of depth 4 = too many parameters
3. **No regularization:** Pure boosting without constraints
4. **Noise sensitivity:** Small data allows memorization of noise

### Physical Insights

1. **Pitzer parameters matter:** B_MX_0 is highly important, validating that these parameters capture ion-ion interactions

2. **Cation dominance:** Cation molecular weight and radius are top features, suggesting cation hydration shell plays key role

3. **Size effects:** Ion radii (r_M) matter more than charge indicators, suggesting steric effects are important

4. **Missing physics:** Zero importance for some features suggests either redundancy or that model hasn't learned underlying physics

### What Works
- Model captures ~40% of variance on training data
- Feature importance aligns with physical intuition
- Training RMSE of 0.012 is reasonable given target range

### What Doesn't Work
- Negative test R² means worse than predicting mean
- Model doesn't generalize to new salts
- Likely memorizing specific salt patterns rather than learning ion interactions

---

## Citation

If you use this code or analysis, please cite:

```
Gradient Boosting Model for Activity Coefficients
AWH ML Discussion Project
Date: January 25, 2026
```

---

## Contact

For questions or issues, please contact the project maintainers.

---

*Last updated: January 25, 2026*
