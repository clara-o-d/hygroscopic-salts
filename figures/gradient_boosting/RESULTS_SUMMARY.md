# Gradient Boosting Model: ln(γ) at 90% RH

## Executive Summary

A gradient boosting regression model was successfully trained to predict the natural logarithm of the activity coefficient (ln(γ)) at 90% relative humidity for various electrolyte solutions using ion properties and Pitzer parameters as features.

**Date:** January 25, 2026

---

## Dataset Overview

- **Total samples:** 108 electrolytes with ln(γ) at 90% RH data
- **Features:** 22 physicochemical properties
  - Ion properties (radii, molecular weights, hydration energies, diffusion coefficients, viscosity coefficients)
  - Pitzer parameters (B_MX_0, B_MX_1, B_MX_2, C_MX_phi)
  - Electrolyte type indicators
- **Target range:** -0.0374 to 0.1079
- **Train/Test split:** 80/20 (86 training, 22 test samples)

---

## Model Configuration

### Hyperparameters
- **Number of trees:** 100
- **Learning rate:** 0.1 (shrinkage parameter)
- **Max tree depth:** 4
- **Min leaf size:** 5 observations
- **Subsample ratio:** 0.8 (stochastic gradient boosting)

### Algorithm
Custom implementation of gradient boosting with:
- Sequential tree fitting to residuals
- Squared loss function (L2 regression)
- Variance reduction for split selection
- Feature subsampling for each tree

---

## Model Performance

### Training Set
- **RMSE:** 0.011822
- **MAE:** 0.006024
- **R²:** 0.3987

### Test Set
- **RMSE:** 0.014595
- **MAE:** 0.009566
- **R²:** -0.5173

### Interpretation

1. **Training performance:** The model achieves reasonable performance on the training set with R² = 0.40, explaining about 40% of the variance in ln(γ).

2. **Test performance:** The negative R² on the test set indicates **significant overfitting**. The model performs worse than simply predicting the mean value for test samples.

3. **RMSE comparison:** Given the target range of ~0.14 units, the test RMSE of 0.0146 represents about 10% of the range, but the negative R² suggests poor generalization.

---

## Feature Importance Analysis

### Top 15 Most Important Features

| Rank | Feature | Importance |
|------|---------|------------|
| 1 | cation_1_molecular_weight | 0.2592 |
| 2 | B_MX_0_original | 0.2368 |
| 3 | r_M_angstrom (cation radius) | 0.0809 |
| 4 | max_molality_original | 0.0770 |
| 5 | C_MX_phi_original | 0.0718 |
| 6 | anion_1_n_atoms | 0.0705 |
| 7 | cation_1_radius_hydrated | 0.0286 |
| 8 | anion_1_diffusion_coefficient | 0.0282 |
| 9 | cation_1_diffusion_coefficient | 0.0262 |
| 10 | electrolyte_type_numeric | 0.0252 |
| 11 | anion_1_viscosity_jones_dole | 0.0238 |
| 12 | B_MX_1_original | 0.0233 |
| 13 | anion_1_delta_G_hydration | 0.0220 |
| 14 | anion_1_molecular_weight | 0.0204 |
| 15 | cation_1_delta_G_hydration | 0.0061 |

### Key Insights

1. **Pitzer parameters dominate:** B_MX_0 and C_MX_phi are among the top 5 most important features, which makes physical sense as these parameters directly capture ion-ion interactions.

2. **Cation properties matter more:** The most important feature is cation molecular weight, and several cation-specific properties rank highly (cation radius, hydrated radius, diffusion coefficient).

3. **Size and charge effects:** Cation radius (r_M_angstrom) is the 3rd most important feature, indicating that ion size plays a crucial role in determining activity coefficients.

4. **Anion contribution:** Anion properties (n_atoms, diffusion coefficient, viscosity coefficient) also contribute but generally rank lower than cation properties.

5. **Unused features:** Several features have zero importance (r_X_angstrom, anion_1_radius_hydrated, cation_1_n_atoms, cation_1_viscosity_jones_dole, B_MX_2_original), suggesting they may be redundant or uninformative for this target.

---

## Learning Curves Analysis

The training iterations show:
- **Training RMSE:** Decreases steadily from 0.0151 (iteration 1) to 0.0118 (iteration 100)
- **Test RMSE:** Increases from 0.0123 (iteration 1) to 0.0146 (iteration 100)

This pattern is a **classic sign of overfitting**:
- The model continues to improve on training data
- Performance on test data degrades after early iterations
- Early stopping around iteration 1-10 would have been optimal

---

## Worst Predictions (Test Set)

| Electrolyte | Observed | Predicted | Residual | Error Type |
|-------------|----------|-----------|----------|------------|
| ZnCl2 | 0.0439 | -0.0093 | 0.0532 | Large underprediction |
| LiClO4 | 0.0136 | -0.0108 | 0.0244 | Large underprediction |
| CsBr | 0.0195 | 0.0042 | 0.0153 | Underprediction |
| ZnI2 | -0.0269 | -0.0082 | -0.0187 | Overprediction |

### Error Patterns

1. **Zinc salts:** ZnCl2 and ZnI2 are poorly predicted, suggesting the model struggles with divalent transition metal cations.

2. **High ln(γ) values:** The largest errors occur for samples with higher observed ln(γ) values (ZnCl2, LiClO4, CsBr), indicating the model may be biased toward predicting values near zero.

3. **Lithium salts:** LiClO4 is poorly predicted despite LiBr being in the training set, suggesting sensitivity to anion type.

---

## Recommendations

### Model Improvements

1. **Reduce overfitting:**
   - Implement early stopping (stop at iteration 5-10)
   - Increase regularization (lower learning rate, increase min_leaf_size)
   - Reduce model complexity (fewer trees, shallower trees)
   - Add L2 regularization to leaf predictions

2. **Feature engineering:**
   - Create interaction terms (cation × anion properties)
   - Add ratios (e.g., r_M/r_X, molecular weight ratio)
   - Include charge-to-radius ratios (charge density)
   - Consider polynomial features for Pitzer parameters

3. **Cross-validation:**
   - Use k-fold CV to tune hyperparameters systematically
   - Test different learning rates (0.01, 0.05, 0.2)
   - Optimize number of trees based on CV performance

4. **Data augmentation:**
   - Collect more data, especially for underrepresented salts
   - Include other RH values as features or multi-task targets
   - Consider semi-supervised learning with unlabeled data

5. **Alternative approaches:**
   - Compare with simpler models (linear regression, regularized regression)
   - Try random forests (may generalize better)
   - Ensemble multiple models (bagging)
   - Consider physics-informed neural networks

### Physics Validation

1. **Check thermodynamic consistency:**
   - Ensure predictions respect Gibbs-Duhem relation
   - Validate against Debye-Hückel limiting law at low concentrations

2. **Investigate outliers:**
   - Why do zinc salts perform poorly? (coordination chemistry, hydrolysis)
   - Are there systematic errors by ion type or charge?

3. **Compare with Pitzer model:**
   - Benchmark against full Pitzer equation predictions
   - Identify where ML adds value vs. classical theory

---

## Conclusions

1. **Modest success:** The gradient boosting model captures some patterns in activity coefficients (training R² = 0.40) but struggles to generalize (test R² = -0.52).

2. **Overfitting is severe:** The model memorizes training data rather than learning generalizable patterns, likely due to:
   - Small dataset (108 samples, 22 features)
   - High model complexity (100 trees)
   - Lack of regularization

3. **Physical insights:**
   - Cation molecular weight and Pitzer parameters are most predictive
   - Ion size (radii) matters significantly
   - Some features are uninformative and could be removed

4. **Path forward:**
   - Simplify the model (fewer trees, early stopping, more regularization)
   - Engineer better features based on physical chemistry
   - Collect more data or use transfer learning from related systems
   - Consider hybrid physics-ML approaches

---

## Files Generated

- `gradient_boosting_results.txt` - Detailed numerical results
- `predictions.csv` - All predictions with residuals
- `gradient_boosting.m` - Model training script
- `RESULTS_SUMMARY.md` - This summary document

---

## Technical Notes

- Implementation: Custom gradient boosting in Octave/MATLAB
- Random seed: 42 (for reproducibility)
- Missing values: None in feature matrix after filtering
- Feature scaling: Standardized (mean=0, std=1)
- Graphics: Not generated (no graphics toolkit available in headless environment)

---

*Analysis completed: January 25, 2026*
