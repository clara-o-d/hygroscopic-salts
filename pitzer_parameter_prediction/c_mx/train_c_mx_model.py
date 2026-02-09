#!/usr/bin/env python3
"""
Train Machine Learning Models to Predict C_MX_phi Pitzer Coefficient
Using Baseline Numeric Data

This script demonstrates a complete ML workflow:
1. Load the baseline numeric dataset
2. Split data into train/test sets
3. Scale features
4. Train multiple models
5. Evaluate and compare performance
6. Save the best model
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import warnings
warnings.filterwarnings('ignore')


def load_data(filepath='../data/baseline_numeric_only.csv'):
    """Load the ML-ready baseline dataset and separate features and targets."""
    print("=" * 80)
    print("Loading Data")
    print("=" * 80)
    
    df = pd.read_csv(filepath)
    print(f"✓ Loaded {len(df)} samples with {len(df.columns)} columns")
    
    # Define target column (C_MX_phi Pitzer coefficient)
    target_col = 'C_MX_phi_original'
    
    # ID columns
    id_cols = ['electrolyte']
    
    # Columns to exclude (target, other Pitzer coefficients, reference info)
    exclude_cols = [
        target_col,
        'B_MX_0_original', 'B_MX_1_original', 'B_MX_2_original',  # Other Pitzer coefficients
    ]
    
    # Get all feature columns
    feature_cols = [col for col in df.columns 
                   if col not in id_cols + exclude_cols]
    
    # Remove rows with missing target values
    df = df.dropna(subset=[target_col])
    
    X = df[feature_cols]
    y = df[target_col]
    molecules = df['electrolyte']
    
    print(f"✓ Features: {X.shape[1]} columns")
    print(f"  Features: {', '.join(feature_cols[:10])}..." if len(feature_cols) > 10 else f"  Features: {', '.join(feature_cols)}")
    print(f"✓ Target: {target_col}")
    print(f"✓ Samples after removing missing targets: {len(X)}")
    print(f"✓ Missing values: {X.isna().sum().sum()} in features, {y.isna().sum()} in target")
    
    return X, y, molecules, feature_cols, target_col


def explore_data(X, y, feature_cols, target_col):
    """Display basic statistics about the dataset."""
    print("\n" + "=" * 80)
    print("Data Exploration")
    print("=" * 80)
    
    print("\nFeatures Statistics (first 10):")
    print(X.describe().T[['mean', 'std', 'min', 'max']].head(10))
    print("...")
    
    print("\nTarget Statistics:")
    print(f"  {target_col}:")
    print(f"    Mean: {y.mean():.6f}")
    print(f"    Std: {y.std():.6f}")
    print(f"    Min: {y.min():.6f}")
    print(f"    Max: {y.max():.6f}")
    print(f"    Non-zero count: {(y != 0).sum()} / {len(y)}")


def split_and_scale_data(X, y, train_size=0.75, val_size=0.18, test_size=0.07, random_state=42):
    """Split data into train/validation/test, impute missing values, and scale features."""
    print("\n" + "=" * 80)
    print("Data Preparation (Train/Validation/Test Split)")
    print("=" * 80)
    
    # First split: separate training from temp (val + test)
    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=(val_size + test_size), random_state=random_state
    )
    
    # Second split: separate validation from test
    val_prop = val_size / (val_size + test_size)
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=(1 - val_prop), random_state=random_state
    )
    
    print(f"✓ Train set: {len(X_train)} samples ({train_size*100:.0f}%)")
    print(f"✓ Validation set: {len(X_val)} samples ({val_size*100:.0f}%)")
    print(f"✓ Test set: {len(X_test)} samples ({test_size*100:.0f}%)")
    print(f"✓ Total: {len(X)} samples")
    
    # Impute missing values using training data statistics
    imputer = SimpleImputer(strategy='mean')
    X_train_imputed = imputer.fit_transform(X_train)
    X_val_imputed = imputer.transform(X_val)
    X_test_imputed = imputer.transform(X_test)
    
    print(f"✓ Missing values imputed (mean strategy) using training data only")
    print(f"  Missing values in train: {np.isnan(X_train_imputed).sum()}")
    print(f"  Missing values in val: {np.isnan(X_val_imputed).sum()}")
    print(f"  Missing values in test: {np.isnan(X_test_imputed).sum()}")
    
    # Scale features using only training data statistics
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_imputed)
    X_val_scaled = scaler.transform(X_val_imputed)
    X_test_scaled = scaler.transform(X_test_imputed)
    
    print(f"✓ Features scaled (mean=0, std=1) using training data only")
    
    return (X_train, X_val, X_test,
            X_train_scaled, X_val_scaled, X_test_scaled,
            y_train, y_val, y_test, scaler, imputer)


def train_models(X_train_scaled, y_train):
    """Train multiple models and return them."""
    print("\n" + "=" * 80)
    print("Training Models")
    print("=" * 80)
    
    models = {}
    
    # 1. Random Forest
    print("\n[1/5] Training Random Forest...")
    print(f"  Using {X_train_scaled.shape[1]} features")
    rf = RandomForestRegressor(
        n_estimators=150,
        max_depth=12,
        min_samples_split=5,
        min_samples_leaf=2,
        random_state=42,
        n_jobs=1
    )
    rf.fit(X_train_scaled, y_train)
    models['Random Forest'] = rf
    print("  ✓ Complete")
    
    # 2. Gradient Boosting
    print("\n[2/5] Training Gradient Boosting...")
    print(f"  Using {X_train_scaled.shape[1]} features")
    gb = GradientBoostingRegressor(
        n_estimators=150,
        max_depth=5,
        learning_rate=0.1,
        random_state=42
    )
    gb.fit(X_train_scaled, y_train)
    models['Gradient Boosting'] = gb
    print("  ✓ Complete")
    
    # 3. Ridge Regression
    print("\n[3/5] Training Ridge Regression...")
    print(f"  Using {X_train_scaled.shape[1]} features")
    ridge = Ridge(alpha=1.0, random_state=42)
    ridge.fit(X_train_scaled, y_train)
    models['Ridge'] = ridge
    print("  ✓ Complete")
    
    # 4. Lasso Regression
    print("\n[4/5] Training Lasso Regression...")
    print(f"  Using {X_train_scaled.shape[1]} features")
    lasso = Lasso(alpha=0.01, max_iter=10000, random_state=42)
    lasso.fit(X_train_scaled, y_train)
    models['Lasso'] = lasso
    print("  ✓ Complete")
    
    # 5. Elastic Net
    print("\n[5/5] Training Elastic Net...")
    print(f"  Using {X_train_scaled.shape[1]} features")
    elastic = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=10000, random_state=42)
    elastic.fit(X_train_scaled, y_train)
    models['Elastic Net'] = elastic
    print("  ✓ Complete")
    
    return models


def evaluate_models(models, X_train_scaled, X_val_scaled, y_train, y_val, dataset_name="Validation"):
    """Evaluate all models on validation set (used for model selection)."""
    print("\n" + "=" * 80)
    print(f"Model Evaluation ({dataset_name} Set - For Model Selection)")
    print("=" * 80)
    
    results = []
    
    for model_name, model in models.items():
        print(f"\n{model_name}:")
        print("-" * 40)
        
        # Predictions
        y_train_pred = model.predict(X_train_scaled)
        y_val_pred = model.predict(X_val_scaled)
        
        # Metrics
        train_r2 = r2_score(y_train, y_train_pred)
        train_rmse = np.sqrt(mean_squared_error(y_train, y_train_pred))
        
        val_r2 = r2_score(y_val, y_val_pred)
        val_rmse = np.sqrt(mean_squared_error(y_val, y_val_pred))
        val_mae = mean_absolute_error(y_val, y_val_pred)
        
        print(f"  Train R² = {train_r2:.4f}, RMSE = {train_rmse:.6f}")
        print(f"  Val   R² = {val_r2:.4f}, RMSE = {val_rmse:.6f}, MAE = {val_mae:.6f}")
        
        results.append({
            'model': model_name,
            'val_r2': val_r2,
            'val_rmse': val_rmse,
            'val_mae': val_mae
        })
    
    return pd.DataFrame(results)


def evaluate_final_model(model, X_test_scaled, y_test, model_name="Final Model"):
    """Evaluate the final selected model on test set (ONLY USED ONCE)."""
    print("\n" + "=" * 80)
    print("FINAL MODEL EVALUATION ON TEST SET")
    print("=" * 80)
    
    y_test_pred = model.predict(X_test_scaled)
    
    test_r2 = r2_score(y_test, y_test_pred)
    test_rmse = np.sqrt(mean_squared_error(y_test, y_test_pred))
    test_mae = mean_absolute_error(y_test, y_test_pred)
    
    print(f"\n{model_name} - Test Set Performance:")
    print("-" * 60)
    print(f"  Test R² = {test_r2:.4f}, RMSE = {test_rmse:.6f}, MAE = {test_mae:.6f}")
    print("\n" + "=" * 80)
    
    results = {
        'test_r2': test_r2,
        'test_rmse': test_rmse,
        'test_mae': test_mae
    }
    
    return results, y_test_pred


def cross_validate_best_model(best_model, X_train_scaled, y_train, cv=5):
    """Perform cross-validation on the best model."""
    print("\n" + "=" * 80)
    print("Cross-Validation (Best Model)")
    print("=" * 80)
    
    kfold = KFold(n_splits=cv, shuffle=True, random_state=42)
    cv_scores = cross_val_score(best_model, X_train_scaled, y_train, cv=kfold, scoring='r2', n_jobs=1)
    
    print(f"\n{cv}-Fold Cross-Validation R² Scores:")
    for i, score in enumerate(cv_scores, 1):
        print(f"  Fold {i}: {score:.4f}")
    
    print(f"\nMean R²: {cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")
    
    return cv_scores


def plot_predictions(models, X_val_scaled, y_val, output_file='c_mx_predictions.png', dataset_name="Validation"):
    """Plot predicted vs actual values for each model on validation set."""
    print("\n" + "=" * 80)
    print(f"Generating Prediction Plots ({dataset_name} Set)")
    print("=" * 80)
    
    n_models = len(models)
    
    fig, axes = plt.subplots(1, n_models, figsize=(7*n_models, 5))
    if n_models == 1:
        axes = [axes]
    
    fig.suptitle(f'C_MX_phi Prediction: Predicted vs Actual ({dataset_name} Set)', 
                fontsize=16, fontweight='bold', y=0.995)
    
    for i, (model_name, model) in enumerate(models.items()):
        ax = axes[i]
        y_pred = model.predict(X_val_scaled)
        
        # Scatter plot
        ax.scatter(y_val, y_pred, alpha=0.6, s=80, edgecolors='k', linewidth=0.5)
        
        # Perfect prediction line
        min_val = min(y_val.min(), y_pred.min())
        max_val = max(y_val.max(), y_pred.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2, label='Perfect')
        
        # Metrics
        r2 = r2_score(y_val, y_pred)
        rmse = np.sqrt(mean_squared_error(y_val, y_pred))
        
        ax.set_xlabel('Actual', fontsize=11)
        ax.set_ylabel('Predicted', fontsize=11)
        title = f'{model_name}\nR² = {r2:.4f}, RMSE = {rmse:.4f}'
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved to '{output_file}'")
    plt.close()
    
    return fig


def plot_best_model_test_predictions(best_model, best_model_name, X_test_scaled, y_test,
                                     output_file='best_model_test_predictions.png'):
    """Plot predicted vs actual values for the best model on test set."""
    print("\n" + "=" * 80)
    print("Generating Best Model Test Set Prediction Plot")
    print("=" * 80)
    
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    fig.suptitle('Test Set Predictions: Best Model', 
                fontsize=16, fontweight='bold', y=0.995)
    
    y_pred = best_model.predict(X_test_scaled)
    
    # Scatter plot
    ax.scatter(y_test, y_pred, alpha=0.7, s=100, edgecolors='k', linewidth=0.5, color='steelblue')
    
    # Perfect prediction line
    min_val = min(y_test.min(), y_pred.min())
    max_val = max(y_test.max(), y_pred.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2.5, label='Perfect Prediction')
    
    # Metrics
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    
    ax.set_xlabel('Actual C_MX_phi', fontsize=12)
    ax.set_ylabel('Predicted C_MX_phi', fontsize=12)
    title = f'{best_model_name}\nR² = {r2:.4f}, RMSE = {rmse:.4f}, MAE = {mae:.4f}'
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved to '{output_file}'")
    plt.close()
    
    return fig


def analyze_feature_importance(model, feature_cols, model_name, top_n=20):
    """Analyze and display feature importance."""
    print("\n" + "=" * 80)
    print(f"Feature Importance Analysis ({model_name})")
    print("=" * 80)
    
    importances = None
    
    # Extract feature importance based on model type
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_
        print(f"\nUsing built-in feature importances")
        
    elif hasattr(model, 'coef_'):
        print(f"\nUsing absolute coefficient values")
        importances = np.abs(model.coef_)
    
    if importances is not None:
        feature_importance = pd.DataFrame({
            'feature': feature_cols,
            'importance': importances
        }).sort_values('importance', ascending=False)
        
        print(f"\nTOP {top_n} MOST IMPORTANT FEATURES:")
        print("-" * 80)
        for rank, (idx, row) in enumerate(feature_importance.head(top_n).iterrows(), 1):
            print(f"  {rank:2d}. {row['feature']:45s} : {row['importance']:.6f}")
        
        return feature_importance
    else:
        print(f"\n⚠️ Feature importance not available for this model type")
        return None


def save_model(model, scaler, imputer, filepath='best_c_mx_model.pkl'):
    """Save the trained model, scaler, and imputer."""
    import pickle
    
    print("\n" + "=" * 80)
    print("Saving Model")
    print("=" * 80)
    
    model_data = {
        'model': model,
        'scaler': scaler,
        'imputer': imputer
    }
    
    with open(filepath, 'wb') as f:
        pickle.dump(model_data, f)
    
    print(f"✓ Model, scaler, and imputer saved to '{filepath}'")


def main():
    """Main execution function."""
    print("\n" + "=" * 80)
    print("C_MX_phi PITZER COEFFICIENT PREDICTION")
    print("Baseline Numeric Data")
    print("=" * 80)
    
    # Set up paths
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, '../data/baseline_numeric_only.csv')
    
    # 1. Load data
    X, y, molecules, feature_cols, target_col = load_data(data_path)
    
    # 2. Explore data
    explore_data(X, y, feature_cols, target_col)
    
    # 3. Split into train/validation/test, impute, and scale
    (X_train, X_val, X_test,
     X_train_scaled, X_val_scaled, X_test_scaled,
     y_train, y_val, y_test, scaler, imputer) = \
        split_and_scale_data(X, y, train_size=0.75, val_size=0.10, test_size=0.15, random_state=42)
    
    # 4. Train models
    models = train_models(X_train_scaled, y_train)
    
    # 5. Evaluate models on VALIDATION set (used for model selection)
    results_df = evaluate_models(models, X_train_scaled, X_val_scaled, y_train, y_val, dataset_name="Validation")
    
    # 6. Compare models based on VALIDATION performance
    print("\n" + "=" * 80)
    print("Model Comparison (Based on Validation Set)")
    print("=" * 80)
    print("\nValidation R² by Model:")
    comparison = results_df[['model', 'val_r2', 'val_rmse']].sort_values('val_r2', ascending=False)
    print(comparison.to_string(index=False))
    
    best_model_name = comparison.iloc[0]['model']
    best_model = models[best_model_name]
    
    print(f"\n🏆 Best Model (Selected on Validation): {best_model_name}")
    print(f"   Validation R²: {comparison.iloc[0]['val_r2']:.4f}")
    
    # 7. Cross-validation on training set only
    cv_scores = cross_validate_best_model(best_model, X_train_scaled, y_train)
    
    # 8. FINAL EVALUATION: Evaluate selected model ONCE on test set
    final_test_results, y_test_pred = evaluate_final_model(
        best_model, X_test_scaled, y_test, best_model_name
    )
    
    # 9. Plot predictions on validation set (all models)
    plot_output = os.path.join(script_dir, 'c_mx_predictions.png')
    plot_predictions(models, X_val_scaled, y_val, plot_output, dataset_name="Validation")
    
    # 9.5. Plot best model predictions on test set
    best_model_plot_output = os.path.join(script_dir, 'best_model_test_predictions.png')
    plot_best_model_test_predictions(
        best_model, best_model_name, X_test_scaled, y_test, best_model_plot_output
    )
    
    # 10. Feature importance
    feature_importance = analyze_feature_importance(
        best_model, feature_cols, best_model_name
    )
    
    # 11. Save best model
    model_path = os.path.join(script_dir, 'best_c_mx_model.pkl')
    save_model(best_model, scaler, imputer, model_path)
    
    # Final summary
    print("\n" + "=" * 80)
    print("TRAINING COMPLETE!")
    print("=" * 80)
    print(f"\nDataset: Baseline Numeric Data")
    print(f"  Samples: {len(X)} ({len(X_train)} train, {len(X_val)} validation, {len(X_test)} test)")
    print(f"  Features: {len(feature_cols)}")
    print(f"  Target: {target_col}")
    
    print(f"\nBest Model: {best_model_name}")
    print(f"  Selected based on Validation R²: {comparison.iloc[0]['val_r2']:.4f}")
    print(f"  Final Test R² (unbiased): {final_test_results['test_r2']:.4f}")
    print(f"  Final Test RMSE: {final_test_results['test_rmse']:.6f}")
    
    print("\nGenerated Files:")
    print("  - best_c_mx_model.pkl (trained model + scaler)")
    print("  - c_mx_predictions.png (all models on validation set)")
    print("  - best_model_test_predictions.png (best model on test set)")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
