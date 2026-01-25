#!/usr/bin/env python3
"""
Create visualizations for gradient boosting results
This script generates plots that couldn't be created in Octave due to graphics toolkit issues.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
OUTPUT_DIR = Path('../../../figures/gradient_boosting/')
PREDICTIONS_FILE = OUTPUT_DIR / 'predictions.csv'

# Set style
plt.style.use('seaborn-v0_8-darkgrid')  # Use built-in style
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 10
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

def load_data():
    """Load predictions and parse results"""
    df = pd.read_csv(PREDICTIONS_FILE)
    
    # Separate train and test
    train_df = df[df['Dataset'] == 'Train'].copy()
    test_df = df[df['Dataset'] == 'Test'].copy()
    
    return train_df, test_df

def calculate_metrics(df):
    """Calculate performance metrics"""
    rmse = np.sqrt(np.mean(df['Residual']**2))
    mae = np.mean(np.abs(df['Residual']))
    ss_res = np.sum(df['Residual']**2)
    ss_tot = np.sum((df['Observed'] - df['Observed'].mean())**2)
    r2 = 1 - ss_res / ss_tot
    
    return {'RMSE': rmse, 'MAE': mae, 'R²': r2}

def plot_predicted_vs_observed(train_df, test_df):
    """Create predicted vs observed plots"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Training set
    train_metrics = calculate_metrics(train_df)
    ax1.scatter(train_df['Observed'], train_df['Predicted'], 
                alpha=0.6, s=60, c='steelblue', edgecolors='black', linewidth=0.5)
    
    # Perfect prediction line
    min_val = min(train_df['Observed'].min(), train_df['Predicted'].min())
    max_val = max(train_df['Observed'].max(), train_df['Predicted'].max())
    ax1.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='Perfect prediction')
    
    ax1.set_xlabel('Observed ln(γ) at 90% RH', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Predicted ln(γ) at 90% RH', fontsize=12, fontweight='bold')
    ax1.set_title(f'Training Set\nR² = {train_metrics["R²"]:.4f}, RMSE = {train_metrics["RMSE"]:.4f}', 
                  fontsize=13, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='box')
    
    # Test set
    test_metrics = calculate_metrics(test_df)
    ax2.scatter(test_df['Observed'], test_df['Predicted'], 
                alpha=0.6, s=60, c='coral', edgecolors='black', linewidth=0.5)
    
    min_val = min(test_df['Observed'].min(), test_df['Predicted'].min())
    max_val = max(test_df['Observed'].max(), test_df['Predicted'].max())
    ax2.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='Perfect prediction')
    
    ax2.set_xlabel('Observed ln(γ) at 90% RH', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Predicted ln(γ) at 90% RH', fontsize=12, fontweight='bold')
    ax2.set_title(f'Test Set\nR² = {test_metrics["R²"]:.4f}, RMSE = {test_metrics["RMSE"]:.4f}', 
                  fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'predicted_vs_observed.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: predicted_vs_observed.png")
    plt.close()

def plot_residuals(train_df, test_df):
    """Create residual plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Residuals vs predicted (train)
    axes[0, 0].scatter(train_df['Predicted'], train_df['Residual'], 
                       alpha=0.6, s=50, c='steelblue', edgecolors='black', linewidth=0.5)
    axes[0, 0].axhline(y=0, color='red', linestyle='--', linewidth=2)
    axes[0, 0].set_xlabel('Predicted ln(γ)', fontsize=11, fontweight='bold')
    axes[0, 0].set_ylabel('Residuals', fontsize=11, fontweight='bold')
    axes[0, 0].set_title('Training Set: Residuals vs Predicted', fontsize=12, fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Residuals vs predicted (test)
    axes[0, 1].scatter(test_df['Predicted'], test_df['Residual'], 
                       alpha=0.6, s=50, c='coral', edgecolors='black', linewidth=0.5)
    axes[0, 1].axhline(y=0, color='red', linestyle='--', linewidth=2)
    axes[0, 1].set_xlabel('Predicted ln(γ)', fontsize=11, fontweight='bold')
    axes[0, 1].set_ylabel('Residuals', fontsize=11, fontweight='bold')
    axes[0, 1].set_title('Test Set: Residuals vs Predicted', fontsize=12, fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Residuals histogram (train)
    axes[1, 0].hist(train_df['Residual'], bins=25, color='steelblue', 
                    edgecolor='black', alpha=0.7)
    axes[1, 0].axvline(x=0, color='red', linestyle='--', linewidth=2)
    axes[1, 0].set_xlabel('Residuals', fontsize=11, fontweight='bold')
    axes[1, 0].set_ylabel('Frequency', fontsize=11, fontweight='bold')
    axes[1, 0].set_title(f'Training Set: Residual Distribution\nMean = {train_df["Residual"].mean():.4f}', 
                         fontsize=12, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3, axis='y')
    
    # Residuals histogram (test)
    axes[1, 1].hist(test_df['Residual'], bins=15, color='coral', 
                    edgecolor='black', alpha=0.7)
    axes[1, 1].axvline(x=0, color='red', linestyle='--', linewidth=2)
    axes[1, 1].set_xlabel('Residuals', fontsize=11, fontweight='bold')
    axes[1, 1].set_ylabel('Frequency', fontsize=11, fontweight='bold')
    axes[1, 1].set_title(f'Test Set: Residual Distribution\nMean = {test_df["Residual"].mean():.4f}', 
                         fontsize=12, fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'residuals_analysis.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: residuals_analysis.png")
    plt.close()

def plot_error_analysis(test_df):
    """Analyze prediction errors by electrolyte"""
    # Sort by absolute error
    test_df['AbsError'] = np.abs(test_df['Residual'])
    test_sorted = test_df.sort_values('AbsError', ascending=False)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Top 10 worst predictions
    top10 = test_sorted.head(10)
    colors = ['coral' if x > 0 else 'steelblue' for x in top10['Residual']]
    
    ax1.barh(range(len(top10)), top10['Residual'], color=colors, edgecolor='black', linewidth=1)
    ax1.set_yticks(range(len(top10)))
    ax1.set_yticklabels(top10['Electrolyte'])
    ax1.axvline(x=0, color='black', linewidth=2)
    ax1.set_xlabel('Residual (Observed - Predicted)', fontsize=12, fontweight='bold')
    ax1.set_title('Top 10 Worst Predictions (Test Set)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.invert_yaxis()
    
    # Error vs observed value
    ax2.scatter(test_df['Observed'], test_df['AbsError'], 
                alpha=0.7, s=100, c=test_df['AbsError'], cmap='Reds', 
                edgecolors='black', linewidth=1)
    
    # Add labels for worst predictions
    for idx, row in top10.head(5).iterrows():
        ax2.annotate(row['Electrolyte'], 
                     xy=(row['Observed'], row['AbsError']),
                     xytext=(5, 5), textcoords='offset points',
                     fontsize=9, alpha=0.8)
    
    ax2.set_xlabel('Observed ln(γ) at 90% RH', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Absolute Error', fontsize=12, fontweight='bold')
    ax2.set_title('Absolute Error vs Observed Value', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'error_analysis.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: error_analysis.png")
    plt.close()

def plot_feature_importance():
    """Plot feature importance from results file"""
    # Parse feature importance from results file
    results_file = OUTPUT_DIR / 'gradient_boosting_results.txt'
    
    with open(results_file, 'r') as f:
        lines = f.readlines()
    
    # Find feature importance section
    features = []
    importances = []
    in_importance_section = False
    
    for line in lines:
        if 'Top 20 Feature Importances:' in line:
            in_importance_section = True
            continue
        if in_importance_section and line.strip():
            parts = line.strip().split()
            if len(parts) >= 3 and parts[0].endswith('.'):
                feature = ' '.join(parts[1:-1])
                importance = float(parts[-1])
                features.append(feature)
                importances.append(importance)
                if len(features) >= 15:  # Get top 15
                    break
    
    if not features:
        print("⚠ Could not parse feature importances")
        return
    
    # Create horizontal bar plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = np.arange(len(features))
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(features)))
    
    ax.barh(y_pos, importances, color=colors, edgecolor='black', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(features, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('Importance', fontsize=12, fontweight='bold')
    ax.set_title('Top 15 Feature Importances', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'feature_importance.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: feature_importance.png")
    plt.close()

def main():
    """Main execution"""
    print("\n" + "="*60)
    print("Creating Gradient Boosting Visualizations")
    print("="*60 + "\n")
    
    # Load data
    print("Loading predictions...")
    train_df, test_df = load_data()
    print(f"  - Training samples: {len(train_df)}")
    print(f"  - Test samples: {len(test_df)}\n")
    
    # Generate plots
    print("Generating visualizations...")
    plot_predicted_vs_observed(train_df, test_df)
    plot_residuals(train_df, test_df)
    plot_error_analysis(test_df)
    plot_feature_importance()
    
    # Summary
    print(f"\n{'='*60}")
    print("Summary Statistics")
    print("="*60)
    
    train_metrics = calculate_metrics(train_df)
    test_metrics = calculate_metrics(test_df)
    
    print("\nTraining Set:")
    for metric, value in train_metrics.items():
        print(f"  {metric}: {value:.6f}")
    
    print("\nTest Set:")
    for metric, value in test_metrics.items():
        print(f"  {metric}: {value:.6f}")
    
    print(f"\n{'='*60}")
    print(f"All visualizations saved to: {OUTPUT_DIR}")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()
