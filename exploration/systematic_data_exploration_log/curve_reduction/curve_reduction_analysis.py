#!/usr/bin/env python3
"""
Curve Features Analysis
This script analyzes the curve features from the curve reduction analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')

# Load data
print("Loading data...")
data = pd.read_csv('../../../figures/curve_reduction/curve_features.csv')

# Remove empty rows
data = data.dropna(subset=['Salt'])

# Extract features
salts = data['Salt'].values
slope = data['Initial_Slope_dln_aw_dm'].values
curvature = data['Curvature_d2ln_aw_dm2'].values
threshold = data['Threshold_Molality_aw_075'].values
integral = data['Integral_ln_aw'].values
max_molality = data['Max_Molality'].values

# Replace Inf with NaN for plotting
threshold = np.where(np.isinf(threshold), np.nan, threshold)

# Define efficiency and capacity
efficiency = np.abs(slope)  # Higher absolute slope = more efficient
capacity = np.abs(integral)  # Higher absolute integral = more capacity

n_salts = len(salts)

print(f"Analyzing {n_salts} salts...")

# Figure 1: Slope vs Curvature
print("Creating Figure 1: Slope vs Curvature...")
fig, ax = plt.subplots(figsize=(12, 9))

scatter = ax.scatter(slope, curvature, s=100, c=capacity, cmap='viridis', 
                     edgecolors='k', linewidth=0.5, alpha=0.8)

# Add text labels
for i, salt in enumerate(salts):
    ax.text(slope[i], curvature[i], f'  {salt}', fontsize=8, 
            ha='left', va='center')

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Capacity (|Integral ln $a_w$|)', fontsize=12)
ax.set_xlabel('Initial Slope (d ln $a_w$ / dm)', fontsize=12)
ax.set_ylabel('Curvature (d² ln $a_w$ / dm²)', fontsize=12)
ax.set_title('Salt Properties: Slope vs Curvature', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/slope_vs_curvature.png', dpi=300, bbox_inches='tight')
print("  Saved: slope_vs_curvature.png")
plt.close()

# Figure 2: Efficiency vs Capacity
print("Creating Figure 2: Efficiency vs Capacity...")
fig, ax = plt.subplots(figsize=(12, 9))

scatter = ax.scatter(capacity, efficiency, s=100, c=max_molality, cmap='viridis',
                     edgecolors='k', linewidth=0.5, alpha=0.8)

# Add text labels
for i, salt in enumerate(salts):
    ax.text(capacity[i], efficiency[i], f'  {salt}', fontsize=8,
            ha='left', va='center')

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Max Molality', fontsize=12)
ax.set_xlabel('Capacity (|Integral ln $a_w$|)', fontsize=12)
ax.set_ylabel('Efficiency (|Initial Slope|)', fontsize=12)
ax.set_title('Salt Performance: Efficiency vs Capacity', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

# Add quadrant lines
ax.axvline(np.median(capacity), color='k', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(np.median(efficiency), color='k', linestyle='--', linewidth=1, alpha=0.5)

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/efficiency_vs_capacity.png', dpi=300, bbox_inches='tight')
print("  Saved: efficiency_vs_capacity.png")
plt.close()

# Figure 3: Parallel Coordinates Plot
print("Creating Figure 3: Parallel Coordinates Plot...")

# Create feature matrix
features_raw = np.column_stack([np.abs(slope), curvature, max_molality, capacity])
feature_names = ['|Slope|', 'Curvature', 'Max Molality', 'Capacity']

# Normalize features to [0, 1]
features_norm = (features_raw - features_raw.min(axis=0)) / (features_raw.max(axis=0) - features_raw.min(axis=0))

fig, ax = plt.subplots(figsize=(14, 8))

# Plot each salt as a polyline
cmap = plt.cm.viridis(np.linspace(0, 1, n_salts))
for i in range(n_salts):
    ax.plot(range(4), features_norm[i, :], '-o', linewidth=1.5, 
            color=cmap[i], alpha=0.6, markersize=6)

ax.set_xticks(range(4))
ax.set_xticklabels(feature_names, fontsize=11)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 1)
ax.set_ylabel('Normalized Value', fontsize=12)
ax.set_title('Parallel Coordinates Plot of Salt Features', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Add legend for top salts by capacity
top_idx = np.argsort(capacity)[::-1][:10]
for idx in top_idx:
    ax.plot([], [], '-o', color=cmap[idx], linewidth=2, label=salts[idx])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/parallel_coordinates.png', dpi=300, bbox_inches='tight')
print("  Saved: parallel_coordinates.png")
plt.close()

# Clustering Analysis
print("\n=== Clustering Analysis ===")

# Determine optimal number of clusters
max_clusters = 6
silhouette_scores = []

for k in range(2, max_clusters + 1):
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
    cluster_labels = kmeans.fit_predict(features_norm)
    silhouette_avg = silhouette_score(features_norm, cluster_labels)
    silhouette_scores.append(silhouette_avg)
    print(f"  k={k}: silhouette = {silhouette_avg:.3f}")

# Find optimal k
optimal_k = np.argmax(silhouette_scores) + 2
best_silhouette = silhouette_scores[optimal_k - 2]

print(f"\nOptimal number of clusters: {optimal_k} (silhouette = {best_silhouette:.3f})")

# Perform final clustering
kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=20)
cluster_idx = kmeans.fit_predict(features_norm)
centroids = kmeans.cluster_centers_

# Display cluster information
for k in range(optimal_k):
    cluster_salts = salts[cluster_idx == k]
    print(f"\nCluster {k+1} ({len(cluster_salts)} salts):")
    for salt in cluster_salts:
        print(f"  - {salt}")

# Figure 4: Silhouette Plot
print("\nCreating Figure 4: Silhouette Analysis...")
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(range(2, max_clusters + 1), silhouette_scores, '-o', linewidth=2, markersize=8)
ax.plot(optimal_k, best_silhouette, 'ro', markersize=12, linewidth=2)
ax.set_xlabel('Number of Clusters', fontsize=12)
ax.set_ylabel('Mean Silhouette Score', fontsize=12)
ax.set_title('Optimal Number of Clusters', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/silhouette_analysis.png', dpi=300, bbox_inches='tight')
print("  Saved: silhouette_analysis.png")
plt.close()

# Figure 5: Clustered Slope vs Curvature
print("Creating Figure 5: Clustered Slope vs Curvature...")
fig, ax = plt.subplots(figsize=(12, 9))

colors = plt.cm.tab10(np.linspace(0, 1, optimal_k))

for k in range(optimal_k):
    mask = cluster_idx == k
    ax.scatter(slope[mask], curvature[mask], s=150, c=[colors[k]], 
              edgecolors='k', linewidth=1.5, label=f'Cluster {k+1}', alpha=0.8)

# Add labels
for i, salt in enumerate(salts):
    ax.text(slope[i], curvature[i], f'  {salt}', fontsize=8,
            ha='left', va='center', color=colors[cluster_idx[i]], fontweight='bold')

ax.set_xlabel('Initial Slope (d ln $a_w$ / dm)', fontsize=12)
ax.set_ylabel('Curvature (d² ln $a_w$ / dm²)', fontsize=12)
ax.set_title(f'Clustered Salts: Slope vs Curvature ({optimal_k} Clusters)', 
            fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10, loc='best')

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/clustered_slope_vs_curvature.png', dpi=300, bbox_inches='tight')
print("  Saved: clustered_slope_vs_curvature.png")
plt.close()

# Figure 6: Clustered Efficiency vs Capacity
print("Creating Figure 6: Clustered Efficiency vs Capacity...")
fig, ax = plt.subplots(figsize=(12, 9))

for k in range(optimal_k):
    mask = cluster_idx == k
    ax.scatter(capacity[mask], efficiency[mask], s=150, c=[colors[k]], 
              edgecolors='k', linewidth=1.5, label=f'Cluster {k+1}', alpha=0.8)

# Add labels
for i, salt in enumerate(salts):
    ax.text(capacity[i], efficiency[i], f'  {salt}', fontsize=8,
            ha='left', va='center', color=colors[cluster_idx[i]], fontweight='bold')

ax.set_xlabel('Capacity (|Integral ln $a_w$|)', fontsize=12)
ax.set_ylabel('Efficiency (|Initial Slope|)', fontsize=12)
ax.set_title(f'Clustered Salts: Efficiency vs Capacity ({optimal_k} Clusters)', 
            fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10, loc='best')

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/clustered_efficiency_vs_capacity.png', dpi=300, bbox_inches='tight')
print("  Saved: clustered_efficiency_vs_capacity.png")
plt.close()

# Figure 7: Parallel Coordinates by Cluster
print("Creating Figure 7: Parallel Coordinates by Cluster...")
nrows = int(np.ceil(optimal_k / 2))
ncols = 2
fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4*nrows))
axes = axes.flatten() if optimal_k > 1 else [axes]

for k in range(optimal_k):
    ax = axes[k]
    mask = cluster_idx == k
    
    # Plot all salts in gray in background
    for i in range(n_salts):
        ax.plot(range(4), features_norm[i, :], '-', linewidth=0.5,
               color='gray', alpha=0.2)
    
    # Plot cluster salts in color
    for i in np.where(mask)[0]:
        ax.plot(range(4), features_norm[i, :], '-o', linewidth=2,
               color=colors[k], markersize=6, alpha=0.8)
    
    # Plot cluster centroid
    ax.plot(range(4), centroids[k, :], 'k-', linewidth=3, label='Centroid')
    
    ax.set_xticks(range(4))
    ax.set_xticklabels(feature_names, fontsize=9)
    ax.set_xlim(-0.5, 3.5)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Normalized Value', fontsize=10)
    ax.set_title(f'Cluster {k+1} ({np.sum(mask)} salts)', fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

# Hide empty subplots
for k in range(optimal_k, len(axes)):
    axes[k].set_visible(False)

plt.suptitle('Parallel Coordinates by Cluster', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/parallel_coordinates_by_cluster.png', dpi=300, bbox_inches='tight')
print("  Saved: parallel_coordinates_by_cluster.png")
plt.close()

# Figure 8: Cluster Characteristics
print("Creating Figure 8: Cluster Centroids...")
fig, ax = plt.subplots(figsize=(12, 8))

x = np.arange(len(feature_names))
width = 0.8 / optimal_k

for k in range(optimal_k):
    offset = (k - optimal_k/2 + 0.5) * width
    ax.bar(x + offset, centroids[k, :], width, label=f'Cluster {k+1}',
          color=colors[k], edgecolor='k', linewidth=1)

ax.set_xticks(x)
ax.set_xticklabels(feature_names, fontsize=11)
ax.set_ylabel('Normalized Value', fontsize=12)
ax.set_title('Cluster Centroids', fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='best')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('../../../figures/curve_reduction/cluster_centroids.png', dpi=300, bbox_inches='tight')
print("  Saved: cluster_centroids.png")
plt.close()

# Summary Statistics
print("\n=== Cluster Summary Statistics ===")
for k in range(optimal_k):
    mask = cluster_idx == k
    print(f"\nCluster {k+1}:")
    print(f"  Mean Efficiency: {np.mean(efficiency[mask]):.4f}")
    print(f"  Mean Capacity: {np.mean(capacity[mask]):.4f}")
    print(f"  Mean Slope: {np.mean(slope[mask]):.4f}")
    print(f"  Mean Curvature: {np.mean(curvature[mask]):.4f}")
    print(f"  Mean Max Molality: {np.mean(max_molality[mask]):.4f}")

print("\n=== Analysis Complete ===")
print("All figures saved to ../../../figures/curve_reduction/")
