#!/usr/bin/env python3
"""
Create a proper pie chart for variance decomposition in mixed effects model.
This ensures the pie chart fills a full 360 degrees with all components.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

# Read variance values from the summary file
summary_file = '../../../figures/mixed_effects/mixed_effects_summary.txt'
output_file = '../../../figures/mixed_effects/variance_decomposition.png'

# Parse the summary file to get variance values
var_fixed = None
var_cation = None
var_anion = None
var_residual = None
var_total = None

with open(summary_file, 'r') as f:
    for line in f:
        if 'Fixed effect (concentration):' in line:
            parts = line.split()
            var_fixed = float(parts[3])
            pct_fixed = float(parts[4].strip('()%'))
        elif 'Random effect (cation):' in line:
            parts = line.split()
            var_cation = float(parts[3])
            pct_cation = float(parts[4].strip('()%'))
        elif 'Random effect (anion):' in line:
            parts = line.split()
            var_anion = float(parts[3])
            pct_anion = float(parts[4].strip('()%'))
        elif 'Residual:' in line and 'Variance Decomposition' not in line:
            parts = line.split()
            var_residual = float(parts[1])
            pct_residual = float(parts[2].strip('()%'))
        elif 'Total variance:' in line:
            parts = line.split()
            var_total = float(parts[2])

# Verify we got all values
if None in [var_fixed, var_cation, var_anion, var_residual, var_total]:
    print("Error: Could not parse all variance values from summary file")
    print(f"Fixed: {var_fixed}, Cation: {var_cation}, Anion: {var_anion}, Residual: {var_residual}, Total: {var_total}")
    exit(1)

# Create variance components array
variance_components = [var_fixed, var_cation, var_anion, var_residual]
labels = ['Fixed Effect (Concentration)', 'Cation Identity', 'Anion Identity', 'Residual']
colors = ['#3366CC', '#CC4D4D', '#4DB34D', '#B3B3B3']  # Blue, Red, Green, Gray

# Calculate percentages
percentages = [100 * v / sum(variance_components) for v in variance_components]

# Verify percentages sum to 100%
print(f"Variance components: {variance_components}")
print(f"Sum: {sum(variance_components):.6f}")
print(f"Percentages: {percentages}")
print(f"Sum of percentages: {sum(percentages):.1f}%")

# Create labels with percentages
labels_with_pct = [f'{label}\n({pct:.1f}%)' for label, pct in zip(labels, percentages)]

# Create the pie chart
fig, ax = plt.subplots(figsize=(10, 8))

# Create pie chart - this will automatically create a full 360 degree circle
wedges, texts, autotexts = ax.pie(
    variance_components,
    labels=labels_with_pct,
    colors=colors,
    autopct='',  # We're using labels_with_pct instead
    startangle=90,  # Start at top
    textprops={'fontsize': 11, 'fontweight': 'bold'}
)

# Make the pie chart fill the figure
ax.set_aspect('equal')

# Add title
plt.title('Variance Decomposition in ln(a_w)', fontsize=14, fontweight='bold', pad=20)

# Add total variance annotation
total_text = f'Total Variance: {var_total:.4f}'
fig.text(0.02, 0.02, total_text, fontsize=10, 
         bbox=dict(boxstyle='round', facecolor='white', edgecolor='black'))

# Save the figure
plt.tight_layout()
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPie chart saved to: {output_file}")
print(f"Chart shows full 360 degrees with {sum(percentages):.1f}% total")

plt.close()
