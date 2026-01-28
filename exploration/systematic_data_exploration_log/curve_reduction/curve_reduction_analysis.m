%% Curve Features Analysis
% This script analyzes the curve features from the curve reduction analysis

clear; clc; close all;

%% Load Data
fprintf('Loading data...\n');
data = readtable('../../../figures/curve_reduction/curve_features.csv');

% Remove any empty rows
data = data(~cellfun(@isempty, data.Salt), :);

% Extract features
salts = data.Salt;
slope = data.Initial_Slope_dln_aw_dm;
curvature = data.Curvature_d2ln_aw_dm2;
threshold = data.Threshold_Molality_aw_075;
integral = data.Integral_ln_aw;
max_molality = data.Max_Molality;

% Replace Inf in threshold with NaN for plotting
threshold(isinf(threshold)) = NaN;

% Define efficiency and capacity
% Efficiency: steepness of response (higher slope magnitude = more efficient)
% Capacity: total water activity reduction (negative integral = more capacity)
efficiency = abs(slope);  % Higher absolute slope = more efficient
capacity = abs(integral);  % Higher absolute integral = more capacity

n_salts = length(salts);

fprintf('Analyzing %d salts...\n', n_salts);

%% Figure 1: Slope vs Curvature
fprintf('Creating Figure 1: Slope vs Curvature...\n');
figure('Position', [100, 100, 1200, 900]);

% Create scatter plot with capacity as color
scatter(slope, curvature, 100, capacity, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
hold on;

% Add text labels for each point
for i = 1:n_salts
    text(slope(i), curvature(i), ['  ' salts{i}], ...
        'FontSize', 8, 'HorizontalAlignment', 'left');
end

% Use viridis-like colormap (parula is similar in MATLAB)
colormap(parula);
cb = colorbar;
cb.Label.String = 'Capacity (|Integral ln a_w|)';
cb.Label.FontSize = 12;

xlabel('Initial Slope (d ln a_w / dm)', 'FontSize', 12);
ylabel('Curvature (d^2 ln a_w / dm^2)', 'FontSize', 12);
title('Salt Properties: Slope vs Curvature', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/slope_vs_curvature.png');
saveas(gcf, '../../../figures/curve_reduction/slope_vs_curvature.fig');
fprintf('  Saved: slope_vs_curvature.png\n');

%% Figure 2: Efficiency vs Capacity
fprintf('Creating Figure 2: Efficiency vs Capacity...\n');
figure('Position', [120, 120, 1200, 900]);

scatter(capacity, efficiency, 100, max_molality, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
hold on;

% Add text labels
for i = 1:n_salts
    text(capacity(i), efficiency(i), ['  ' salts{i}], ...
        'FontSize', 8, 'HorizontalAlignment', 'left');
end

colormap(parula);
cb = colorbar;
cb.Label.String = 'Max Molality';
cb.Label.FontSize = 12;

xlabel('Capacity (|Integral ln a_w|)', 'FontSize', 12);
ylabel('Efficiency (|Initial Slope|)', 'FontSize', 12);
title('Salt Performance: Efficiency vs Capacity', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add quadrant lines
xlims = xlim;
ylims = ylim;
plot([median(capacity), median(capacity)], ylims, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(xlims, [median(efficiency), median(efficiency)], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Save figure
saveas(gcf, '../../../figures/curve_reduction/efficiency_vs_capacity.png');
saveas(gcf, '../../../figures/curve_reduction/efficiency_vs_capacity.fig');
fprintf('  Saved: efficiency_vs_capacity.png\n');

%% Figure 3: Parallel Coordinates Plot
fprintf('Creating Figure 3: Parallel Coordinates Plot...\n');

% Create feature matrix (normalized)
features_raw = [abs(slope), curvature, max_molality, capacity];
feature_names = {'|Slope|', 'Curvature', 'Max Molality', 'Capacity'};

% Normalize features to [0, 1]
features_norm = zeros(size(features_raw));
for i = 1:size(features_raw, 2)
    feat = features_raw(:, i);
    features_norm(:, i) = (feat - min(feat)) / (max(feat) - min(feat));
end

figure('Position', [140, 140, 1600, 900]);

% Create distinguishable colors - combine multiple colormaps for variety
cmap_distinct = distinguishable_colors(n_salts);

% Plot each salt as a polyline with legend
hold on;
legend_handles = zeros(n_salts, 1);
for i = 1:n_salts
    legend_handles(i) = plot(1:4, features_norm(i, :), '-o', 'LineWidth', 1.8, ...
        'Color', cmap_distinct(i, :), 'MarkerSize', 7, ...
        'MarkerFaceColor', cmap_distinct(i, :), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    
    % Add text label on the line - position it at different points along the line
    % to distribute labels and reduce overlap
    % Use modulo to distribute labels across different segments
    segment = mod(i-1, 3) + 1;  % Cycle through 3 segments
    x_positions = [1.5, 2.5, 3.5];  % Between axes 1-2, 2-3, 3-4
    x_label = x_positions(segment);
    
    % Interpolate to find y position at the chosen x
    y_label = interp1(1:4, features_norm(i, :), x_label, 'linear');
    
    % Add small perpendicular offset to avoid overlap with line
    % Offset direction alternates
    offset_y = 0.015 * (-1)^i;  % Alternate up/down
    
    text(x_label, y_label + offset_y, salts{i}, ...
        'FontSize', 6.5, 'HorizontalAlignment', 'center', ...
        'Color', cmap_distinct(i, :), 'FontWeight', 'bold', ...
        'BackgroundColor', [1, 1, 1, 0.8], 'EdgeColor', cmap_distinct(i, :), ...
        'LineWidth', 0.5, 'Margin', 0.5);
end

% Customize axes
set(gca, 'XTick', 1:4, 'XTickLabel', feature_names);
xlim([0.5, 4.5]);
ylim([0, 1]);
ylabel('Normalized Value', 'FontSize', 12);
title('Parallel Coordinates Plot of Salt Features', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add full legend with all salts
legend(legend_handles, salts, 'Location', 'eastoutside', 'FontSize', 8, 'NumColumns', 1);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/parallel_coordinates.png');
saveas(gcf, '../../../figures/curve_reduction/parallel_coordinates.fig');
fprintf('  Saved: parallel_coordinates.png\n');

%% Helper function for distinguishable colors
function colors = distinguishable_colors(n)
    % Generate n maximally distinguishable colors
    % Based on generating colors in HSV space with varied hue
    if n <= 20
        % Use predefined distinct colors for small n
        base_colors = [
            0.00, 0.45, 0.74;  % blue
            0.85, 0.33, 0.10;  % red-orange
            0.93, 0.69, 0.13;  % yellow
            0.49, 0.18, 0.56;  % purple
            0.47, 0.67, 0.19;  % green
            0.30, 0.75, 0.93;  % cyan
            0.64, 0.08, 0.18;  % dark red
            0.74, 0.74, 0.13;  % yellow-green
            0.00, 0.62, 0.45;  % teal
            0.90, 0.17, 0.00;  % red
            0.50, 0.50, 0.50;  % gray
            0.00, 0.00, 0.00;  % black
            0.80, 0.40, 0.80;  % pink-purple
            0.20, 0.60, 0.80;  % light blue
            0.60, 0.40, 0.20;  % brown
            0.90, 0.60, 0.00;  % orange
            0.40, 0.20, 0.60;  % deep purple
            0.00, 0.50, 0.00;  % dark green
            0.80, 0.00, 0.40;  % magenta
            0.20, 0.20, 0.80;  % dark blue
        ];
        if n <= size(base_colors, 1)
            colors = base_colors(1:n, :);
        else
            colors = [base_colors; hsv(n - size(base_colors, 1))];
        end
    else
        % For larger n, use HSV with varied saturation and value
        colors = zeros(n, 3);
        hues = linspace(0, 1, n+1);
        hues = hues(1:n);
        
        for i = 1:n
            % Vary saturation and value to make colors more distinct
            sat = 0.6 + 0.4 * mod(i, 3) / 2;
            val = 0.7 + 0.3 * mod(i, 2);
            colors(i, :) = hsv2rgb([hues(i), sat, val]);
        end
    end
end

%% Clustering Analysis
fprintf('\n=== Clustering Analysis ===\n');

% Use normalized features for clustering
X = features_norm;

% Determine optimal number of clusters using silhouette analysis
max_clusters = 6;
silhouette_scores = zeros(max_clusters - 1, 1);

for k = 2:max_clusters
    idx_temp = kmeans(X, k, 'Replicates', 10, 'Distance', 'sqeuclidean');
    silhouette_scores(k-1) = mean(silhouette(X, idx_temp));
    fprintf('  k=%d: silhouette = %.3f\n', k, silhouette_scores(k-1));
end

% Find optimal k
[best_silhouette, best_k_idx] = max(silhouette_scores);
optimal_k = best_k_idx + 1;

fprintf('\nOptimal number of clusters: %d (silhouette = %.3f)\n', optimal_k, best_silhouette);

% Perform clustering with optimal k
[cluster_idx, centroids] = kmeans(X, optimal_k, 'Replicates', 20, 'Distance', 'sqeuclidean');

% Display cluster information
for k = 1:optimal_k
    cluster_salts = salts(cluster_idx == k);
    fprintf('\nCluster %d (%d salts):\n', k, length(cluster_salts));
    for i = 1:length(cluster_salts)
        fprintf('  - %s\n', cluster_salts{i});
    end
end

%% Figure 4: Silhouette Plot
fprintf('\nCreating Figure 4: Silhouette Analysis...\n');
figure('Position', [160, 160, 1000, 700]);

plot(2:max_clusters, silhouette_scores, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2 0.4 0.8]);
hold on;
plot(optimal_k, best_silhouette, 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'r');

xlabel('Number of Clusters', 'FontSize', 12);
ylabel('Mean Silhouette Score', 'FontSize', 12);
title('Optimal Number of Clusters', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/silhouette_analysis.png');
saveas(gcf, '../../../figures/curve_reduction/silhouette_analysis.fig');
fprintf('  Saved: silhouette_analysis.png\n');

%% Figure 5: Clustered Slope vs Curvature
fprintf('Creating Figure 5: Clustered Slope vs Curvature...\n');
figure('Position', [180, 180, 1200, 900]);

cmap_clusters = lines(optimal_k);

for k = 1:optimal_k
    idx_k = cluster_idx == k;
    scatter(slope(idx_k), curvature(idx_k), 150, cmap_clusters(k, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end

% Add labels
for i = 1:n_salts
    text(slope(i), curvature(i), ['  ' salts{i}], ...
        'FontSize', 8, 'HorizontalAlignment', 'left', ...
        'Color', cmap_clusters(cluster_idx(i), :), 'FontWeight', 'bold');
end

xlabel('Initial Slope (d ln a_w / dm)', 'FontSize', 12);
ylabel('Curvature (d^2 ln a_w / dm^2)', 'FontSize', 12);
title(sprintf('Clustered Salts: Slope vs Curvature (%d Clusters)', optimal_k), ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add legend
legend_entries = cell(optimal_k, 1);
for k = 1:optimal_k
    legend_entries{k} = sprintf('Cluster %d', k);
end
legend(legend_entries, 'Location', 'best', 'FontSize', 10);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/clustered_slope_vs_curvature.png');
saveas(gcf, '../../../figures/curve_reduction/clustered_slope_vs_curvature.fig');
fprintf('  Saved: clustered_slope_vs_curvature.png\n');

%% Figure 6: Clustered Efficiency vs Capacity
fprintf('Creating Figure 6: Clustered Efficiency vs Capacity...\n');
figure('Position', [200, 200, 1200, 900]);

for k = 1:optimal_k
    idx_k = cluster_idx == k;
    scatter(capacity(idx_k), efficiency(idx_k), 150, cmap_clusters(k, :), ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
end

% Add labels
for i = 1:n_salts
    text(capacity(i), efficiency(i), ['  ' salts{i}], ...
        'FontSize', 8, 'HorizontalAlignment', 'left', ...
        'Color', cmap_clusters(cluster_idx(i), :), 'FontWeight', 'bold');
end

xlabel('Capacity (|Integral ln a_w|)', 'FontSize', 12);
ylabel('Efficiency (|Initial Slope|)', 'FontSize', 12);
title(sprintf('Clustered Salts: Efficiency vs Capacity (%d Clusters)', optimal_k), ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
legend(legend_entries, 'Location', 'best', 'FontSize', 10);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/clustered_efficiency_vs_capacity.png');
saveas(gcf, '../../../figures/curve_reduction/clustered_efficiency_vs_capacity.fig');
fprintf('  Saved: clustered_efficiency_vs_capacity.png\n');

%% Figure 7: Parallel Coordinates by Cluster
fprintf('Creating Figure 7: Parallel Coordinates by Cluster...\n');

nrows = ceil(optimal_k/2);
ncols = 2;
figure('Position', [220, 220, 1400, 400*nrows]);

for k = 1:optimal_k
    subplot(nrows, ncols, k);
    hold on;
    
    idx_k = cluster_idx == k;
    
    % Plot all salts in gray in background
    for i = 1:n_salts
        plot(1:4, features_norm(i, :), '-', 'LineWidth', 0.5, ...
            'Color', [0.8, 0.8, 0.8, 0.3], 'HandleVisibility', 'off');
    end
    
    % Plot cluster salts in color
    salts_in_cluster = find(idx_k);
    for i = salts_in_cluster'
        plot(1:4, features_norm(i, :), '-o', 'LineWidth', 2, ...
            'Color', cmap_clusters(k, :), 'MarkerSize', 6, ...
            'MarkerFaceColor', cmap_clusters(k, :), 'HandleVisibility', 'off');
    end
    
    % Plot cluster centroid
    plot(1:4, centroids(k, :), 'k-', 'LineWidth', 3, 'DisplayName', 'Centroid');
    
    set(gca, 'XTick', 1:4, 'XTickLabel', feature_names);
    xlim([0.5, 4.5]);
    ylim([0, 1]);
    ylabel('Normalized Value', 'FontSize', 10);
    title(sprintf('Cluster %d (%d salts)', k, sum(idx_k)), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 9);
    if k == 1
        legend('Location', 'best');
    end
end

% Add main title
sgtitle('Parallel Coordinates by Cluster', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, '../../../figures/curve_reduction/parallel_coordinates_by_cluster.png');
saveas(gcf, '../../../figures/curve_reduction/parallel_coordinates_by_cluster.fig');
fprintf('  Saved: parallel_coordinates_by_cluster.png\n');

%% Figure 8: Cluster Characteristics
fprintf('Creating Figure 8: Cluster Centroids...\n');
figure('Position', [240, 240, 1200, 800]);

% Create bar plot of cluster centroids
x_pos = 1:size(centroids, 2);
bar_width = 0.8 / optimal_k;

hold on;
for k = 1:optimal_k
    offset = (k - (optimal_k+1)/2) * bar_width;
    bar(x_pos + offset, centroids(k, :), bar_width, ...
        'FaceColor', cmap_clusters(k, :), 'EdgeColor', 'k', 'LineWidth', 1);
end

set(gca, 'XTick', 1:length(feature_names), 'XTickLabel', feature_names);
ylabel('Normalized Value', 'FontSize', 12);
title('Cluster Centroids', 'FontSize', 14, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% Save figure
saveas(gcf, '../../../figures/curve_reduction/cluster_centroids.png');
saveas(gcf, '../../../figures/curve_reduction/cluster_centroids.fig');
fprintf('  Saved: cluster_centroids.png\n');

%% Summary Statistics
fprintf('\n=== Cluster Summary Statistics ===\n');
for k = 1:optimal_k
    idx_k = cluster_idx == k;
    fprintf('\nCluster %d:\n', k);
    fprintf('  Mean Efficiency: %.4f\n', mean(efficiency(idx_k)));
    fprintf('  Mean Capacity: %.4f\n', mean(capacity(idx_k)));
    fprintf('  Mean Slope: %.4f\n', mean(slope(idx_k)));
    fprintf('  Mean Curvature: %.4f\n', mean(curvature(idx_k)));
    fprintf('  Mean Max Molality: %.4f\n', mean(max_molality(idx_k)));
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('All figures saved to ../../../figures/curve_reduction/\n');
