%% Curve Reduction: Extract Interpretable Scalar Features from Water Activity Curves
% This script reduces full aw(m) curves to physically meaningful scalar features
% for clustering, regression, and ranking salts for AWH uptake.
%
% Features extracted:
%   1. Initial slope: d(ln aw)/dm at m→0
%   2. Curvature: d²(ln aw)/dm²  
%   3. Threshold molality where aw < 0.75 (AWH target)
%   4. Integral: ∫ln(aw) dm (total water binding strength)
%
% Author: Generated for AWH ML Discussion
% Date: 2026-01-25

clear; clc; close all;

%% Configuration
RH_THRESHOLD = 0.75;  % Target RH for AWH applications
MIN_POINTS = 5;       % Minimum data points required for fitting
OUTPUT_DIR = '../../../figures/curve_reduction/';
DATA_FILE = '../../data/water_activity_all_salts_combined.csv';

% Create output directory if it doesn't exist
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('=== Water Activity Curve Reduction ===\n\n');

%% Load Data
fprintf('Loading water activity data...\n');
data = readtable(DATA_FILE);

% Get unique salts
salts = unique(data.Salt);
n_salts = length(salts);
fprintf('Found %d unique salts\n\n', n_salts);

%% Initialize Results Structure
results = struct();
feature_names = {'Salt', 'Initial_Slope_dln_aw_dm', 'Curvature_d2ln_aw_dm2', ...
                'Threshold_Molality_aw_075', 'Integral_ln_aw', ...
                'Max_Molality', 'N_Data_Points', 'Fit_RMSE'};

results_table = cell(n_salts, length(feature_names));

%% Process Each Salt
fprintf('Processing salts:\n');
fprintf('--------------------------------------------------\n');

for i = 1:n_salts
    salt_name = salts{i};
    fprintf('[%d/%d] %s ... ', i, n_salts, salt_name);
    
    % Extract data for this salt
    idx = strcmp(data.Salt, salt_name);
    salt_data = data(idx, :);
    
    % Sort by molality
    [m_sorted, sort_idx] = sort(salt_data.Molality_mol_per_kg);
    aw_sorted = salt_data.RH_Water_Activity(sort_idx);
    ln_aw_sorted = log(aw_sorted);
    
    % Store basic info
    results_table{i, 1} = salt_name;
    results_table{i, 6} = max(m_sorted);
    results_table{i, 7} = length(m_sorted);
    
    % Check if we have enough points
    if length(m_sorted) < MIN_POINTS
        fprintf('SKIP (insufficient data: %d points)\n', length(m_sorted));
        results_table{i, 2:5} = {NaN, NaN, NaN, NaN};
        results_table{i, 8} = NaN;
        continue;
    end
    
    %% Fit Smooth Monotone Model
    % Use piecewise cubic Hermite interpolating polynomial (PCHIP)
    % which preserves monotonicity and shape of the data
    try
        % Create fine grid for evaluation
        m_fine = linspace(min(m_sorted), max(m_sorted), 500)';
        
        % PCHIP interpolation (shape-preserving, monotone)
        ln_aw_fit = pchip(m_sorted, ln_aw_sorted, m_fine);
        
        % Calculate RMSE on original data points
        ln_aw_fit_orig = pchip(m_sorted, ln_aw_sorted, m_sorted);
        rmse = sqrt(mean((ln_aw_sorted - ln_aw_fit_orig).^2));
        results_table{i, 8} = rmse;
        
        %% Feature 1: Initial Slope at m→0
        % Use finite difference on the fitted curve near m=0
        if min(m_sorted) < 0.1
            % If we have data near zero, use that region
            idx_near_zero = find(m_fine < 0.5, 20, 'first');
        else
            % Otherwise use the first few points
            idx_near_zero = 1:min(20, length(m_fine));
        end
        
        if length(idx_near_zero) >= 2
            % Linear fit to initial region
            p_init = polyfit(m_fine(idx_near_zero), ln_aw_fit(idx_near_zero), 1);
            initial_slope = p_init(1);
        else
            initial_slope = NaN;
        end
        results_table{i, 2} = initial_slope;
        
        %% Feature 2: Average Curvature
        % Calculate numerical second derivative
        dm = m_fine(2) - m_fine(1);
        d2ln_aw = diff(ln_aw_fit, 2) / (dm^2);
        
        % Use mean absolute curvature across the full range
        avg_curvature = mean(abs(d2ln_aw));
        results_table{i, 3} = avg_curvature;
        
        %% Feature 3: Threshold Molality where aw < RH_THRESHOLD
        % Convert back to aw space
        aw_fit = exp(ln_aw_fit);
        
        % Find where aw crosses threshold
        idx_below = find(aw_fit < RH_THRESHOLD, 1, 'first');
        if ~isempty(idx_below) && idx_below > 1
            % Interpolate to get exact crossing point
            m1 = m_fine(idx_below-1);
            m2 = m_fine(idx_below);
            aw1 = aw_fit(idx_below-1);
            aw2 = aw_fit(idx_below);
            
            % Linear interpolation
            threshold_molality = m1 + (m2-m1) * (RH_THRESHOLD - aw1) / (aw2 - aw1);
        elseif min(aw_fit) > RH_THRESHOLD
            % Never reaches threshold - use max molality + marker
            threshold_molality = Inf;
        else
            % Starts below threshold
            threshold_molality = min(m_fine);
        end
        results_table{i, 4} = threshold_molality;
        
        %% Feature 4: Integral of ln(aw) over molality range
        % ∫ln(aw) dm from 0 (or min m) to max m
        % This represents total "water binding strength"
        integral_ln_aw = trapz(m_fine, ln_aw_fit);
        results_table{i, 5} = integral_ln_aw;
        
        fprintf('OK (RMSE=%.4f)\n', rmse);
        
        %% Optional: Create diagnostic plots for selected salts
        % Plot first 6 salts and a few interesting ones
        if i <= 6 || ismember(salt_name, {'LiCl', 'CaCl2', 'LiBr', 'MgCl2'})
            fig = figure('Visible', 'off');
            set(fig, 'Position', [100, 100, 1200, 400]);
            
            % Plot 1: Data and fit
            subplot(1, 3, 1);
            plot(m_sorted, ln_aw_sorted, 'ko', 'MarkerSize', 6, 'LineWidth', 1.5); hold on;
            plot(m_fine, ln_aw_fit, 'b-', 'LineWidth', 2);
            xlabel('Molality (mol/kg)');
            ylabel('ln(a_w)');
            title(sprintf('%s: Fitted Curve', strrep(salt_name, '_', '\_')));
            legend('Data', 'PCHIP Fit', 'Location', 'best');
            grid on;
            
            % Plot 2: First derivative
            subplot(1, 3, 2);
            dln_aw = diff(ln_aw_fit) / dm;
            plot(m_fine(1:end-1), dln_aw, 'r-', 'LineWidth', 2);
            hold on;
            % Mark initial slope
            if ~isnan(initial_slope)
                plot([0, m_fine(20)], initial_slope*[0, m_fine(20)], 'g--', 'LineWidth', 2);
            end
            xlabel('Molality (mol/kg)');
            ylabel('d(ln a_w) / dm');
            title('First Derivative (Slope)');
            legend('Derivative', 'Initial Slope', 'Location', 'best');
            grid on;
            
            % Plot 3: Threshold visualization
            subplot(1, 3, 3);
            plot(m_fine, aw_fit, 'b-', 'LineWidth', 2); hold on;
            yline(RH_THRESHOLD, 'r--', 'LineWidth', 2, 'Label', sprintf('RH = %.2f', RH_THRESHOLD));
            if isfinite(threshold_molality)
                plot(threshold_molality, RH_THRESHOLD, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
                text(threshold_molality, RH_THRESHOLD-0.05, ...
                     sprintf('m = %.2f', threshold_molality), ...
                     'HorizontalAlignment', 'center');
            end
            xlabel('Molality (mol/kg)');
            ylabel('Water Activity (a_w)');
            title('Threshold Detection');
            grid on;
            ylim([0, 1]);
            
            % Save figure
            saveas(fig, fullfile(OUTPUT_DIR, sprintf('curve_fit_%s.png', salt_name)));
            close(fig);
        end
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        results_table{i, 2:5} = {NaN, NaN, NaN, NaN};
        results_table{i, 8} = NaN;
    end
end

fprintf('--------------------------------------------------\n\n');

%% Convert to Table and Save
fprintf('Creating results table...\n');
results_final = cell2table(results_table, 'VariableNames', feature_names);

% Sort by integral (water binding strength) - lower is better for AWH
results_final = sortrows(results_final, 'Integral_ln_aw', 'ascend');

% Save to CSV
output_file = fullfile(OUTPUT_DIR, 'curve_features.csv');
writetable(results_final, output_file);
fprintf('Saved results to: %s\n\n', output_file);

%% Display Summary Statistics
fprintf('=== Feature Summary Statistics ===\n');
fprintf('Initial Slope (d ln_aw/dm at m→0):\n');
initial_slopes = results_final.Initial_Slope_dln_aw_dm;
if iscell(initial_slopes), initial_slopes = cell2mat(initial_slopes); end
initial_slopes_clean = initial_slopes(~isnan(initial_slopes));
fprintf('  Mean: %.4f, Std: %.4f\n', ...
        mean(initial_slopes_clean), std(initial_slopes_clean));

fprintf('\nCurvature (mean |d²ln_aw/dm²|):\n');
curvatures = results_final.Curvature_d2ln_aw_dm2;
if iscell(curvatures), curvatures = cell2mat(curvatures); end
curvatures_clean = curvatures(~isnan(curvatures));
fprintf('  Mean: %.4f, Std: %.4f\n', ...
        mean(curvatures_clean), std(curvatures_clean));

fprintf('\nThreshold Molality (aw < 0.75):\n');
threshold_vals = results_final.Threshold_Molality_aw_075;
if iscell(threshold_vals), threshold_vals = cell2mat(threshold_vals); end
threshold_finite = threshold_vals(isfinite(threshold_vals));
fprintf('  Mean (finite): %.4f, Std: %.4f\n', mean(threshold_finite), std(threshold_finite));
fprintf('  %d salts never reach threshold\n', sum(isinf(threshold_vals)));

fprintf('\nIntegral (∫ln_aw dm):\n');
integrals = results_final.Integral_ln_aw;
if iscell(integrals), integrals = cell2mat(integrals); end
integrals_clean = integrals(~isnan(integrals));
fprintf('  Mean: %.4f, Std: %.4f\n', ...
        mean(integrals_clean), std(integrals_clean));

fprintf('\nFit Quality (RMSE):\n');
rmses = results_final.Fit_RMSE;
if iscell(rmses), rmses = cell2mat(rmses); end
rmses_clean = rmses(~isnan(rmses));
fprintf('  Mean: %.6f, Std: %.6f\n\n', ...
        mean(rmses_clean), std(rmses_clean));

%% Top 10 Salts by Water Binding Strength
fprintf('=== Top 10 Salts by Water Binding Strength (lowest ∫ln_aw dm) ===\n');
fprintf('%-12s %12s %12s %12s %12s\n', 'Salt', 'Init_Slope', 'Curvature', 'Threshold_m', 'Integral');
fprintf('------------------------------------------------------------------------\n');
for i = 1:min(10, height(results_final))
    salt_name_i = results_final.Salt(i);
    if iscell(salt_name_i), salt_name_i = salt_name_i{1}; end
    
    init_slope_i = results_final.Initial_Slope_dln_aw_dm(i);
    if iscell(init_slope_i), init_slope_i = init_slope_i{1}; end
    
    curvature_i = results_final.Curvature_d2ln_aw_dm2(i);
    if iscell(curvature_i), curvature_i = curvature_i{1}; end
    
    threshold_i = results_final.Threshold_Molality_aw_075(i);
    if iscell(threshold_i), threshold_i = threshold_i{1}; end
    
    integral_i = results_final.Integral_ln_aw(i);
    if iscell(integral_i), integral_i = integral_i{1}; end
    
    fprintf('%-12s %12.4f %12.6f %12.4f %12.4f\n', ...
            salt_name_i, init_slope_i, curvature_i, threshold_i, integral_i);
end
fprintf('\n');

%% Create Summary Visualizations
fprintf('Creating summary visualizations...\n');

% Figure 1: Feature distributions
fig = figure('Visible', 'off');
set(fig, 'Position', [100, 100, 1400, 900]);

subplot(2, 3, 1);
slopes_plot = initial_slopes(~isnan(initial_slopes));
histogram(slopes_plot, 20);
xlabel('Initial Slope d(ln a_w)/dm');
ylabel('Count');
title('Distribution: Initial Slope');
grid on;

subplot(2, 3, 2);
curvatures_plot = curvatures(~isnan(curvatures));
histogram(curvatures_plot, 20);
xlabel('Mean |Curvature|');
ylabel('Count');
title('Distribution: Curvature');
grid on;

subplot(2, 3, 3);
histogram(threshold_finite, 20);
xlabel('Threshold Molality (mol/kg)');
ylabel('Count');
title('Distribution: Threshold (aw < 0.75)');
grid on;

subplot(2, 3, 4);
integrals_plot = integrals(~isnan(integrals));
histogram(integrals_plot, 20);
xlabel('∫ln(a_w) dm');
ylabel('Count');
title('Distribution: Integral (Water Binding)');
grid on;

subplot(2, 3, 5);
scatter(initial_slopes, integrals, 50, 'filled');
xlabel('Initial Slope');
ylabel('Integral ln(a_w)');
title('Initial Slope vs Integral');
grid on;

subplot(2, 3, 6);
threshold_finite_idx = isfinite(threshold_vals);
scatter(threshold_vals(threshold_finite_idx), ...
        integrals(threshold_finite_idx), 50, 'filled');
xlabel('Threshold Molality');
ylabel('Integral ln(a_w)');
title('Threshold vs Integral');
grid on;

saveas(fig, fullfile(OUTPUT_DIR, 'feature_distributions.png'));
close(fig);

fprintf('Saved summary visualizations\n');

%% Feature Correlation Analysis
fprintf('\n=== Feature Correlation Matrix ===\n');
feature_matrix = [initial_slopes, curvatures, threshold_vals, integrals];

% Remove rows with NaN or Inf
valid_rows = all(isfinite(feature_matrix), 2);
feature_matrix_clean = feature_matrix(valid_rows, :);

if size(feature_matrix_clean, 1) > 2
    % Manual correlation calculation (no Statistics Toolbox required)
    n_features = size(feature_matrix_clean, 2);
    corr_matrix = zeros(n_features, n_features);
    
    for ii = 1:n_features
        for jj = 1:n_features
            x = feature_matrix_clean(:, ii);
            y = feature_matrix_clean(:, jj);
            
            % Correlation formula: r = cov(x,y) / (std(x) * std(y))
            x_centered = x - mean(x);
            y_centered = y - mean(y);
            
            corr_matrix(ii, jj) = sum(x_centered .* y_centered) / ...
                                  sqrt(sum(x_centered.^2) * sum(y_centered.^2));
        end
    end
    
    fprintf('\n');
    fprintf('                Init_Slope  Curvature  Threshold   Integral\n');
    fprintf('Init_Slope      %8.3f   %8.3f  %8.3f  %8.3f\n', corr_matrix(1,:));
    fprintf('Curvature       %8.3f   %8.3f  %8.3f  %8.3f\n', corr_matrix(2,:));
    fprintf('Threshold       %8.3f   %8.3f  %8.3f  %8.3f\n', corr_matrix(3,:));
    fprintf('Integral        %8.3f   %8.3f  %8.3f  %8.3f\n', corr_matrix(4,:));
    
    % Plot correlation matrix
    fig = figure('Visible', 'off');
    imagesc(corr_matrix);
    colorbar;
    caxis([-1, 1]);
    colormap(redblue(256));
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Init_Slope', 'Curvature', 'Threshold', 'Integral'});
    set(gca, 'YTick', 1:4, 'YTickLabel', {'Init_Slope', 'Curvature', 'Threshold', 'Integral'});
    title('Feature Correlation Matrix');
    
    % Add correlation values as text
    for ii = 1:4
        for jj = 1:4
            text(jj, ii, sprintf('%.2f', corr_matrix(ii,jj)), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
    
    saveas(fig, fullfile(OUTPUT_DIR, 'feature_correlations.png'));
    close(fig);
end

fprintf('\n=== Curve Reduction Complete ===\n');
fprintf('Results saved to: %s\n', OUTPUT_DIR);
fprintf('Key output file: curve_features.csv\n\n');

%% Helper function for red-blue colormap
function cmap = redblue(n)
    if nargin < 1
        n = 256;
    end
    r = [(0:n/2-1)'/(n/2-1); ones(n/2,1)];
    g = [(0:n/2-1)'/(n/2-1); (n/2-1:-1:0)'/(n/2-1)];
    b = [ones(n/2,1); (n/2-1:-1:0)'/(n/2-1)];
    cmap = [r g b];
end
