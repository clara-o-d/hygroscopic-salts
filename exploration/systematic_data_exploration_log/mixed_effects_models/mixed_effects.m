%% Mixed-Effects Model for Water Activity: Hierarchical Modeling of Ion Effects
% This script implements a hierarchical (mixed-effects) model to decompose
% water activity into fixed concentration effects and random ion effects.
%
% Model: ln(aw_ij) = f(m_ij) + u_cation(i) + v_anion(i) + ε_ij
%
% Where:
%   - f(m) is a fixed effect (concentration dependence)
%   - u_cation and v_anion are random effects (ion-specific deviations)
%   - ε is residual error
%
% Goals:
%   1. Quantify variance from cations vs anions
%   2. Identify consistent ion effects across different pairings
%   3. Enable prediction for unseen salt combinations
%
% Author: Generated for AWH ML Discussion
% Date: 2026-01-25

clear; clc; close all;

%% Configuration
DATA_FILE = '../../../data/water_activity_all_salts_combined.csv';
OUTPUT_DIR = '../../../figures/mixed_effects/';
MIN_POINTS_PER_SALT = 5;
MOLALITY_GRID = linspace(0, 10, 100)';  % Grid for fixed effect estimation
GENERATE_FIGURES = true;  % Set to false if no graphics toolkit available

% Create output directory
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('=== Mixed-Effects Model for Water Activity ===\n\n');

%% Load and Prepare Data
fprintf('Loading data...\n');

% Read CSV file
fid = fopen(DATA_FILE, 'r');
header = fgetl(fid);
fclose(fid);

% Parse header to find column indices
header_parts = strsplit(header, ',');
salt_col = find(strcmp(header_parts, 'Salt'));
molality_col = find(strcmp(header_parts, 'Molality_mol_per_kg'));
aw_col = find(strcmp(header_parts, 'RH_Water_Activity'));

% Read data using dlmread (skip header)
raw_data = {};
fid = fopen(DATA_FILE, 'r');
fgetl(fid); % Skip header

line_count = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        line_count = line_count + 1;
        parts = strsplit(line, ',');
        raw_data{line_count, 1} = parts{salt_col};
        raw_data{line_count, 2} = str2double(parts{molality_col});
        raw_data{line_count, 3} = str2double(parts{aw_col});
    end
end
fclose(fid);

% Convert to structure
data = struct();
data.Salt = raw_data(:, 1);
data.Molality = cell2mat(raw_data(:, 2));
data.aw = cell2mat(raw_data(:, 3));
data.ln_aw = log(data.aw);

n_total = length(data.Salt);

% Filter salts with sufficient data
salts = unique(data.Salt);
valid_salts = {};
for i = 1:length(salts)
    if sum(strcmp(data.Salt, salts{i})) >= MIN_POINTS_PER_SALT
        valid_salts{end+1} = salts{i};
    end
end

% Filter data
valid_idx = false(n_total, 1);
for i = 1:length(valid_salts)
    valid_idx = valid_idx | strcmp(data.Salt, valid_salts{i});
end

data.Salt = data.Salt(valid_idx);
data.Molality = data.Molality(valid_idx);
data.aw = data.aw(valid_idx);
data.ln_aw = data.ln_aw(valid_idx);

n_obs = length(data.Salt);

fprintf('Total observations: %d\n', n_obs);
fprintf('Unique salts: %d\n', length(valid_salts));

%% Parse Salt Names to Extract Ions
fprintf('\nParsing salt formulas to extract cations and anions...\n');
[cations, anions] = parse_salt_formulas(data.Salt);

% Add to data table
data.Cation = cations;
data.Anion = anions;

% Get unique ions
unique_cations = unique(cations);
unique_anions = unique(anions);

fprintf('Unique cations: %d\n', length(unique_cations));
fprintf('Unique anions: %d\n', length(unique_anions));

% Display ion statistics
fprintf('\nCations: ');
fprintf('%s ', unique_cations{:});
fprintf('\n');
fprintf('Anions: ');
fprintf('%s ', unique_anions{:});
fprintf('\n');

%% Create Ion Encoding
fprintf('\nCreating ion encoding matrices...\n');

% Create dummy variables for cations and anions
n_cations = length(unique_cations);
n_anions = length(unique_anions);

% Cation matrix (n_obs × n_cations)
X_cation = zeros(n_obs, n_cations);
for i = 1:n_cations
    X_cation(:, i) = strcmp(data.Cation, unique_cations{i});
end

% Anion matrix (n_obs × n_anions)
X_anion = zeros(n_obs, n_anions);
for i = 1:n_anions
    X_anion(:, i) = strcmp(data.Anion, unique_anions{i});
end

%% Step 1: Estimate Fixed Effect f(m) - Average Concentration Dependence
fprintf('\n=== Step 1: Estimating Fixed Effect f(m) ===\n');

% Use local regression (loess) on all data to get average trend
% This captures the typical ln(aw) vs molality relationship
m_all = data.Molality;
ln_aw_all = data.ln_aw;

% Bin data for smoother estimation
m_bins = 0:0.5:ceil(max(m_all));
f_m_bins = zeros(length(m_bins), 1);
f_m_std = zeros(length(m_bins), 1);

for i = 1:length(m_bins)
    if i < length(m_bins)
        idx = (m_all >= m_bins(i)) & (m_all < m_bins(i+1));
    else
        idx = m_all >= m_bins(i);
    end
    
    if sum(idx) > 0
        f_m_bins(i) = median(ln_aw_all(idx));
        f_m_std(i) = std(ln_aw_all(idx));
    else
        f_m_bins(i) = NaN;
        f_m_std(i) = NaN;
    end
end

% Interpolate to create smooth f(m)
valid_bins = ~isnan(f_m_bins);
m_bins_valid = m_bins(valid_bins);
f_m_bins_valid = f_m_bins(valid_bins);

% Use linear interpolation with no extrapolation (returns NaN for out-of-range)
% Then fill NaN with nearest neighbor
f_m_values = interp1(m_bins_valid, f_m_bins_valid, data.Molality, 'pchip');

% For any NaN values (extrapolation), use nearest neighbor
nan_idx = isnan(f_m_values);
if any(nan_idx)
    fprintf('Warning: %d points outside interpolation range, using nearest neighbor\n', sum(nan_idx));
    f_m_values(nan_idx) = interp1(m_bins_valid, f_m_bins_valid, data.Molality(nan_idx), 'nearest', 'extrap');
end

% Create interpolation function for later use
f_m_interp = @(m) interp1_with_nearest(m_bins_valid, f_m_bins_valid, m);
data.f_m = f_m_values;

% Initialize other data fields
data.residual = zeros(n_obs, 1);
data.ln_aw_pred = zeros(n_obs, 1);
data.residual_final = zeros(n_obs, 1);

fprintf('Fixed effect f(m) estimated using binned medians\n');
fprintf('Concentration range: %.2f - %.2f mol/kg\n', min(m_all), max(m_all));

%% Step 2: Estimate Random Effects (Ion Deviations)
fprintf('\n=== Step 2: Estimating Random Ion Effects ===\n');

% Residuals after removing fixed effect
for i = 1:n_obs
    data.residual(i) = data.ln_aw(i) - data.f_m(i);
end

% Estimate cation effects (marginal means)
cation_effects = zeros(n_cations, 1);
cation_std = zeros(n_cations, 1);
cation_n = zeros(n_cations, 1);

for i = 1:n_cations
    idx = false(n_obs, 1);
    for j = 1:n_obs
        if strcmp(data.Cation{j}, unique_cations{i})
            idx(j) = true;
        end
    end
    residual_subset = data.residual(idx);
    cation_effects(i) = mean(residual_subset);
    cation_std(i) = std(residual_subset);
    cation_n(i) = sum(idx);
end

% Estimate anion effects (marginal means)
anion_effects = zeros(n_anions, 1);
anion_std = zeros(n_anions, 1);
anion_n = zeros(n_anions, 1);

for i = 1:n_anions
    idx = false(n_obs, 1);
    for j = 1:n_obs
        if strcmp(data.Anion{j}, unique_anions{i})
            idx(j) = true;
        end
    end
    residual_subset = data.residual(idx);
    anion_effects(i) = mean(residual_subset);
    anion_std(i) = std(residual_subset);
    anion_n(i) = sum(idx);
end

% Re-center to sum to zero (identifiability constraint)
cation_effects = cation_effects - mean(cation_effects);
anion_effects = anion_effects - mean(anion_effects);

fprintf('Cation effects estimated (centered)\n');
fprintf('Anion effects estimated (centered)\n');

%% Step 3: Mixed-Effects Model Prediction
fprintf('\n=== Step 3: Model Predictions ===\n');

% Predict ln(aw) using mixed model
for i = 1:n_obs
    data.ln_aw_pred(i) = data.f_m(i);
    
    % Add cation effect
    for c = 1:n_cations
        if strcmp(data.Cation{i}, unique_cations{c})
            data.ln_aw_pred(i) = data.ln_aw_pred(i) + cation_effects(c);
            break;
        end
    end
    
    % Add anion effect
    for a = 1:n_anions
        if strcmp(data.Anion{i}, unique_anions{a})
            data.ln_aw_pred(i) = data.ln_aw_pred(i) + anion_effects(a);
            break;
        end
    end
end

% Calculate model fit
for i = 1:n_obs
    data.residual_final(i) = data.ln_aw(i) - data.ln_aw_pred(i);
end
rmse_total = sqrt(mean(data.residual_final.^2));
r2_total = 1 - sum(data.residual_final.^2) / sum((data.ln_aw - mean(data.ln_aw)).^2);

fprintf('\nModel Performance:\n');
fprintf('  RMSE: %.4f\n', rmse_total);
fprintf('  R²: %.4f\n', r2_total);

%% Step 4: Variance Decomposition
fprintf('\n=== Step 4: Variance Decomposition ===\n');

% Total variance in ln(aw)
var_total = var(data.ln_aw);

% Variance explained by fixed effect
% Check for NaN values first
if any(isnan(data.f_m))
    fprintf('Warning: NaN values found in f_m\n');
    fprintf('Number of NaN: %d\n', sum(isnan(data.f_m)));
end
var_fixed = var(data.f_m(~isnan(data.f_m)));

% Variance from cation effects (weighted by occurrence)
cation_var_contributions = zeros(n_obs, 1);
for i = 1:n_obs
    for c = 1:n_cations
        if strcmp(data.Cation{i}, unique_cations{c})
            cation_var_contributions(i) = cation_effects(c);
            break;
        end
    end
end
var_cation = var(cation_var_contributions);

% Variance from anion effects
anion_var_contributions = zeros(n_obs, 1);
for i = 1:n_obs
    for a = 1:n_anions
        if strcmp(data.Anion{i}, unique_anions{a})
            anion_var_contributions(i) = anion_effects(a);
            break;
        end
    end
end
var_anion = var(anion_var_contributions);

% Residual variance
var_residual = var(data.residual_final);

% Compute percentages
pct_fixed = 100 * var_fixed / var_total;
pct_cation = 100 * var_cation / var_total;
pct_anion = 100 * var_anion / var_total;
pct_residual = 100 * var_residual / var_total;

fprintf('\nVariance Decomposition:\n');
fprintf('  Total variance: %.4f\n', var_total);
fprintf('  Fixed effect (concentration): %.4f (%.1f%%)\n', var_fixed, pct_fixed);
fprintf('  Random effect (cation): %.4f (%.1f%%)\n', var_cation, pct_cation);
fprintf('  Random effect (anion): %.4f (%.1f%%)\n', var_anion, pct_anion);
fprintf('  Residual: %.4f (%.1f%%)\n', var_residual, pct_residual);

% Ratio of cation to anion variance
cation_anion_ratio = var_cation / var_anion;
fprintf('\nCation:Anion variance ratio: %.2f\n', cation_anion_ratio);
if cation_anion_ratio > 1
    fprintf('  → Cation identity explains %.1fx more variance than anion\n', cation_anion_ratio);
else
    fprintf('  → Anion identity explains %.1fx more variance than cation\n', 1/cation_anion_ratio);
end

%% Step 5: Cross-Validation - Predict Unseen Salts
fprintf('\n=== Step 5: Cross-Validation for Unseen Salts ===\n');

% Leave-one-salt-out cross-validation
unique_salts = unique(data.Salt);
n_salts = length(unique_salts);

cv_predictions = cell(n_salts, 1);
cv_rmse = zeros(n_salts, 1);
cv_r2 = zeros(n_salts, 1);

fprintf('Performing leave-one-salt-out CV...\n');
fprintf('Progress: ');

for s = 1:n_salts
    if mod(s, 10) == 0 || s == n_salts
        fprintf('%d/%d ', s, n_salts);
    end
    
    % Split data
    test_idx = false(n_obs, 1);
    for i = 1:n_obs
        if strcmp(data.Salt{i}, unique_salts{s})
            test_idx(i) = true;
        end
    end
    train_idx = ~test_idx;
    
    test_salt = unique_salts{s};
    
    % Find first test observation to get cation/anion
    first_test = find(test_idx, 1);
    test_cation = data.Cation{first_test};
    test_anion = data.Anion{first_test};
    
    % Re-estimate fixed effect on training data
    m_train = data.Molality(train_idx);
    ln_aw_train = data.ln_aw(train_idx);
    
    % Bin training data
    f_m_bins_train = zeros(length(m_bins), 1);
    for i = 1:length(m_bins)
        if i < length(m_bins)
            idx = (m_train >= m_bins(i)) & (m_train < m_bins(i+1));
        else
            idx = m_train >= m_bins(i);
        end
        
        if sum(idx) > 0
            f_m_bins_train(i) = median(ln_aw_train(idx));
        else
            f_m_bins_train(i) = NaN;
        end
    end
    
    valid_bins_train = ~isnan(f_m_bins_train);
    m_bins_train_valid = m_bins(valid_bins_train);
    f_m_bins_train_valid = f_m_bins_train(valid_bins_train);
    
    % Re-estimate random effects on training data
    f_m_train_values = interp1_with_nearest(m_bins_train_valid, f_m_bins_train_valid, m_train);
    residual_train = ln_aw_train - f_m_train_values;
    
    % Get training cations and anions
    train_cations = data.Cation(train_idx);
    train_anions = data.Anion(train_idx);
    
    % Cation effects (training)
    cation_effects_train = zeros(n_cations, 1);
    for c = 1:n_cations
        idx = false(length(train_cations), 1);
        for j = 1:length(train_cations)
            if strcmp(train_cations{j}, unique_cations{c})
                idx(j) = true;
            end
        end
        if sum(idx) > 0
            cation_effects_train(c) = mean(residual_train(idx));
        end
    end
    cation_effects_train = cation_effects_train - mean(cation_effects_train);
    
    % Anion effects (training)
    anion_effects_train = zeros(n_anions, 1);
    for a = 1:n_anions
        idx = false(length(train_anions), 1);
        for j = 1:length(train_anions)
            if strcmp(train_anions{j}, unique_anions{a})
                idx(j) = true;
            end
        end
        if sum(idx) > 0
            anion_effects_train(a) = mean(residual_train(idx));
        end
    end
    anion_effects_train = anion_effects_train - mean(anion_effects_train);
    
    % Predict test salt
    m_test = data.Molality(test_idx);
    ln_aw_test = data.ln_aw(test_idx);
    
    f_m_test = interp1_with_nearest(m_bins_train_valid, f_m_bins_train_valid, m_test);
    
    cation_idx = find(strcmp(unique_cations, test_cation));
    anion_idx = find(strcmp(unique_anions, test_anion));
    
    ln_aw_pred_test = f_m_test + cation_effects_train(cation_idx) + anion_effects_train(anion_idx);
    
    % Store results
    cv_predictions{s} = struct('Salt', test_salt, ...
                               'Cation', test_cation, ...
                               'Anion', test_anion, ...
                               'Molality', m_test, ...
                               'ln_aw_true', ln_aw_test, ...
                               'ln_aw_pred', ln_aw_pred_test);
    
    % Calculate metrics
    cv_rmse(s) = sqrt(mean((ln_aw_test - ln_aw_pred_test).^2));
    cv_r2(s) = 1 - sum((ln_aw_test - ln_aw_pred_test).^2) / sum((ln_aw_test - mean(ln_aw_test)).^2);
end

fprintf('\n');

fprintf('\nCross-Validation Results:\n');
fprintf('  Mean RMSE: %.4f ± %.4f\n', mean(cv_rmse), std(cv_rmse));
fprintf('  Mean R²: %.4f ± %.4f\n', mean(cv_r2), std(cv_r2));
fprintf('  Median RMSE: %.4f\n', median(cv_rmse));
fprintf('  Median R²: %.4f\n', median(cv_r2));

% Best and worst predictions
[~, best_idx] = min(cv_rmse);
[~, worst_idx] = max(cv_rmse);

fprintf('\nBest prediction: %s (RMSE=%.4f, R²=%.4f)\n', ...
    unique_salts{best_idx}, cv_rmse(best_idx), cv_r2(best_idx));
fprintf('Worst prediction: %s (RMSE=%.4f, R²=%.4f)\n', ...
    unique_salts{worst_idx}, cv_rmse(worst_idx), cv_r2(worst_idx));

%% Visualization
fprintf('\n=== Generating Visualizations ===\n');

% Check if graphics are available
graphics_available = true;
try
    figure('visible', 'off');
    close;
catch
    graphics_available = false;
    fprintf('Warning: No graphics toolkit available. Skipping figure generation.\n');
end

if ~graphics_available || ~GENERATE_FIGURES
    fprintf('Skipping figure generation (graphics not available or disabled)\n');
else

%% Figure 1: Variance Decomposition Pie Chart
try
fig1 = figure('Position', [100, 100, 800, 600], 'visible', 'off');
variance_components = [var_fixed, var_cation, var_anion, var_residual];
labels = {'Fixed Effect (Concentration)', 'Cation Identity', 'Anion Identity', 'Residual'};
colors = [0.2 0.4 0.8; 0.8 0.3 0.3; 0.3 0.7 0.3; 0.7 0.7 0.7];

pie(variance_components);
colormap(colors);
legend(labels, 'Location', 'eastoutside', 'FontSize', 11);
title('Variance Decomposition in ln(a_w)', 'FontSize', 14, 'FontWeight', 'bold');

% Add percentage labels
text_objs = findobj(fig1, 'Type', 'text');
for i = 1:length(text_objs)
    curr_text = text_objs(i).String;
    if contains(curr_text, '%')
        text_objs(i).FontSize = 12;
        text_objs(i).FontWeight = 'bold';
    end
end

saveas(fig1, fullfile(OUTPUT_DIR, 'variance_decomposition.png'));
fprintf('Saved: variance_decomposition.png\n');
catch e
    fprintf('Error with figure 1: %s\n', e.message);
end

%% Figure 2: Ion Effects
try
fig2 = figure('Position', [100, 100, 1200, 500], 'visible', 'off');

% Sort ions by effect size
[cation_effects_sorted, cation_sort_idx] = sort(cation_effects);
[anion_effects_sorted, anion_sort_idx] = sort(anion_effects);

subplot(1, 2, 1);
bar(cation_effects_sorted, 'FaceColor', [0.8 0.3 0.3]);
hold on;
% Error bars
errorbar(1:n_cations, cation_effects_sorted, cation_std(cation_sort_idx) / sqrt(mean(cation_n)), ...
    'k.', 'LineWidth', 1.5);
set(gca, 'XTick', 1:n_cations, 'XTickLabel', unique_cations(cation_sort_idx), ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Random Effect on ln(a_w)', 'FontSize', 12);
title('Cation Effects (u_{cation})', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xl = xlim;
line([xl(1) xl(2)], [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

subplot(1, 2, 2);
bar(anion_effects_sorted, 'FaceColor', [0.3 0.7 0.3]);
hold on;
errorbar(1:n_anions, anion_effects_sorted, anion_std(anion_sort_idx) / sqrt(mean(anion_n)), ...
    'k.', 'LineWidth', 1.5);
set(gca, 'XTick', 1:n_anions, 'XTickLabel', unique_anions(anion_sort_idx), ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Random Effect on ln(a_w)', 'FontSize', 12);
title('Anion Effects (v_{anion})', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xl = xlim;
line([xl(1) xl(2)], [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

saveas(fig2, fullfile(OUTPUT_DIR, 'ion_random_effects.png'));
fprintf('Saved: ion_random_effects.png\n');
catch e
    fprintf('Error with figure 2: %s\n', e.message);
end

%% Figure 3: Fixed Effect f(m)
try
fig3 = figure('Position', [100, 100, 800, 600], 'visible', 'off');

% Plot individual salt curves (gray)
hold on;
for s = 1:min(20, n_salts)  % Plot first 20 salts for clarity
    salt_idx = false(n_obs, 1);
    for i = 1:n_obs
        if strcmp(data.Salt{i}, unique_salts{s})
            salt_idx(i) = true;
        end
    end
    m_salt = data.Molality(salt_idx);
    ln_aw_salt = data.ln_aw(salt_idx);
    
    [m_sorted, sort_idx] = sort(m_salt);
    plot(m_sorted, ln_aw_salt(sort_idx), '-', 'Color', [0.7 0.7 0.7], ...
        'LineWidth', 0.5);
end

% Plot fixed effect
m_plot = linspace(0, max(m_bins(valid_bins)), 200);
f_m_plot = f_m_interp(m_plot);
plot(m_plot, f_m_plot, 'b-', 'LineWidth', 3);

% Plot bins with error bars
errorbar(m_bins(valid_bins), f_m_bins(valid_bins), f_m_std(valid_bins), ...
    'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.5);

xlabel('Molality (mol/kg)', 'FontSize', 12);
ylabel('ln(a_w)', 'FontSize', 12);
title('Fixed Effect: Average Concentration Dependence f(m)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Individual Salts', '', '', '', '', '', '', '', '', '', '', ...
       '', '', '', '', '', '', '', '', '', 'Fixed Effect f(m)', ...
       'Binned Medians', 'Location', 'southwest', 'FontSize', 10);
grid on;

saveas(fig3, fullfile(OUTPUT_DIR, 'fixed_effect_f_m.png'));
fprintf('Saved: fixed_effect_f_m.png\n');
catch e
    fprintf('Error with figure 3: %s\n', e.message);
end

%% Figure 4: Model Fit (Observed vs Predicted)
try
fig4 = figure('Position', [100, 100, 800, 600], 'visible', 'off');

scatter(data.ln_aw, data.ln_aw_pred, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
hold on;
xl = xlim;
plot([xl(1) xl(2)], [xl(1) xl(2)], 'k--', 'LineWidth', 2);

xlabel('Observed ln(a_w)', 'FontSize', 12);
ylabel('Predicted ln(a_w)', 'FontSize', 12);
title(sprintf('Mixed-Effects Model Fit (R² = %.3f, RMSE = %.4f)', r2_total, rmse_total), ...
    'FontSize', 13, 'FontWeight', 'bold');
axis equal;
grid on;

% Add text box with statistics
text_str = sprintf('RMSE: %.4f\nR²: %.4f\nN: %d', rmse_total, r2_total, n_obs);
text(0.05, 0.95, text_str, 'Units', 'normalized', 'FontSize', 11, ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'VerticalAlignment', 'top');

saveas(fig4, fullfile(OUTPUT_DIR, 'model_fit_observed_vs_predicted.png'));
fprintf('Saved: model_fit_observed_vs_predicted.png\n');
catch e
    fprintf('Error with figure 4: %s\n', e.message);
end

%% Figure 5: Cross-Validation Results
try
fig5 = figure('Position', [100, 100, 1000, 400], 'visible', 'off');

subplot(1, 2, 1);
histogram(cv_rmse, 20, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'k');
hold on;
yl = ylim;
plot([mean(cv_rmse) mean(cv_rmse)], [yl(1) yl(2)], 'r--', 'LineWidth', 2);
plot([median(cv_rmse) median(cv_rmse)], [yl(1) yl(2)], 'g--', 'LineWidth', 2);
xlabel('RMSE', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Cross-Validation RMSE Distribution', 'FontSize', 13, 'FontWeight', 'bold');
legend('RMSE', 'Mean', 'Median', 'Location', 'best');
grid on;

subplot(1, 2, 2);
histogram(cv_r2, 20, 'FaceColor', [0.8 0.5 0.3], 'EdgeColor', 'k');
hold on;
yl = ylim;
plot([mean(cv_r2) mean(cv_r2)], [yl(1) yl(2)], 'r--', 'LineWidth', 2);
plot([median(cv_r2) median(cv_r2)], [yl(1) yl(2)], 'g--', 'LineWidth', 2);
xlabel('R²', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Cross-Validation R² Distribution', 'FontSize', 13, 'FontWeight', 'bold');
legend('R²', 'Mean', 'Median', 'Location', 'best');
grid on;

saveas(fig5, fullfile(OUTPUT_DIR, 'cross_validation_results.png'));
fprintf('Saved: cross_validation_results.png\n');
catch e
    fprintf('Error with figure 5: %s\n', e.message);
end

%% Figure 6: Example Predictions for Specific Salts
try
fig6 = figure('Position', [100, 100, 1200, 800], 'visible', 'off');

% Select 6 salts to showcase (best, worst, and some random)
showcase_salts = [best_idx, worst_idx, ...
                  randi(n_salts), randi(n_salts), randi(n_salts), randi(n_salts)];
showcase_salts = unique(showcase_salts);
showcase_salts = showcase_salts(1:min(6, length(showcase_salts)));

for i = 1:length(showcase_salts)
    subplot(2, 3, i);
    s_idx = showcase_salts(i);
    
    cv_data = cv_predictions{s_idx};
    
    scatter(cv_data.Molality, cv_data.ln_aw_true, 50, 'b', 'filled');
    hold on;
    plot(cv_data.Molality, cv_data.ln_aw_pred, 'r-', 'LineWidth', 2);
    
    xlabel('Molality (mol/kg)', 'FontSize', 10);
    ylabel('ln(a_w)', 'FontSize', 10);
    title(sprintf('%s\nRMSE=%.4f, R²=%.3f', cv_data.Salt, cv_rmse(s_idx), cv_r2(s_idx)), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
    legend('Observed', 'Predicted', 'Location', 'best', 'FontSize', 9);
    grid on;
end

% Note: sgtitle may not be available in all Octave versions
try
    sgtitle('Cross-Validation: Example Salt Predictions', 'FontSize', 14, 'FontWeight', 'bold');
catch
    % If sgtitle doesn't exist, use annotation instead
    annotation('textbox', [0 0.9 1 0.1], 'String', 'Cross-Validation: Example Salt Predictions', ...
               'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', ...
               'EdgeColor', 'none');
end

saveas(fig6, fullfile(OUTPUT_DIR, 'cv_example_predictions.png'));
fprintf('Saved: cv_example_predictions.png\n');
catch e
    fprintf('Error generating figure: %s\n', e.message);
end

end  % End of if graphics_available block

%% Save Results
fprintf('\n=== Saving Results ===\n');

% Ensure output directory exists
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('Created output directory: %s\n', OUTPUT_DIR);
end

% Save cation effects
fid = fopen(fullfile(OUTPUT_DIR, 'cation_effects.csv'), 'w');
if fid == -1
    error('Cannot open file for writing: %s', fullfile(OUTPUT_DIR, 'cation_effects.csv'));
end
fprintf(fid, 'Cation,Cation_Effect,Cation_StdDev,Cation_N_Obs\n');
for i = 1:n_cations
    fprintf(fid, '%s,%.6f,%.6f,%d\n', unique_cations{i}, cation_effects(i), cation_std(i), cation_n(i));
end
fclose(fid);
fprintf('Saved: cation_effects.csv\n');

% Save anion effects
fid = fopen(fullfile(OUTPUT_DIR, 'anion_effects.csv'), 'w');
if fid == -1
    error('Cannot open file for writing: %s', fullfile(OUTPUT_DIR, 'anion_effects.csv'));
end
fprintf(fid, 'Anion,Anion_Effect,Anion_StdDev,Anion_N_Obs\n');
for i = 1:n_anions
    fprintf(fid, '%s,%.6f,%.6f,%d\n', unique_anions{i}, anion_effects(i), anion_std(i), anion_n(i));
end
fclose(fid);
fprintf('Saved: anion_effects.csv\n');

% Save CV results
fid = fopen(fullfile(OUTPUT_DIR, 'cross_validation_results.csv'), 'w');
if fid == -1
    error('Cannot open file for writing: %s', fullfile(OUTPUT_DIR, 'cross_validation_results.csv'));
end
fprintf(fid, 'Salt,RMSE,R2\n');
for i = 1:n_salts
    fprintf(fid, '%s,%.6f,%.6f\n', unique_salts{i}, cv_rmse(i), cv_r2(i));
end
fclose(fid);
fprintf('Saved: cross_validation_results.csv\n');

% Save summary statistics
fid = fopen(fullfile(OUTPUT_DIR, 'mixed_effects_summary.txt'), 'w');
if fid == -1
    error('Cannot open file for writing: %s', fullfile(OUTPUT_DIR, 'mixed_effects_summary.txt'));
end
fprintf(fid, '=== Mixed-Effects Model Summary ===\n\n');
fprintf(fid, 'Model: ln(aw_ij) = f(m_ij) + u_cation(i) + v_anion(i) + e_ij\n\n');
fprintf(fid, 'Data:\n');
fprintf(fid, '  Total observations: %d\n', n_obs);
fprintf(fid, '  Unique salts: %d\n', n_salts);
fprintf(fid, '  Unique cations: %d\n', n_cations);
fprintf(fid, '  Unique anions: %d\n', n_anions);
fprintf(fid, '\nModel Performance:\n');
fprintf(fid, '  RMSE: %.4f\n', rmse_total);
fprintf(fid, '  R²: %.4f\n', r2_total);
fprintf(fid, '\nVariance Decomposition:\n');
fprintf(fid, '  Total variance: %.4f\n', var_total);
fprintf(fid, '  Fixed effect (concentration): %.4f (%.1f%%)\n', var_fixed, pct_fixed);
fprintf(fid, '  Random effect (cation): %.4f (%.1f%%)\n', var_cation, pct_cation);
fprintf(fid, '  Random effect (anion): %.4f (%.1f%%)\n', var_anion, pct_anion);
fprintf(fid, '  Residual: %.4f (%.1f%%)\n', var_residual, pct_residual);
fprintf(fid, '\nCation:Anion variance ratio: %.2f\n', cation_anion_ratio);
fprintf(fid, '\nCross-Validation Results:\n');
fprintf(fid, '  Mean RMSE: %.4f ± %.4f\n', mean(cv_rmse), std(cv_rmse));
fprintf(fid, '  Mean R²: %.4f ± %.4f\n', mean(cv_r2), std(cv_r2));
fprintf(fid, '  Median RMSE: %.4f\n', median(cv_rmse));
fprintf(fid, '  Median R²: %.4f\n', median(cv_r2));
fclose(fid);

fprintf('Saved: mixed_effects_summary.txt\n');

fprintf('\n=== Analysis Complete ===\n');
fprintf('All results saved to: %s\n', OUTPUT_DIR);
