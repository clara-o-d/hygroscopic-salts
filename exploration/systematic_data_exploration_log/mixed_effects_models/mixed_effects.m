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

% Check if we have enough points for interpolation
if length(m_bins_valid) < 2
    error('Not enough valid bins for interpolation. Need at least 2 valid bins.');
end

% Use interpolation with robust fallback
if length(m_bins_valid) >= 4
    % Use pchip for smooth interpolation (requires at least 4 points)
    f_m_values = interp1(m_bins_valid, f_m_bins_valid, data.Molality, 'pchip');
else
    % Use linear interpolation if we have 2-3 points
    f_m_values = interp1(m_bins_valid, f_m_bins_valid, data.Molality, 'linear');
end

% For any NaN values (extrapolation or interpolation failure), use nearest neighbor
nan_idx = isnan(f_m_values);
if any(nan_idx)
    fprintf('Warning: %d points outside interpolation range, using nearest neighbor\n', sum(nan_idx));
    % Use nearest neighbor with extrapolation
    f_m_values(nan_idx) = interp1(m_bins_valid, f_m_bins_valid, data.Molality(nan_idx), 'nearest', 'extrap');
end

% Final check: if any values are still NaN, use a constant fallback
if any(isnan(f_m_values))
    fprintf('Warning: %d points still have NaN after interpolation, using median value\n', sum(isnan(f_m_values)));
    f_m_values(isnan(f_m_values)) = median(f_m_bins_valid);
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

% Calculate model fit (filter out NaN values)
valid_pred_idx = ~isnan(data.ln_aw_pred) & ~isnan(data.ln_aw);
if ~all(valid_pred_idx)
    fprintf('Warning: %d predictions have NaN values, excluding from RMSE/R² calculation\n', sum(~valid_pred_idx));
end

for i = 1:n_obs
    data.residual_final(i) = data.ln_aw(i) - data.ln_aw_pred(i);
end

% Calculate RMSE only on valid predictions
if any(valid_pred_idx)
    rmse_total = sqrt(mean(data.residual_final(valid_pred_idx).^2));
else
    rmse_total = NaN;
    fprintf('Error: No valid predictions for RMSE calculation\n');
end

% R² calculation with numerical stability check (only on valid predictions)
if any(valid_pred_idx)
    ss_res_total = sum(data.residual_final(valid_pred_idx).^2);
    ln_aw_valid = data.ln_aw(valid_pred_idx);
    ss_tot_total = sum((ln_aw_valid - mean(ln_aw_valid)).^2);
    if ss_tot_total < 1e-10
        r2_total = NaN;
    else
        r2_total = 1 - ss_res_total / ss_tot_total;
    end
else
    r2_total = NaN;
end

fprintf('\nModel Performance:\n');
fprintf('  RMSE: %.4f\n', rmse_total);
fprintf('  R²: %.4f\n', r2_total);

%% Step 4: Variance Decomposition
fprintf('\n=== Step 4: Variance Decomposition ===\n');

% Total variance in ln(aw)
var_total = var(data.ln_aw);

% For proper variance decomposition, we need to calculate the variance
% contribution of each component. Since the model is:
% ln(aw) = f(m) + u_cation + v_anion + ε
% and these components are (approximately) uncorrelated, we can decompose
% the variance of the predicted values.

% Get contributions for each component (only on valid predictions)
if any(valid_pred_idx)
    % Fixed effect contribution
    f_m_valid = data.f_m(valid_pred_idx);
    var_fixed = var(f_m_valid);
    
    % Cation effects contribution
    cation_var_contributions = zeros(sum(valid_pred_idx), 1);
    valid_idx_list = find(valid_pred_idx);
    for j = 1:length(valid_idx_list)
        i = valid_idx_list(j);
        for c = 1:n_cations
            if strcmp(data.Cation{i}, unique_cations{c})
                cation_var_contributions(j) = cation_effects(c);
                break;
            end
        end
    end
    var_cation = var(cation_var_contributions);
    
    % Anion effects contribution
    anion_var_contributions = zeros(sum(valid_pred_idx), 1);
    for j = 1:length(valid_idx_list)
        i = valid_idx_list(j);
        for a = 1:n_anions
            if strcmp(data.Anion{i}, unique_anions{a})
                anion_var_contributions(j) = anion_effects(a);
                break;
            end
        end
    end
    var_anion = var(anion_var_contributions);
    
    % Residual variance
    var_residual = var(data.residual_final(valid_pred_idx));
    
    % Variance of predicted values (explained variance)
    var_predicted = var(data.ln_aw_pred(valid_pred_idx));
    var_explained = var_predicted;
    
    % The variance components are correlated, so they don't add up directly.
    % We normalize them so they sum to the explained variance.
    % The correct decomposition: var_total = var_explained + var_residual
    % where var_explained = var(f(m) + u + v)
    
    % Normalize component variances to sum to explained variance
    var_sum_components = var_fixed + var_cation + var_anion;
    if var_sum_components > 0 && var_explained > 0
        % Scale factors to make components sum to explained variance
        scale_factor = var_explained / var_sum_components;
        var_fixed = var_fixed * scale_factor;
        var_cation = var_cation * scale_factor;
        var_anion = var_anion * scale_factor;
    end
    
    % Now compute percentages - normalize all components to sum to 100% of total variance
    % Total variance = explained + residual (approximately)
    var_total_check = var_explained + var_residual;
    
    % Sum of all variance components (for normalization)
    var_sum_all = var_fixed + var_cation + var_anion + var_residual;
    
    % Normalize so percentages sum to 100%
    if var_sum_all > 0
        % Scale all components proportionally to sum to total variance
        normalization_factor = var_total / var_sum_all;
        var_fixed_normalized = var_fixed * normalization_factor;
        var_cation_normalized = var_cation * normalization_factor;
        var_anion_normalized = var_anion * normalization_factor;
        var_residual_normalized = var_residual * normalization_factor;
        
        % Update for reporting
        var_fixed = var_fixed_normalized;
        var_cation = var_cation_normalized;
        var_anion = var_anion_normalized;
        var_residual = var_residual_normalized;
        
        % Debug check
        fprintf('After normalization: sum=%.4f, total=%.4f\n', ...
            var_fixed + var_cation + var_anion + var_residual, var_total);
    end
    
    % Compute percentages relative to total variance (should sum to 100%)
    pct_fixed = 100 * var_fixed / var_total;
    pct_cation = 100 * var_cation / var_total;
    pct_anion = 100 * var_anion / var_total;
    pct_residual = 100 * var_residual / var_total;
    
    fprintf('Variance of predicted (explained): %.4f\n', var_explained);
    fprintf('Residual variance: %.4f\n', var_residual);
    fprintf('Sum of all components: %.4f (should equal total: %.4f)\n', ...
        var_fixed + var_cation + var_anion + var_residual, var_total);
else
    var_fixed = NaN;
    var_cation = NaN;
    var_anion = NaN;
    var_residual = NaN;
    pct_fixed = NaN;
    pct_cation = NaN;
    pct_anion = NaN;
    pct_residual = NaN;
end

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
    
    % Check if we have enough valid bins
    if length(m_bins_train_valid) < 2
        % Not enough bins, skip this CV fold
        cv_rmse(s) = NaN;
        cv_r2(s) = NaN;
        cv_predictions{s} = struct('Salt', test_salt, ...
                                   'Cation', test_cation, ...
                                   'Anion', test_anion, ...
                                   'Molality', m_test, ...
                                   'ln_aw_true', ln_aw_test, ...
                                   'ln_aw_pred', NaN(size(m_test)));
        continue;
    end
    
    % Re-estimate random effects on training data
    f_m_train_values = interp1_with_nearest(m_bins_train_valid, f_m_bins_train_valid, m_train);
    
    % Check for NaN values in f_m_train_values
    if any(isnan(f_m_train_values))
        fprintf('Warning: CV fold %d has NaN in f_m_train_values\n', s);
        f_m_train_values(isnan(f_m_train_values)) = median(f_m_bins_train_valid);
    end
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
    
    % Check for NaN in f_m_test
    if any(isnan(f_m_test))
        f_m_test(isnan(f_m_test)) = median(f_m_bins_train_valid);
    end
    
    cation_idx = find(strcmp(unique_cations, test_cation));
    anion_idx = find(strcmp(unique_anions, test_anion));
    
    % Check if cation/anion were in training set
    if isempty(cation_idx) || isempty(anion_idx)
        % Ion not in training set, use zero effect (already centered)
        cation_effect_test = 0;
        anion_effect_test = 0;
    else
        cation_effect_test = cation_effects_train(cation_idx);
        anion_effect_test = anion_effects_train(anion_idx);
    end
    
    ln_aw_pred_test = f_m_test + cation_effect_test + anion_effect_test;
    
    % Store results
    cv_predictions{s} = struct('Salt', test_salt, ...
                               'Cation', test_cation, ...
                               'Anion', test_anion, ...
                               'Molality', m_test, ...
                               'ln_aw_true', ln_aw_test, ...
                               'ln_aw_pred', ln_aw_pred_test);
    
    % Calculate metrics (check for NaN)
    valid_test_idx = ~isnan(ln_aw_test) & ~isnan(ln_aw_pred_test);
    if any(valid_test_idx)
        cv_rmse(s) = sqrt(mean((ln_aw_test(valid_test_idx) - ln_aw_pred_test(valid_test_idx)).^2));
    else
        cv_rmse(s) = NaN;
    end
    
    % R² calculation with numerical stability check (only on valid predictions)
    if any(valid_test_idx)
        ss_res = sum((ln_aw_test(valid_test_idx) - ln_aw_pred_test(valid_test_idx)).^2);
        ln_aw_test_valid = ln_aw_test(valid_test_idx);
        ss_tot = sum((ln_aw_test_valid - mean(ln_aw_test_valid)).^2);
        
        if ss_tot < 1e-10  % Very small variance in test data
            % If test data has almost no variance, R² is undefined
            cv_r2(s) = NaN;
        else
            cv_r2(s) = 1 - ss_res / ss_tot;
            % Cap R² at reasonable bounds to avoid extreme values
            if cv_r2(s) < -10
                cv_r2(s) = -10;  % Cap at -10 for very poor predictions
            end
        end
    else
        cv_r2(s) = NaN;
    end
end

fprintf('\n');

fprintf('\nCross-Validation Results:\n');
fprintf('  Mean RMSE: %.4f ± %.4f\n', mean(cv_rmse), std(cv_rmse));
% Handle NaN values in R²
valid_r2 = cv_r2(~isnan(cv_r2));
if ~isempty(valid_r2)
    fprintf('  Mean R²: %.4f ± %.4f\n', mean(valid_r2), std(valid_r2));
    fprintf('  Median R²: %.4f\n', median(valid_r2));
    fprintf('  Valid R² values: %d/%d\n', length(valid_r2), length(cv_r2));
else
    fprintf('  Mean R²: NaN (all values invalid)\n');
    fprintf('  Median R²: NaN\n');
end
fprintf('  Median RMSE: %.4f\n', median(cv_rmse));

% Best and worst predictions
[~, best_idx] = min(cv_rmse);
[~, worst_idx] = max(cv_rmse);

fprintf('\nBest prediction: %s (RMSE=%.4f', unique_salts{best_idx}, cv_rmse(best_idx));
if ~isnan(cv_r2(best_idx))
    fprintf(', R²=%.4f', cv_r2(best_idx));
else
    fprintf(', R²=NaN');
end
fprintf(')\n');

fprintf('Worst prediction: %s (RMSE=%.4f', unique_salts{worst_idx}, cv_rmse(worst_idx));
if ~isnan(cv_r2(worst_idx))
    fprintf(', R²=%.4f', cv_r2(worst_idx));
else
    fprintf(', R²=NaN');
end
fprintf(')\n');

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

% For the pie chart, show the variance components that sum to total variance
% Include residual to show complete decomposition
variance_components = [var_fixed, var_cation, var_anion, var_residual];
labels_full = {'Fixed Effect (Concentration)', 'Cation Identity', 'Anion Identity', 'Residual'};
colors = [0.2 0.4 0.8; 0.8 0.3 0.3; 0.3 0.7 0.3; 0.7 0.7 0.7];

% DON'T filter out any components - include ALL to get full 360 degree pie
if all(~isnan(variance_components))
    % Calculate percentages
    pcts = 100 * variance_components / sum(variance_components);
    
    % Create labels with percentages
    labels_with_pct = cell(size(labels_full));
    for i = 1:length(labels_full)
        labels_with_pct{i} = sprintf('%s (%.1f%%)', labels_full{i}, pcts(i));
    end
    
    % Create pie chart with ALL components (this ensures 360 degree circle)
    pie(variance_components, labels_with_pct);
    colormap(colors);
    
    % Get current axes and adjust to fill figure
    ax = gca;
    set(ax, 'Position', [0.05 0.05 0.9 0.9]);  % Fill most of the figure
    axis equal;
    
    title('Variance Decomposition in ln(a_w)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add note about total variance
    annotation_text = sprintf('Total Variance: %.4f', var_total);
    text(0.02, 0.02, annotation_text, 'Units', 'normalized', 'FontSize', 10, ...
         'BackgroundColor', 'w', 'EdgeColor', 'k');
else
    text(0.5, 0.5, 'No valid variance components', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Variance Decomposition in ln(a_w)', 'FontSize', 14, 'FontWeight', 'bold');
end

% Add percentage labels (Octave-compatible)
try
    text_objs = findobj(fig1, 'Type', 'text');
    for i = 1:length(text_objs)
        try
            curr_text = get(text_objs(i), 'String');
            if ~isempty(strfind(curr_text, '%'))
                set(text_objs(i), 'FontSize', 12, 'FontWeight', 'bold');
            end
        catch
            % Skip if property access fails
        end
    end
catch
    % Skip if text object manipulation fails
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
% Error bars (with fallback for Octave compatibility)
try
    error_y = cation_std(cation_sort_idx) / sqrt(mean(cation_n));
    errorbar(1:n_cations, cation_effects_sorted, error_y, 'k.', 'LineWidth', 1.5);
catch
    % Skip error bars if not supported
end
set(gca, 'XTick', 1:n_cations, 'XTickLabel', unique_cations(cation_sort_idx), ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Random Effect on ln(a_w)', 'FontSize', 12);
title('Cation Effects (u_{cation})', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xl = xlim;
% Use plot instead of line to avoid variable name conflicts
plot([xl(1) xl(2)], [0 0], 'k--', 'LineWidth', 1.5);

subplot(1, 2, 2);
bar(anion_effects_sorted, 'FaceColor', [0.3 0.7 0.3]);
hold on;
% Error bars (with fallback for Octave compatibility)
try
    error_y_anion = anion_std(anion_sort_idx) / sqrt(mean(anion_n));
    errorbar(1:n_anions, anion_effects_sorted, error_y_anion, 'k.', 'LineWidth', 1.5);
catch
    % Skip error bars if not supported
end
set(gca, 'XTick', 1:n_anions, 'XTickLabel', unique_anions(anion_sort_idx), ...
    'XTickLabelRotation', 45, 'FontSize', 10);
ylabel('Random Effect on ln(a_w)', 'FontSize', 12);
title('Anion Effects (v_{anion})', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
xl = xlim;
% Use plot instead of line to avoid variable name conflicts
plot([xl(1) xl(2)], [0 0], 'k--', 'LineWidth', 1.5);

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

% Plot bins with error bars (Octave-compatible)
try
    errorbar(m_bins(valid_bins), f_m_bins(valid_bins), f_m_std(valid_bins), ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.5);
catch
    % Fallback: plot without error bars if errorbar fails
    plot(m_bins(valid_bins), f_m_bins(valid_bins), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.5);
end

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

% Use scatter with Octave-compatible syntax (MarkerFaceAlpha not supported)
try
    scatter(data.ln_aw, data.ln_aw_pred, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
catch
    scatter(data.ln_aw, data.ln_aw_pred, 20, 'b', 'filled');
end
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
% Use hist for Octave compatibility (histogram is MATLAB-only)
[n_rmse, edges_rmse] = hist(cv_rmse, 20);
bar(edges_rmse, n_rmse, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'k');
hold on;
yl = ylim;
plot([mean(cv_rmse) mean(cv_rmse)], [yl(1) yl(2)], 'r--', 'LineWidth', 2);
plot([median(cv_rmse) median(cv_rmse)], [yl(1) yl(2)], 'g--', 'LineWidth', 2);
xlabel('RMSE', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Cross-Validation RMSE Distribution', 'FontSize', 13, 'FontWeight', 'bold');
legend('RMSE', 'Mean', 'Median', 'Location', 'northeast');
grid on;

subplot(1, 2, 2);
% Filter out NaN values for R² histogram
valid_r2_plot = cv_r2(~isnan(cv_r2));
if ~isempty(valid_r2_plot)
    [n_r2, edges_r2] = hist(valid_r2_plot, 20);
    bar(edges_r2, n_r2, 'FaceColor', [0.8 0.5 0.3], 'EdgeColor', 'k');
    hold on;
    yl = ylim;
    plot([mean(valid_r2_plot) mean(valid_r2_plot)], [yl(1) yl(2)], 'r--', 'LineWidth', 2);
    plot([median(valid_r2_plot) median(valid_r2_plot)], [yl(1) yl(2)], 'g--', 'LineWidth', 2);
    xlabel('R²', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title('Cross-Validation R² Distribution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('R²', 'Mean', 'Median', 'Location', 'northeast');
else
    text(0.5, 0.5, 'No valid R² values', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('Cross-Validation R² Distribution', 'FontSize', 13, 'FontWeight', 'bold');
end
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
% Handle NaN values in R²
valid_r2 = cv_r2(~isnan(cv_r2));
if ~isempty(valid_r2)
    fprintf(fid, '  Mean R²: %.4f ± %.4f\n', mean(valid_r2), std(valid_r2));
    fprintf(fid, '  Median R²: %.4f\n', median(valid_r2));
    fprintf(fid, '  Valid R² values: %d/%d\n', length(valid_r2), length(cv_r2));
else
    fprintf(fid, '  Mean R²: NaN (all values invalid)\n');
    fprintf(fid, '  Median R²: NaN\n');
end
fprintf(fid, '  Median RMSE: %.4f\n', median(cv_rmse));
fclose(fid);

fprintf('Saved: mixed_effects_summary.txt\n');

fprintf('\n=== Analysis Complete ===\n');
fprintf('All results saved to: %s\n', OUTPUT_DIR);
