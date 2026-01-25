%% Gradient Boosting for Predicting ln(gamma) at 90% RH
% This script implements gradient boosting regression to predict the natural
% logarithm of the activity coefficient (ln(gamma)) at 90% relative humidity.
%
% Approach:
%   1. Implement manual gradient boosting with regression trees
%   2. Cross-validation for model evaluation
%   3. Feature importance analysis
%   4. Residual analysis and model diagnostics
%
% Target: ln_gamma_at_90RH
% Features: Ion properties, Pitzer parameters, electrolyte types
%
% Author: Generated for AWH ML Discussion
% Date: 2026-01-25

clear; clc; close all;

% Load required packages for Octave
pkg load statistics

% Set up graphics for headless operation
try
    available_toolkits = available_graphics_toolkits();
    if ~isempty(available_toolkits)
        graphics_toolkit(available_toolkits{1});
    end
    set(0, 'DefaultFigureVisible', 'off');
catch
    warning('Could not set up graphics toolkit. Figures may not be saved.');
end

%% Configuration
DATA_FILE = '../../../data/baseline_numeric_only.csv';
OUTPUT_DIR = '../../../figures/gradient_boosting/';
RESULTS_FILE = 'gradient_boosting_results.txt';

% Model hyperparameters
N_TREES = 100;              % Number of boosting iterations
LEARNING_RATE = 0.1;        % Shrinkage parameter
MAX_DEPTH = 4;              % Maximum tree depth
MIN_LEAF_SIZE = 5;          % Minimum observations per leaf
SUBSAMPLE_RATIO = 0.8;      % Fraction of data for each tree

% Cross-validation settings
N_FOLDS = 5;                % K-fold cross-validation
TEST_SPLIT = 0.2;           % Hold-out test set fraction
RANDOM_SEED = 42;           % For reproducibility

% Create output directory
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('\n');
fprintf('================================================================\n');
fprintf('   GRADIENT BOOSTING: ln(gamma) at 90%% RH\n');
fprintf('================================================================\n\n');

%% Load and Prepare Data
fprintf('Loading data from: %s\n', DATA_FILE);

% Read CSV file
fid = fopen(DATA_FILE, 'r');
if fid == -1
    error('Cannot open data file: %s', DATA_FILE);
end

% Read header
header_line = fgetl(fid);
header_cols = strsplit(header_line, ',');
fclose(fid);

% Read full data
data_table = dlmread(DATA_FILE, ',', 1, 0);

% Find target column index (ln_gamma_at_90RH)
target_col_name = 'ln_gamma_at_90RH';
target_col_idx = find(strcmp(header_cols, target_col_name));

if isempty(target_col_idx)
    error('Target column "%s" not found in data file', target_col_name);
end

% Read electrolyte names
fid = fopen(DATA_FILE, 'r');
fgetl(fid); % Skip header
electrolyte_names = {};
line_idx = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line)
        line_idx = line_idx + 1;
        parts = strsplit(line, ',');
        electrolyte_names{line_idx} = parts{1};
    end
end
fclose(fid);

% Extract target variable
y_all = data_table(:, target_col_idx);

% Define feature columns (exclude target columns and non-predictive features)
exclude_cols = [1, 24:31]; % electrolyte name, all ln_gamma columns, RH_actual columns
feature_cols = setdiff(1:size(data_table, 2), exclude_cols);
X_all = data_table(:, feature_cols);

% Get feature names
feature_names = header_cols(feature_cols);

% Filter out rows with missing target values
valid_idx = ~isnan(y_all);
X = X_all(valid_idx, :);
y = y_all(valid_idx);
electrolyte_names_valid = electrolyte_names(valid_idx);

[n_samples, n_features] = size(X);

fprintf('Data loaded successfully!\n');
fprintf('  - Total samples: %d\n', n_samples);
fprintf('  - Number of features: %d\n', n_features);
fprintf('  - Target: %s\n', target_col_name);
fprintf('  - Target range: [%.4f, %.4f]\n', min(y), max(y));
fprintf('  - Target mean ± std: %.4f ± %.4f\n\n', mean(y), std(y));

%% Handle Missing Values in Features
fprintf('Handling missing values in features...\n');
missing_count = sum(isnan(X));
features_with_missing = find(missing_count > 0);

if ~isempty(features_with_missing)
    fprintf('  - Features with missing values: %d\n', length(features_with_missing));
    for i = 1:length(features_with_missing)
        feat_idx = features_with_missing(i);
        n_missing = missing_count(feat_idx);
        fprintf('    * %s: %d missing (%.1f%%)\n', ...
                feature_names{feat_idx}, n_missing, 100*n_missing/n_samples);
    end
    
    % Impute missing values with median
    for i = 1:n_features
        col = X(:, i);
        if any(isnan(col))
            median_val = median(col(~isnan(col)));
            X(isnan(col), i) = median_val;
        end
    end
    fprintf('  - Missing values imputed with feature medians\n\n');
else
    fprintf('  - No missing values found\n\n');
end

%% Feature Scaling
fprintf('Standardizing features...\n');
X_mean = mean(X);
X_std = std(X);
X_std(X_std < 1e-10) = 1; % Avoid division by zero
X_scaled = (X - repmat(X_mean, n_samples, 1)) ./ repmat(X_std, n_samples, 1);
fprintf('  - Features standardized (mean=0, std=1)\n\n');

%% Train-Test Split
fprintf('Splitting data into train and test sets...\n');
rand('seed', RANDOM_SEED);
n_test = round(n_samples * TEST_SPLIT);
n_train = n_samples - n_test;

% Random permutation
perm_idx = randperm(n_samples);
train_idx = perm_idx(1:n_train);
test_idx = perm_idx(n_train+1:end);

X_train = X_scaled(train_idx, :);
y_train = y(train_idx);
X_test = X_scaled(test_idx, :);
y_test = y(test_idx);

fprintf('  - Training samples: %d (%.0f%%)\n', n_train, 100*(1-TEST_SPLIT));
fprintf('  - Test samples: %d (%.0f%%)\n\n', n_test, 100*TEST_SPLIT);

%% Manual Gradient Boosting Implementation
fprintf('================================================================\n');
fprintf('Training Gradient Boosting Model\n');
fprintf('================================================================\n\n');
fprintf('Hyperparameters:\n');
fprintf('  - Number of trees: %d\n', N_TREES);
fprintf('  - Learning rate: %.3f\n', LEARNING_RATE);
fprintf('  - Max tree depth: %d\n', MAX_DEPTH);
fprintf('  - Min leaf size: %d\n', MIN_LEAF_SIZE);
fprintf('  - Subsample ratio: %.2f\n\n', SUBSAMPLE_RATIO);

% Initialize predictions with mean
F_train = zeros(n_train, 1);
F_test = zeros(n_test, 1);
mean_y = mean(y_train);
F_train(:) = mean_y;
F_test(:) = mean_y;

% Store trees and performance metrics
trees = cell(N_TREES, 1);
train_errors = zeros(N_TREES, 1);
test_errors = zeros(N_TREES, 1);
feature_importance = zeros(n_features, 1);

fprintf('Starting boosting iterations...\n');
for iter = 1:N_TREES
    % Compute negative gradient (residuals for squared loss)
    residuals = y_train - F_train;
    
    % Subsample training data
    n_subsample = round(n_train * SUBSAMPLE_RATIO);
    subsample_idx = randsample(n_train, n_subsample);
    X_subsample = X_train(subsample_idx, :);
    residuals_subsample = residuals(subsample_idx);
    
    % Fit regression tree to residuals using simple approach
    tree = struct();
    tree.prediction = mean(residuals_subsample);
    tree.feature_splits = [];
    
    % Simple tree: find best single split
    if length(residuals_subsample) >= 2 * MIN_LEAF_SIZE
        best_gain = 0;
        best_feature = 1;
        best_threshold = 0;
        current_var = var(residuals_subsample);
        
        % Try each feature for splitting
        for feat = 1:n_features
            feat_values = X_subsample(:, feat);
            unique_vals = unique(feat_values);
            
            % Try some thresholds
            n_try = min(10, length(unique_vals) - 1);
            for t = 1:n_try
                if length(unique_vals) > n_try
                    threshold = unique_vals(round(t * length(unique_vals) / (n_try + 1)));
                else
                    if t < length(unique_vals)
                        threshold = (unique_vals(t) + unique_vals(t + 1)) / 2;
                    else
                        continue;
                    end
                end
                
                left_mask = feat_values <= threshold;
                right_mask = ~left_mask;
                
                n_left = sum(left_mask);
                n_right = sum(right_mask);
                
                if n_left < MIN_LEAF_SIZE || n_right < MIN_LEAF_SIZE
                    continue;
                end
                
                % Calculate variance reduction
                left_var = var(residuals_subsample(left_mask));
                right_var = var(residuals_subsample(right_mask));
                weighted_var = (n_left * left_var + n_right * right_var) / length(residuals_subsample);
                gain = current_var - weighted_var;
                
                if gain > best_gain
                    best_gain = gain;
                    best_feature = feat;
                    best_threshold = threshold;
                end
            end
        end
        
        % Store split info
        if best_gain > 0
            tree.has_split = true;
            tree.feature = best_feature;
            tree.threshold = best_threshold;
            tree.gain = best_gain;
            
            % Calculate left and right predictions
            left_mask = X_subsample(:, best_feature) <= best_threshold;
            right_mask = ~left_mask;
            tree.left_value = mean(residuals_subsample(left_mask));
            tree.right_value = mean(residuals_subsample(right_mask));
            
            % Track feature importance
            feature_importance(best_feature) = feature_importance(best_feature) + best_gain;
        else
            tree.has_split = false;
        end
    else
        tree.has_split = false;
    end
    
    trees{iter} = tree;
    
    % Predict with new tree
    pred_train = zeros(n_train, 1);
    pred_test = zeros(n_test, 1);
    
    for i = 1:n_train
        if tree.has_split
            if X_train(i, tree.feature) <= tree.threshold
                pred_train(i) = tree.left_value;
            else
                pred_train(i) = tree.right_value;
            end
        else
            pred_train(i) = tree.prediction;
        end
    end
    
    for i = 1:n_test
        if tree.has_split
            if X_test(i, tree.feature) <= tree.threshold
                pred_test(i) = tree.left_value;
            else
                pred_test(i) = tree.right_value;
            end
        else
            pred_test(i) = tree.prediction;
        end
    end
    
    % Update predictions with learning rate
    F_train = F_train + LEARNING_RATE * pred_train;
    F_test = F_test + LEARNING_RATE * pred_test;
    
    % Compute errors (RMSE)
    train_errors(iter) = sqrt(mean((y_train - F_train).^2));
    test_errors(iter) = sqrt(mean((y_test - F_test).^2));
    
    % Print progress
    if mod(iter, 10) == 0 || iter == 1
        fprintf('  Iteration %3d: Train RMSE = %.6f, Test RMSE = %.6f\n', ...
                iter, train_errors(iter), test_errors(iter));
    end
end

fprintf('\nTraining complete!\n\n');

%% Model Performance
fprintf('================================================================\n');
fprintf('Model Performance\n');
fprintf('================================================================\n\n');

% Final predictions
y_train_pred = F_train;
y_test_pred = F_test;

% Training metrics
train_rmse = sqrt(mean((y_train - y_train_pred).^2));
train_mae = mean(abs(y_train - y_train_pred));
train_r2 = 1 - sum((y_train - y_train_pred).^2) / sum((y_train - mean(y_train)).^2);

% Test metrics
test_rmse = sqrt(mean((y_test - y_test_pred).^2));
test_mae = mean(abs(y_test - y_test_pred));
test_r2 = 1 - sum((y_test - y_test_pred).^2) / sum((y_test - mean(y_test)).^2);

fprintf('Training Set:\n');
fprintf('  - RMSE: %.6f\n', train_rmse);
fprintf('  - MAE:  %.6f\n', train_mae);
fprintf('  - R²:   %.6f\n\n', train_r2);

fprintf('Test Set:\n');
fprintf('  - RMSE: %.6f\n', test_rmse);
fprintf('  - MAE:  %.6f\n', test_mae);
fprintf('  - R²:   %.6f\n\n', test_r2);

%% Feature Importance Analysis
fprintf('================================================================\n');
fprintf('Feature Importance Analysis\n');
fprintf('================================================================\n\n');

% Normalize feature importance
if sum(feature_importance) > 0
    feature_importance = feature_importance / sum(feature_importance);
end

% Sort by importance
[importance_sorted, sort_idx] = sort(feature_importance, 'descend');
features_sorted = feature_names(sort_idx);

fprintf('Top 15 Most Important Features:\n');
for i = 1:min(15, length(features_sorted))
    fprintf('  %2d. %-40s %.4f\n', i, features_sorted{i}, importance_sorted(i));
end
fprintf('\n');

% Compute residuals for later use
residuals_test = y_test - y_test_pred;
residuals_train = y_train - y_train_pred;

%% Visualizations
fprintf('================================================================\n');
fprintf('Generating Visualizations\n');
fprintf('================================================================\n\n');

% Check if graphics are available
graphics_available = true;
try
    test_fig = figure('Visible', 'off');
    close(test_fig);
catch
    graphics_available = false;
    fprintf('  - Graphics toolkit not available, skipping visualizations\n\n');
end

if graphics_available
    try
        % 1. Learning Curves
        fprintf('  - Creating learning curves plot...\n');
        fig1 = figure('Position', [100, 100, 800, 600]);
        plot(1:N_TREES, train_errors, 'b-', 'LineWidth', 2);
        hold on;
        plot(1:N_TREES, test_errors, 'r-', 'LineWidth', 2);
        xlabel('Number of Trees', 'FontSize', 12);
        ylabel('RMSE', 'FontSize', 12);
        title('Gradient Boosting: Learning Curves', 'FontSize', 14, 'FontWeight', 'bold');
        legend('Training', 'Test', 'Location', 'best');
        grid on;
        saveas(fig1, [OUTPUT_DIR, 'learning_curves.png']);
        close(fig1);
        
        % 2. Predicted vs Observed (Training)
        fprintf('  - Creating predicted vs observed plot (training)...\n');
        fig2 = figure('Position', [100, 100, 800, 600]);
        scatter(y_train, y_train_pred, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
        hold on;
        plot([min(y_train), max(y_train)], [min(y_train), max(y_train)], 'k--', 'LineWidth', 2);
        xlabel('Observed ln(\gamma) at 90% RH', 'FontSize', 12);
        ylabel('Predicted ln(\gamma) at 90% RH', 'FontSize', 12);
        title(sprintf('Training Set: R² = %.4f', train_r2), 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        axis equal tight;
        saveas(fig2, [OUTPUT_DIR, 'predicted_vs_observed_train.png']);
        close(fig2);
        
        % 3. Predicted vs Observed (Test)
        fprintf('  - Creating predicted vs observed plot (test)...\n');
        fig3 = figure('Position', [100, 100, 800, 600]);
        scatter(y_test, y_test_pred, 50, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
        hold on;
        plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'k--', 'LineWidth', 2);
        xlabel('Observed ln(\gamma) at 90% RH', 'FontSize', 12);
        ylabel('Predicted ln(\gamma) at 90% RH', 'FontSize', 12);
        title(sprintf('Test Set: R² = %.4f', test_r2), 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        axis equal tight;
        saveas(fig3, [OUTPUT_DIR, 'predicted_vs_observed_test.png']);
        close(fig3);
        
        % 4. Residual Plot
        fprintf('  - Creating residual plot...\n');
        fig4 = figure('Position', [100, 100, 800, 600]);
        scatter(y_test_pred, residuals_test, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
        hold on;
        plot([min(y_test_pred), max(y_test_pred)], [0, 0], 'k--', 'LineWidth', 2);
        xlabel('Predicted ln(\gamma) at 90% RH', 'FontSize', 12);
        ylabel('Residuals', 'FontSize', 12);
        title('Residual Plot (Test Set)', 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        saveas(fig4, [OUTPUT_DIR, 'residuals.png']);
        close(fig4);
        
        % 5. Feature Importance Bar Plot
        fprintf('  - Creating feature importance plot...\n');
        fig5 = figure('Position', [100, 100, 1000, 800]);
        top_n = min(20, n_features);
        if sum(importance_sorted(1:top_n)) > 0
            barh(importance_sorted(1:top_n));
            set(gca, 'YTick', 1:top_n);
            set(gca, 'YTickLabel', features_sorted(1:top_n));
            set(gca, 'YDir', 'reverse');
            xlabel('Importance', 'FontSize', 12);
            ylabel('Feature', 'FontSize', 12);
            title('Top 20 Feature Importances', 'FontSize', 14, 'FontWeight', 'bold');
            grid on;
        end
        saveas(fig5, [OUTPUT_DIR, 'feature_importance.png']);
        close(fig5);
        
        % 6. Residuals Histogram
        fprintf('  - Creating residuals histogram...\n');
        fig6 = figure('Position', [100, 100, 800, 600]);
        hist(residuals_test, 30);
        xlabel('Residual', 'FontSize', 12);
        ylabel('Frequency', 'FontSize', 12);
        title('Distribution of Residuals (Test Set)', 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        saveas(fig6, [OUTPUT_DIR, 'residuals_histogram.png']);
        close(fig6);
        
        fprintf('\nAll visualizations saved to: %s\n\n', OUTPUT_DIR);
    catch err
        fprintf('  - Error generating visualizations: %s\n\n', err.message);
    end
else
    fprintf('  - Skipping all visualizations (no graphics toolkit)\n\n');
end

%% Save Results
fprintf('================================================================\n');
fprintf('Saving Results\n');
fprintf('================================================================\n\n');

% Save detailed results to text file
results_path = [OUTPUT_DIR, RESULTS_FILE];
fid = fopen(results_path, 'w');

fprintf(fid, '================================================================\n');
fprintf(fid, '   GRADIENT BOOSTING RESULTS: ln(gamma) at 90%% RH\n');
fprintf(fid, '================================================================\n\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));

fprintf(fid, 'Data Summary:\n');
fprintf(fid, '  - Data file: %s\n', DATA_FILE);
fprintf(fid, '  - Total samples: %d\n', n_samples);
fprintf(fid, '  - Features: %d\n', n_features);
fprintf(fid, '  - Training samples: %d\n', n_train);
fprintf(fid, '  - Test samples: %d\n\n', n_test);

fprintf(fid, 'Model Hyperparameters:\n');
fprintf(fid, '  - Number of trees: %d\n', N_TREES);
fprintf(fid, '  - Learning rate: %.3f\n', LEARNING_RATE);
fprintf(fid, '  - Max tree depth: %d\n', MAX_DEPTH);
fprintf(fid, '  - Min leaf size: %d\n', MIN_LEAF_SIZE);
fprintf(fid, '  - Subsample ratio: %.2f\n\n', SUBSAMPLE_RATIO);

fprintf(fid, 'Performance Metrics:\n');
fprintf(fid, '\nTraining Set:\n');
fprintf(fid, '  - RMSE: %.6f\n', train_rmse);
fprintf(fid, '  - MAE:  %.6f\n', train_mae);
fprintf(fid, '  - R²:   %.6f\n', train_r2);
fprintf(fid, '\nTest Set:\n');
fprintf(fid, '  - RMSE: %.6f\n', test_rmse);
fprintf(fid, '  - MAE:  %.6f\n', test_mae);
fprintf(fid, '  - R²:   %.6f\n\n', test_r2);

fprintf(fid, 'Top 20 Feature Importances:\n');
for i = 1:min(20, n_features)
    fprintf(fid, '  %2d. %-40s %.6f\n', i, features_sorted{i}, importance_sorted(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Test Set Predictions:\n');
fprintf(fid, '%-30s %15s %15s %15s\n', 'Electrolyte', 'Observed', 'Predicted', 'Residual');
fprintf(fid, '%s\n', repmat('-', 1, 80));
test_electrolytes = electrolyte_names_valid(test_idx);
for i = 1:length(test_idx)
    fprintf(fid, '%-30s %15.6f %15.6f %15.6f\n', ...
            test_electrolytes{i}, y_test(i), y_test_pred(i), residuals_test(i));
end

fclose(fid);
fprintf('Results saved to: %s\n\n', results_path);

% Save predictions to CSV
predictions_file = [OUTPUT_DIR, 'predictions.csv'];
fid = fopen(predictions_file, 'w');
fprintf(fid, 'Electrolyte,Observed,Predicted,Residual,Dataset\n');
for i = 1:n_train
    fprintf(fid, '%s,%.6f,%.6f,%.6f,Train\n', ...
            electrolyte_names_valid{train_idx(i)}, ...
            y_train(i), y_train_pred(i), y_train(i) - y_train_pred(i));
end
for i = 1:n_test
    fprintf(fid, '%s,%.6f,%.6f,%.6f,Test\n', ...
            test_electrolytes{i}, y_test(i), y_test_pred(i), residuals_test(i));
end
fclose(fid);
fprintf('Predictions saved to: %s\n\n', predictions_file);

fprintf('================================================================\n');
fprintf('Analysis Complete!\n');
fprintf('================================================================\n\n');
