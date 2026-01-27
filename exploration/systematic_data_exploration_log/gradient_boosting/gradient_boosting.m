%% Gradient Boosting (High Capacity & Fixed) for ln(gamma) at 90% RH
% Improvements:
%   1. FIXED: Loop iteration bug that caused "Logical Scalar" error
%   2. High Capacity: 1000 Trees, Depth 6, Low LR
%   3. Interaction Features: A*B terms
%
% Author: Generated for AWH ML Discussion
% Date: 2026-01-27

clear; clc; close all;

% Check for Toolbox
has_stats_toolbox = ~isempty(ver('stats'));

% Graphics setup
try
    graphics_toolkit('gnuplot'); 
    set(0, 'DefaultFigureVisible', 'off');
catch
end

%% Configuration
DATA_FILE = '../../../data/baseline_numeric_only.csv';
OUTPUT_DIR = '../../../figures/gradient_boosting/gradient_boosting_expanded';

% --- HYPERPARAMETERS ---
N_TREES = 1000;              % Long training
LEARNING_RATE = 0.01;        % Slow learning
MAX_DEPTH = 6;               % Deep trees
MIN_LEAF_SIZE = 5;           
SUBSAMPLE_RATIO = 0.65;      

% Feature Engineering
CREATE_INTERACTIONS = true;  

% Cross-validation
TEST_SPLIT = 0.2;
RANDOM_SEED = 999;           

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

fprintf('\n=== HIGH CAPACITY GRADIENT BOOSTING (FIXED) ===\n');

%% 1. Load Data
fprintf('Loading data...\n');

% 1.1 Read Header
fid = fopen(DATA_FILE);
header_line = fgetl(fid);
fclose(fid);
header_parts = strsplit(header_line, ',');

% 1.2 Read Electrolyte Names
fid = fopen(DATA_FILE);
fgetl(fid); 
names_cell = textscan(fid, '%s %*[^\n]', 'Delimiter', ',');
electrolyte_names = names_cell{1};
fclose(fid);

% 1.3 Read Numeric Data
try
    raw_data = dlmread(DATA_FILE, ',', 1, 1);
catch
    error('Error reading CSV. Ensure the file has a header row and text in col 1.');
end

% 1.4 Re-align Headers
data_headers = header_parts(2:end); 

% Identify Targets and Features
target_name = 'ln_gamma_at_90RH';
target_idx = find(strcmp(data_headers, target_name));

if isempty(target_idx)
    error('Target "%s" not found in headers.', target_name);
end

y_all = raw_data(:, target_idx);

% Feature Selection
exclude_patterns = {'ln_gamma', 'RH_actual', 'Dataset'}; 
feat_mask = true(1, size(raw_data, 2));

for i = 1:length(data_headers)
    if i == target_idx
        feat_mask(i) = false;
        continue;
    end
    for p = 1:length(exclude_patterns)
        if contains(data_headers{i}, exclude_patterns{p})
            feat_mask(i) = false;
        end
    end
end

X_all = raw_data(:, feat_mask);
feature_names = data_headers(feat_mask);

% Remove NaNs
valid_rows = ~isnan(y_all);
X = X_all(valid_rows, :);
y = y_all(valid_rows);
valid_names = electrolyte_names(valid_rows);

fprintf('Samples: %d | Base Features: %d\n', length(y), length(feature_names));

%% 2. Interaction Feature Engineering
if CREATE_INTERACTIONS
    fprintf('Generating interaction features (Product Terms)...\n');
    n_base = size(X, 2);
    X_inter = [];
    new_names = {};
    
    limit_cols = min(n_base, 10); 
    
    count = 0;
    for i = 1:limit_cols
        for j = i:limit_cols
            col_prod = X(:, i) .* X(:, j);
            if var(col_prod) > 1e-10
                X_inter = [X_inter, col_prod]; %#ok<AGROW>
                new_names{end+1} = [feature_names{i} '_x_' feature_names{j}];
                count = count + 1;
            end
        end
    end
    
    X = [X, X_inter];
    feature_names = [feature_names, new_names];
    fprintf('  -> Added %d interaction features. Total Features: %d\n', count, size(X, 2));
end

%% 3. Scaling
fprintf('Preprocessing...\n');
mu = mean(X);
sigma = std(X);
sigma(sigma==0) = 1;
X_scaled = (X - mu) ./ sigma;

%% 4. Train-Test Split
rand('seed', RANDOM_SEED);
n_total = length(y);
idx = randperm(n_total);
n_test = round(n_total * TEST_SPLIT);

idx_test = idx(1:n_test);
idx_train = idx(n_test+1:end);

X_train = X_scaled(idx_train, :);
y_train = y(idx_train);
X_test = X_scaled(idx_test, :);
y_test = y(idx_test);
names_test = valid_names(idx_test);

%% 5. Model Training Functions

    % --- Recursive Tree Builder ---
    function node = build_tree(X_curr, y_curr, depth, config)
        node = struct('is_leaf', true, 'prediction', mean(y_curr), ...
                      'feature', [], 'threshold', [], ...
                      'left', [], 'right', [], 'gain', 0);
        
        if depth >= config.max_depth || length(y_curr) < 2*config.min_leaf || var(y_curr) < 1e-9
            return;
        end
        
        best_gain = 0;
        best_feat = -1;
        best_thresh = 0;
        base_var = var(y_curr) * length(y_curr); 
        n_samples = length(y_curr);
        
        % Random Feature Subsampling
        n_feats = size(X_curr, 2);
        n_check = round(sqrt(n_feats)); 
        
        if config.use_toolbox
            feats_to_check = randsample(n_feats, n_check);
        else
            p = randperm(n_feats);
            feats_to_check = p(1:n_check);
        end
        
        % --- FIXED LOOP: Iterate by index to avoid vector confusion ---
        for k = 1:length(feats_to_check)
            f = feats_to_check(k); % Explicit scalar indexing
            
            vals = X_curr(:, f);
            
            % Robust Quantiles
            if length(vals) > 15
                if config.use_toolbox
                    thresholds = unique(prctile(vals, 20:20:80));
                else
                    thresholds = unique(quantile(vals, 0.2:0.2:0.8));
                end
            else
                u = unique(vals);
                if length(u) < 2, continue; end
                thresholds = (u(1:end-1) + u(2:end)) / 2;
            end

            for t_idx = 1:length(thresholds)
                t_val = thresholds(t_idx);
                mask_left = vals <= t_val;
                
                % Now mask_left is guaranteed to be a column vector
                if ~any(mask_left) || all(mask_left), continue; end
                
                n_left = sum(mask_left);
                n_right = n_samples - n_left;
                
                if n_left < config.min_leaf || n_right < config.min_leaf
                    continue;
                end
                
                y_l = y_curr(mask_left);
                y_r = y_curr(~mask_left);
                
                current_split_sse = var(y_l)*n_left + var(y_r)*n_right;
                gain = base_var - current_split_sse;
                
                if gain > best_gain
                    best_gain = gain;
                    best_feat = f;
                    best_thresh = t_val;
                end
            end
        end
        
        if best_gain > 0
            node.is_leaf = false;
            node.feature = best_feat;
            node.threshold = best_thresh;
            node.gain = best_gain;
            mask_l = X_curr(:, best_feat) <= best_thresh;
            node.left = build_tree(X_curr(mask_l,:), y_curr(mask_l), depth+1, config);
            node.right = build_tree(X_curr(~mask_l,:), y_curr(~mask_l), depth+1, config);
        end
    end

    function pred = predict_single(node, sample)
        if node.is_leaf
            pred = node.prediction;
        else
            if sample(node.feature) <= node.threshold
                pred = predict_single(node.left, sample);
            else
                pred = predict_single(node.right, sample);
            end
        end
    end

%% 6. Main Boosting Loop
fprintf('Training %d trees (Max Depth: %d, LR: %.2f)...\n', N_TREES, MAX_DEPTH, LEARNING_RATE);

tree_config.max_depth = MAX_DEPTH;
tree_config.min_leaf = MIN_LEAF_SIZE;
tree_config.use_toolbox = has_stats_toolbox;

% Initialize with mean
f0 = mean(y_train);
F_train = repmat(f0, length(y_train), 1);
F_test = repmat(f0, length(y_test), 1);

models = {};
train_rmse_hist = [];
test_rmse_hist = [];
feat_importance = zeros(size(X,2), 1);

t_start = tic;

for iter = 1:N_TREES
    residuals = y_train - F_train;
    
    n_sub = round(length(y_train) * SUBSAMPLE_RATIO);
    if has_stats_toolbox
        idx_sub = randsample(length(y_train), n_sub);
    else
        p = randperm(length(y_train));
        idx_sub = p(1:n_sub);
    end
    
    tree = build_tree(X_train(idx_sub, :), residuals(idx_sub), 0, tree_config);
    models{iter} = tree;
    
    % Prediction Update
    p_train = zeros(length(y_train), 1);
    for i=1:length(y_train), p_train(i) = predict_single(tree, X_train(i,:)); end
    
    p_test = zeros(length(y_test), 1);
    for i=1:length(y_test), p_test(i) = predict_single(tree, X_test(i,:)); end
    
    F_train = F_train + LEARNING_RATE * p_train;
    F_test = F_test + LEARNING_RATE * p_test;
    
    rmse_tr = sqrt(mean((y_train - F_train).^2));
    rmse_te = sqrt(mean((y_test - F_test).^2));
    train_rmse_hist = [train_rmse_hist; rmse_tr];
    test_rmse_hist = [test_rmse_hist; rmse_te];
    
    if ~tree.is_leaf
        feat_importance(tree.feature) = feat_importance(tree.feature) + tree.gain;
    end

    if mod(iter, 50) == 0
        fprintf('Iter %4d | Train RMSE: %.4f | Test RMSE: %.4f\n', iter, rmse_tr, rmse_te);
    end
end
toc(t_start);

%% 7. Evaluation & Visualization

residuals_final = y_test - F_test;
r2_test = 1 - sum(residuals_final.^2) / sum((y_test - mean(y_test)).^2);
fprintf('\nFinal Test R^2: %.4f\n', r2_test);

% --- Plot 1: Learning Curve ---
fig1 = figure;
plot(train_rmse_hist, 'b-', 'LineWidth', 2); hold on;
plot(test_rmse_hist, 'r-', 'LineWidth', 2);
legend('Train', 'Test'); xlabel('Trees'); ylabel('RMSE');
title(sprintf('GBM Learning Curve (Best Test RMSE: %.4f)', min(test_rmse_hist))); grid on;
saveas(fig1, [OUTPUT_DIR 'learning_curve.png']);

% --- Plot 2: Obs vs Pred ---
fig2 = figure('Position', [100, 100, 1000, 500]);
subplot(1,2,1);
scatter(y_train, F_train, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5); hold on;
plot([min(y_all), max(y_all)], [min(y_all), max(y_all)], 'k--');
xlabel('Observed'); ylabel('Predicted'); title('Training Set'); grid on; axis square;

subplot(1,2,2);
scatter(y_test, F_test, 40, 'r', 'filled', 'MarkerFaceAlpha', 0.6); hold on;
plot([min(y_all), max(y_all)], [min(y_all), max(y_all)], 'k--');
xlabel('Observed'); ylabel('Predicted'); 
title(sprintf('Test Set (R^2 = %.2f)', r2_test)); grid on; axis square;
saveas(fig2, [OUTPUT_DIR 'obs_vs_pred_final.png']);

% --- Plot 3: Feature Importance ---
[sorted_imp, idx_imp] = sort(feat_importance, 'descend');
top_k = min(20, length(sorted_imp));

fig3 = figure;
if top_k > 0
    barh(sorted_imp(1:top_k));
    set(gca, 'YTickLabel', feature_names(idx_imp(1:top_k)), 'YDir', 'reverse', 'TickLabelInterpreter', 'none');
    xlabel('Split Gain');
    title('Top Feature Importances (Inc. Interactions)');
end
saveas(fig3, [OUTPUT_DIR 'importance.png']);

% Save Predictions
fid_csv = fopen([OUTPUT_DIR 'predictions_test.csv'], 'w');
fprintf(fid_csv, 'Electrolyte,Observed,Predicted,Residual\n');
for i = 1:length(y_test)
    fprintf(fid_csv, '%s,%.6f,%.6f,%.6f\n', ...
            names_test{i}, y_test(i), F_test(i), residuals_final(i));
end
fclose(fid_csv);
fprintf('Predictions saved to %s\n', [OUTPUT_DIR 'predictions_test.csv']);