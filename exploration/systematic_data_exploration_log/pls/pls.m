% Partial Least Squares Regression Analysis for Multiple RH Levels
% Target: ln_gamma at different RH levels (30%, 50%, 70%, 90%)
% Dataset: baseline_numeric_only.csv

clear all; close all; clc;

%% Load Data
data_file = '../../data/baseline_numeric_only.csv';
fprintf('Loading data from: %s\n', data_file);

% Try MATLAB's readtable first, fall back to manual parsing for Octave
try
    % MATLAB approach
    data_table = readtable(data_file);
    use_table = true;
    fprintf('Data loaded using readtable (MATLAB)\n');
catch
    % Octave/manual approach
    use_table = false;
    fprintf('readtable not available, using manual parsing\n');
end

%% Process data based on loading method
if use_table
    % MATLAB path using table
    fprintf('Data shape: %d rows x %d columns\n', height(data_table), width(data_table));
    
    % Get all column names
    all_cols = data_table.Properties.VariableNames;
    
    % Extract all data for processing
    data_struct = struct();
    for i = 1:length(all_cols)
        data_struct.(all_cols{i}) = data_table.(all_cols{i});
    end
    
else
    % Octave path - manual CSV parsing
    % Read header
    fid = fopen(data_file, 'r');
    header_line = fgetl(fid);
    fclose(fid);
    
    % Parse header
    header = strsplit(header_line, ',');
    
    % Read all lines manually
    fid = fopen(data_file, 'r');
    fgetl(fid); % Skip header
    
    data_cell = {};
    row = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            row = row + 1;
            data_cell{row} = strsplit(line, ',');
        end
    end
    fclose(fid);
    
    n_rows = row;
    n_cols = length(header);
    
    fprintf('Data shape: %d rows x %d columns\n', n_rows, n_cols);
    
    % Build data structure
    data_struct = struct();
    all_cols = header;
    
    for col = 1:n_cols
        col_name = header{col};
        col_data = zeros(n_rows, 1);
        
        for r = 1:n_rows
            val_str = data_cell{r}{col};
            if ~isempty(val_str)
                col_data(r) = str2double(val_str);
            else
                col_data(r) = NaN;
            end
        end
        
        data_struct.(col_name) = col_data;
    end
end

%% Define RH levels to analyze
rh_levels = [90, 70, 50, 30];
target_cols = {'ln_gamma_at_90RH', 'ln_gamma_at_70RH', 'ln_gamma_at_50RH', 'ln_gamma_at_30RH'};

% Define columns to exclude from features
exclude_col_names = {'electrolyte', 'molecule_formula', ...
                     'ln_gamma_at_90RH', 'RH_actual_for_ln_gamma_90', ...
                     'ln_gamma_at_70RH', 'RH_actual_for_ln_gamma_70', ...
                     'ln_gamma_at_50RH', 'RH_actual_for_ln_gamma_50', ...
                     'ln_gamma_at_30RH', 'RH_actual_for_ln_gamma_30'};

% Identify feature columns
feature_cols = {};
for i = 1:length(all_cols)
    col_name = all_cols{i};
    if ~ismember(col_name, exclude_col_names)
        % Check if it's numeric and not all NaN
        col_data = data_struct.(col_name);
        if isnumeric(col_data) && ~all(isnan(col_data))
            feature_cols{end+1} = col_name;
        end
    end
end

fprintf('\nTotal feature columns identified: %d\n', length(feature_cols));

%% Run PLS for each RH level
results_all = struct();

% Check if plsregress is available (MATLAB)
has_plsregress = exist('plsregress', 'file') == 2;

if has_plsregress
    fprintf('\nUsing MATLAB plsregress function\n');
else
    fprintf('\nUsing NIPALS algorithm (Octave compatible)\n');
end

for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    target_col = target_cols{rh_idx};
    
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('   ANALYSIS FOR %d%% RH\n', rh);
    fprintf('========================================\n');
    
    % Get target data
    if ~isfield(data_struct, target_col)
        fprintf('Warning: Target column %s not found. Skipping.\n', target_col);
        continue;
    end
    
    target = data_struct.(target_col);
    
    % Find valid rows (non-NaN target)
    valid_idx = ~isnan(target);
    n_valid = sum(valid_idx);
    
    fprintf('Samples with valid data: %d\n', n_valid);
    
    if n_valid < 10
        fprintf('Warning: Too few samples (%d). Skipping.\n', n_valid);
        continue;
    end
    
    % Extract target for valid samples
    y = target(valid_idx);
    
    % Build feature matrix for valid samples
    X = [];
    valid_feature_names = {};
    
    for i = 1:length(feature_cols)
        col_data = data_struct.(feature_cols{i});
        col_data_valid = col_data(valid_idx);
        
        % Only include if no NaN in this subset
        if ~any(isnan(col_data_valid))
            X = [X, col_data_valid];
            valid_feature_names{end+1} = feature_cols{i};
        end
    end
    
    fprintf('Features used: %d\n', size(X, 2));
    
    if size(X, 2) < 5
        fprintf('Warning: Too few features. Skipping.\n');
        continue;
    end
    
    % Standardize features
    X_mean = mean(X);
    X_std = std(X);
    X_std(X_std == 0) = 1;
    X_scaled = bsxfun(@minus, X, X_mean);
    X_scaled = bsxfun(@rdivide, X_scaled, X_std);
    
    % Standardize target
    y_mean = mean(y);
    y_std = std(y);
    y_scaled = (y - y_mean) / y_std;
    
    % Determine number of components
    max_components = min([7, size(X_scaled, 2), size(X_scaled, 1) - 1]);
    opt_ncomp = max_components;
    
    fprintf('Using %d PLS components\n', opt_ncomp);
    
    % Run PLS
    if has_plsregress
        % MATLAB path
        [XL, YL, XS, YS, BETA, PCTVAR_final, ~, stats] = plsregress(X_scaled, y_scaled, opt_ncomp);
        
        % Make predictions
        y_pred_scaled = [ones(size(X_scaled, 1), 1), X_scaled] * BETA;
        y_pred = y_pred_scaled * y_std + y_mean;
        
    else
        % Octave path - NIPALS algorithm
        n_samples = size(X_scaled, 1);
        n_features = size(X_scaled, 2);
        
        XL = zeros(n_features, opt_ncomp);
        YL = zeros(1, opt_ncomp);
        XS = zeros(n_samples, opt_ncomp);
        W = zeros(n_features, opt_ncomp);
        
        X_temp = X_scaled;
        y_temp = y_scaled;
        
        for comp = 1:opt_ncomp
            w = X_temp' * y_temp;
            w = w / norm(w);
            
            for iter = 1:100
                w_old = w;
                t = X_temp * w;
                q = y_temp' * t / (t' * t);
                w = X_temp' * (t * q) / (t' * t * q^2);
                w = w / norm(w);
                if norm(w - w_old) < 1e-6
                    break;
                end
            end
            
            t = X_temp * w;
            p = X_temp' * t / (t' * t);
            q = y_temp' * t / (t' * t);
            
            W(:, comp) = w;
            XS(:, comp) = t;
            XL(:, comp) = p;
            YL(comp) = q;
            
            X_temp = X_temp - t * p';
            y_temp = y_temp - t * q;
        end
        
        % Calculate regression coefficients
        P = XL;
        Q = YL';
        B_scaled = W * inv(P' * W) * Q;
        
        B0 = y_mean - (X_mean ./ X_std) * B_scaled * y_std;
        B = (B_scaled * y_std) ./ X_std';
        
        y_pred = X * B + B0;
        
        stats = struct();
        stats.W = W;
        PCTVAR_final = [];
    end
    
    % Calculate performance metrics
    SS_tot = sum((y - mean(y)).^2);
    SS_res = sum((y - y_pred).^2);
    R2 = 1 - SS_res / SS_tot;
    RMSE = sqrt(mean((y - y_pred).^2));
    MAE = mean(abs(y - y_pred));
    
    fprintf('\nModel Performance:\n');
    fprintf('  R² = %.4f\n', R2);
    fprintf('  RMSE = %.4f\n', RMSE);
    fprintf('  MAE = %.4f\n', MAE);
    
    % Calculate VIP scores
    W0 = stats.W;
    for i = 1:size(W0, 2)
        W0(:, i) = W0(:, i) / norm(W0(:, i));
    end
    
    SS = sum(XS.^2, 1) .* (YL.^2);
    VIP = sqrt(size(XL, 2) * sum(bsxfun(@times, SS, W0.^2), 2) / sum(SS));
    
    [VIP_sorted, sort_idx] = sort(VIP, 'descend');
    
    % Display top VIP features
    fprintf('\nTop 10 Features by VIP:\n');
    fprintf('  %-45s %10s\n', 'Feature', 'VIP Score');
    fprintf('  %s\n', repmat('-', 1, 58));
    for i = 1:min(10, length(VIP_sorted))
        marker = '';
        if VIP_sorted(i) > 1
            marker = ' *';
        end
        fprintf('  %-45s %10.4f%s\n', valid_feature_names{sort_idx(i)}, VIP_sorted(i), marker);
    end
    
    % Component breakdown
    fprintf('\n  Component Breakdown:\n');
    
    if has_plsregress && ~isempty(PCTVAR_final)
        X_var_explained = PCTVAR_final(1, 1:opt_ncomp)';
        Y_var_explained = PCTVAR_final(2, 1:opt_ncomp)';
    else
        X_var_explained = zeros(opt_ncomp, 1);
        Y_var_explained = zeros(opt_ncomp, 1);
        
        for comp = 1:opt_ncomp
            X_var_explained(comp) = 100 * sum(XS(:, comp).^2) / sum(sum(X_scaled.^2));
            Y_var_explained(comp) = 100 * sum((XS(:, comp) * YL(comp)).^2) / sum(y_scaled.^2);
        end
    end
    
    Y_cumvar = cumsum(Y_var_explained);
    
    fprintf('  %-10s %15s %15s\n', 'Component', 'Y Var (%%)', 'Y Cumul (%%)');
    fprintf('  %s\n', repmat('-', 1, 43));
    for comp = 1:opt_ncomp
        fprintf('  %-10d %15.2f %15.2f\n', comp, Y_var_explained(comp), Y_cumvar(comp));
    end
    
    % Store results
    results_all.(sprintf('RH_%d', rh)) = struct();
    results_all.(sprintf('RH_%d', rh)).n_samples = n_valid;
    results_all.(sprintf('RH_%d', rh)).n_features = size(X, 2);
    results_all.(sprintf('RH_%d', rh)).feature_names = valid_feature_names;
    results_all.(sprintf('RH_%d', rh)).R2 = R2;
    results_all.(sprintf('RH_%d', rh)).RMSE = RMSE;
    results_all.(sprintf('RH_%d', rh)).MAE = MAE;
    results_all.(sprintf('RH_%d', rh)).VIP = VIP;
    results_all.(sprintf('RH_%d', rh)).VIP_sorted = VIP_sorted;
    results_all.(sprintf('RH_%d', rh)).VIP_idx = sort_idx;
    results_all.(sprintf('RH_%d', rh)).Y_var_explained = Y_var_explained;
    results_all.(sprintf('RH_%d', rh)).Y_cumvar = Y_cumvar;
    results_all.(sprintf('RH_%d', rh)).XL = XL;
    results_all.(sprintf('RH_%d', rh)).YL = YL;
    results_all.(sprintf('RH_%d', rh)).XS = XS;
    results_all.(sprintf('RH_%d', rh)).W = stats.W;
    results_all.(sprintf('RH_%d', rh)).y = y;
    results_all.(sprintf('RH_%d', rh)).y_pred = y_pred;
    results_all.(sprintf('RH_%d', rh)).opt_ncomp = opt_ncomp;
end

%% Summary comparison across RH levels
fprintf('\n');
fprintf('========================================\n');
fprintf('   SUMMARY ACROSS RH LEVELS\n');
fprintf('========================================\n\n');

fprintf('%-10s %10s %10s %10s %10s %10s\n', 'RH Level', 'Samples', 'Features', 'R²', 'RMSE', 'MAE');
fprintf('%s\n', repmat('-', 1, 62));

for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    field_name = sprintf('RH_%d', rh);
    
    if isfield(results_all, field_name)
        res = results_all.(field_name);
        fprintf('%-10s %10d %10d %10.4f %10.4f %10.4f\n', ...
                sprintf('%d%%', rh), res.n_samples, res.n_features, ...
                res.R2, res.RMSE, res.MAE);
    end
end

%% Save comprehensive results
fprintf('\n=== Saving Results ===\n');

fid = fopen('pls_results_all_RH.txt', 'w');
fprintf(fid, 'Partial Least Squares Regression Results\n');
fprintf(fid, 'Multi-RH Analysis\n');
fprintf(fid, '=========================================\n\n');
fprintf(fid, 'Dataset: %s\n\n', data_file);

% Summary table
fprintf(fid, 'Summary Across RH Levels:\n');
fprintf(fid, '-------------------------\n');
fprintf(fid, '%-10s %10s %10s %10s %10s %10s\n', 'RH Level', 'Samples', 'Features', 'R²', 'RMSE', 'MAE');
fprintf(fid, '%s\n', repmat('-', 1, 62));

for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    field_name = sprintf('RH_%d', rh);
    
    if isfield(results_all, field_name)
        res = results_all.(field_name);
        fprintf(fid, '%-10s %10d %10d %10.4f %10.4f %10.4f\n', ...
                sprintf('%d%%', rh), res.n_samples, res.n_features, ...
                res.R2, res.RMSE, res.MAE);
    end
end

% Detailed results for each RH
for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    field_name = sprintf('RH_%d', rh);
    
    if ~isfield(results_all, field_name)
        continue;
    end
    
    res = results_all.(field_name);
    
    fprintf(fid, '\n\n');
    fprintf(fid, '=========================================\n');
    fprintf(fid, '   DETAILED RESULTS FOR %d%% RH\n', rh);
    fprintf(fid, '=========================================\n\n');
    
    fprintf(fid, 'Samples: %d\n', res.n_samples);
    fprintf(fid, 'Features: %d\n', res.n_features);
    fprintf(fid, 'Components: %d\n\n', res.opt_ncomp);
    
    fprintf(fid, 'Model Performance:\n');
    fprintf(fid, '  R² = %.4f\n', res.R2);
    fprintf(fid, '  RMSE = %.4f\n', res.RMSE);
    fprintf(fid, '  MAE = %.4f\n\n', res.MAE);
    
    fprintf(fid, 'Component Variance Explained:\n');
    fprintf(fid, '%-10s %15s %15s\n', 'Component', 'Y Var (%%)', 'Y Cumul (%%)');
    fprintf(fid, '%s\n', repmat('-', 1, 43));
    for comp = 1:res.opt_ncomp
        fprintf(fid, '%-10d %15.2f %15.2f\n', ...
                comp, res.Y_var_explained(comp), res.Y_cumvar(comp));
    end
    
    fprintf(fid, '\n\nVariable Importance (VIP > 1.0):\n');
    fprintf(fid, '%-45s %10s\n', 'Feature', 'VIP Score');
    fprintf(fid, '%s\n', repmat('-', 1, 58));
    
    for i = 1:length(res.VIP_sorted)
        if res.VIP_sorted(i) > 1.0
            fprintf(fid, '%-45s %10.4f\n', ...
                    res.feature_names{res.VIP_idx(i)}, res.VIP_sorted(i));
        end
    end
    
    fprintf(fid, '\n\nTop Features per Component:\n');
    for comp = 1:res.opt_ncomp
        fprintf(fid, '\nComponent %d (%.1f%% Y variance):\n', comp, res.Y_var_explained(comp));
        
        comp_loadings = abs(res.XL(:, comp));
        [~, load_idx] = sort(comp_loadings, 'descend');
        
        fprintf(fid, '  %-45s %12s %12s\n', 'Feature', 'Loading', 'Weight');
        fprintf(fid, '  %s\n', repmat('-', 1, 72));
        
        for i = 1:min(5, length(res.feature_names))
            idx = load_idx(i);
            fprintf(fid, '  %-45s %12.4f %12.4f\n', ...
                    res.feature_names{idx}, res.XL(idx, comp), res.W(idx, comp));
        end
    end
end

fclose(fid);
fprintf('Comprehensive results saved to pls_results_all_RH.txt\n');

% Save VIP comparison across RH levels
fid = fopen('pls_vip_comparison.csv', 'w');
fprintf(fid, 'Feature');
for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    field_name = sprintf('RH_%d', rh);
    if isfield(results_all, field_name)
        fprintf(fid, ',VIP_%dRH', rh);
    end
end
fprintf(fid, '\n');

% Get all unique features
all_features = {};
for rh_idx = 1:length(rh_levels)
    rh = rh_levels(rh_idx);
    field_name = sprintf('RH_%d', rh);
    if isfield(results_all, field_name)
        for f = results_all.(field_name).feature_names
            if ~ismember(f{1}, all_features)
                all_features{end+1} = f{1};
            end
        end
    end
end

for feat = all_features
    fprintf(fid, '%s', feat{1});
    for rh_idx = 1:length(rh_levels)
        rh = rh_levels(rh_idx);
        field_name = sprintf('RH_%d', rh);
        if isfield(results_all, field_name)
            res = results_all.(field_name);
            feat_idx = find(strcmp(res.feature_names, feat{1}));
            if ~isempty(feat_idx)
                fprintf(fid, ',%.6f', res.VIP(feat_idx));
            else
                fprintf(fid, ',');
            end
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('VIP comparison saved to pls_vip_comparison.csv\n');

fprintf('\n=== Analysis Complete ===\n');

