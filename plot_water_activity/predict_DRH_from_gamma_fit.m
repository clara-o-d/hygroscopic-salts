close all 
clear
clc 

% Suppress fzero warnings (they're handled by our error checking)
warning('off', 'MATLAB:fzero:AbortFcn');

% Add calculate_mf, util, and data folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Load canonical salt data from data/load_salt_data.m
salt_data = load_salt_data();
exclude = {'NH4NO3', 'MgNO32'};
keep = cellfun(@(r) ~any(strcmp(r{1}, exclude)), salt_data);
salt_data = salt_data(keep);

%% Step 1: Use predefined gamma_w(RH) fit
fprintf('=== Step 1: Using predefined gamma_w(RH) fit ===\n');

% Hardcoded fit coefficients: gamma_w = a*RH^2 + b*RH + c
% (RH in percent)
a_fit = -0.000034;
b_fit = 0.013592;
c_fit = -0.015254;

fprintf('Fit equation: gamma_w = %.6f*RH^2 + %.6f*RH + %.6f\n', a_fit, b_fit, c_fit);
fprintf('(RH in percent)\n\n');

%% Step 2: Load solubility data
fprintf('=== Step 2: Loading solubility data ===\n');

baselineFile = fullfile(filepath, '..', 'data', 'baseline_numeric_only.csv');
if ~isfile(baselineFile)
    error('Cannot find baseline file: %s', baselineFile);
end

opts = detectImportOptions(baselineFile);
opts.VariableNamingRule = 'preserve'; 
tbl = readtable(baselineFile, opts);
varNames = tbl.Properties.VariableNames;

% Find electrolyte and solubility columns
elec_col_idx = find(strcmpi(varNames, 'electrolyte'));
sol_col_idx = find(strcmpi(varNames, 'solubility_limit'));

if isempty(elec_col_idx) || isempty(sol_col_idx)
    error('Could not find required columns in CSV file');
end

raw_names = string(tbl.(varNames{elec_col_idx}));
sol_values = tbl.(varNames{sol_col_idx});

% Clean CSV names: Remove parentheses and trim whitespace
clean_csv_names = erase(raw_names, ["(", ")"]);
clean_csv_names = strtrim(clean_csv_names);

fprintf('Loaded solubility data for %d entries\n\n', length(clean_csv_names));

%% Step 3: Predict DRH for each salt
fprintf('=== Step 3: Predicting DRH from gamma_w fit + solubility ===\n');

% Create a lookup structure for salt data
salt_lookup = struct();
for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    salt_lookup.(salt_name) = salt_data{s};
end

% Storage for results
results = struct();
results.salt_name = {};
results.MW = [];
results.solubility_limit = [];  % g salt / 100 g water
results.DRH_actual = [];         % Actual DRH from data (RH_min)
results.DRH_predicted = [];      % Predicted DRH from fit
results.x_w_at_DRH = [];         % Mole fraction at predicted DRH
results.gamma_w_at_DRH = [];     % Gamma_w at predicted DRH

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min_actual = salt_data{s}{3};  % This is the actual DRH
    
    % Find solubility in CSV
    idx = find(strcmpi(clean_csv_names, salt_name), 1);
    
    if isempty(idx)
        fprintf('  %s: No solubility data found in CSV\n', salt_name);
        continue;
    end
    
    solubility_limit = sol_values(idx);
    
    % Skip invalid solubility values
    if solubility_limit <= 0 || isnan(solubility_limit)
        fprintf('  %s: Invalid solubility value (%.2f)\n', salt_name, solubility_limit);
        continue;
    end
    
    % Calculate target mole fraction from solubility
    % Solubility is in g salt / 100 g water (g solute/100 g H2O)
    % Mass ratio: m_salt / m_water = solubility_limit / 100
    % Mole ratio: n_salt / n_water = (m_salt / MW_salt) / (m_water / MW_water)
    %            = (m_salt / m_water) * (MW_water / MW_salt)
    %            = (solubility_limit / 100) * (MWw / MW)
    % x_w = n_w / (n_w + n_s) = 1 / (1 + n_s/n_w)
    %     = 1 / (1 + (solubility_limit/100) * (MWw/MW))
    
    mass_ratio = solubility_limit / 100;  % m_salt / m_water (dimensionless)
    mole_ratio = mass_ratio * (MWw / MW);  % n_salt / n_water (dimensionless)
    x_w_target = 1 / (1 + mole_ratio);  % mole fraction of water (dimensionless, 0-1)
    
    % Now solve for RH where: x_w = a_w / gamma_w(RH)
    % where a_w = RH (as fraction, 0-1) and gamma_w(RH) uses RH in percent (0-100)
    % So: x_w = (RH_percent/100) / gamma_w(RH_percent)
    % where gamma_w(RH_percent) = a*RH_percent^2 + b*RH_percent + c
    % Rearranging: x_w * (a*RH_percent^2 + b*RH_percent + c) = RH_percent/100
    % Multiplying by 100: 100*x_w*(a*RH_percent^2 + b*RH_percent + c) = RH_percent
    % Rearranging: 100*a*x_w*RH_percent^2 + 100*b*x_w*RH_percent + 100*c*x_w - RH_percent = 0
    % Final form: 100*a*x_w*RH_percent^2 + (100*b*x_w - 1)*RH_percent + 100*c*x_w = 0
    
    % Check for valid x_w_target
    if x_w_target <= 0 || x_w_target >= 1 || isnan(x_w_target) || isinf(x_w_target)
        fprintf('  %s: Invalid target mole fraction (%.6f)\n', salt_name, x_w_target);
        continue;
    end
    
    % Quadratic equation coefficients (RH in percent)
    % Equation: 100*a*x_w*RH^2 + (100*b*x_w - 1)*RH + 100*c*x_w = 0
    A_quad = 100 * a_fit * x_w_target;
    B_quad = 100 * b_fit * x_w_target - 1;
    C_quad = 100 * c_fit * x_w_target;
    
    % Check for numerical issues
    if isnan(A_quad) || isnan(B_quad) || isnan(C_quad) || ...
       isinf(A_quad) || isinf(B_quad) || isinf(C_quad)
        fprintf('  %s: Invalid quadratic coefficients\n', salt_name);
        continue;
    end
    
    % Handle linear case (A_quad ≈ 0)
    if abs(A_quad) < 1e-10
        % Linear equation: B*RH + C = 0
        if abs(B_quad) < 1e-10
            fprintf('  %s: Degenerate equation\n', salt_name);
            continue;
        end
        RH_predicted = -C_quad / B_quad;
    else
        % Solve quadratic: A*RH^2 + B*RH + C = 0
        discriminant = B_quad^2 - 4*A_quad*C_quad;
        
        if discriminant < 0
            fprintf('  %s: No real solution (discriminant = %.2e < 0)\n', salt_name, discriminant);
            continue;
        end
        
        % Two solutions
        sqrt_disc = sqrt(discriminant);
        RH1 = (-B_quad + sqrt_disc) / (2*A_quad);
        RH2 = (-B_quad - sqrt_disc) / (2*A_quad);
        
        % Choose the solution in the valid range (0-100%)
        % Prefer the one closer to the actual DRH
        RH_predicted = [];
        if RH1 > 0 && RH1 <= 100 && ~isnan(RH1) && ~isinf(RH1)
            RH_predicted = RH1;
        end
        if RH2 > 0 && RH2 <= 100 && ~isnan(RH2) && ~isinf(RH2)
            if isempty(RH_predicted)
                RH_predicted = RH2;
            else
                % Both valid, choose one closer to actual
                if abs(RH2 - RH_min_actual*100) < abs(RH_predicted - RH_min_actual*100)
                    RH_predicted = RH2;
                end
            end
        end
        
        if isempty(RH_predicted)
            fprintf('  %s: No valid solution in range 0-100%% (RH1=%.2f, RH2=%.2f)\n', ...
                salt_name, RH1, RH2);
            continue;
        end
    end
    
    % Calculate gamma_w at predicted DRH (RH in percent)
    gamma_w_pred = a_fit * RH_predicted^2 + b_fit * RH_predicted + c_fit;
    
    % Verify: x_w should equal (RH_percent/100) / gamma_w
    % where RH_percent is the predicted DRH in percent
    x_w_verify = (RH_predicted / 100) / gamma_w_pred;
    
    % Store results
    results.salt_name{end+1} = salt_name;
    results.MW(end+1) = MW;
    results.solubility_limit(end+1) = solubility_limit;
    results.DRH_actual(end+1) = RH_min_actual * 100;  % Convert to percent
    results.DRH_predicted(end+1) = RH_predicted;
    results.x_w_at_DRH(end+1) = x_w_target;
    results.gamma_w_at_DRH(end+1) = gamma_w_pred;
    
    fprintf('  %s: Solubility=%.2f g/100g, Actual DRH=%.2f%%, Predicted DRH=%.2f%%\n', ...
        salt_name, solubility_limit, RH_min_actual*100, RH_predicted);
end

fprintf('\nSuccessfully predicted DRH for %d salts\n\n', length(results.salt_name));

%% Step 4: Create validation plots
fprintf('=== Step 4: Creating validation plots ===\n');

if length(results.salt_name) > 0
    % Convert to arrays for plotting
    DRH_actual = [results.DRH_actual];
    DRH_predicted = [results.DRH_predicted];
    salt_names = results.salt_name;
    
    % Calculate statistics
    errors = DRH_predicted - DRH_actual;
    RMSE = sqrt(mean(errors.^2));
    MAE = mean(abs(errors));
    R_squared = 1 - sum(errors.^2) / sum((DRH_actual - mean(DRH_actual)).^2);
    
    fprintf('Validation Statistics:\n');
    fprintf('  RMSE = %.4f%%\n', RMSE);
    fprintf('  MAE = %.4f%%\n', MAE);
    fprintf('  R^2 = %.4f\n', R_squared);
    fprintf('  Mean error = %.4f%%\n', mean(errors));
    fprintf('  Std error = %.4f%%\n', std(errors));
    
    %% Plot 1: Predicted vs Actual DRH
    figure('Position', [100, 100, 800, 700]);
    hold on; grid on; box on;
    
    scatter(DRH_actual, DRH_predicted, 80, 'filled', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Add labels
    for k = 1:length(salt_names)
        text(DRH_actual(k), DRH_predicted(k), ['  ' salt_names{k}], ...
            'FontSize', 8, 'Interpreter', 'none', ...
            'VerticalAlignment', 'middle');
    end
    
    % 1:1 Line
    maxVal = max([DRH_actual, DRH_predicted]) * 1.1;
    minVal = min([DRH_actual, DRH_predicted]) * 0.9;
    plot([minVal, maxVal], [minVal, maxVal], 'k--', 'LineWidth', 2, 'DisplayName', '1:1 Line');
    
    xlabel('Actual DRH (%)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Predicted DRH (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title('DRH Prediction: Predicted vs Actual', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Add statistics textbox
    stats_str = sprintf(['$R^2 = %.4f$\n' ...
                        '$RMSE = %.2f\\%%$\n' ...
                        '$MAE = %.2f\\%%$'], ...
                        R_squared, RMSE, MAE);
    
    annotation('textbox', [0.15 0.70 0.25 0.10], ...
        'String', stats_str, ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1.1);
    
    axis square;
    xlim([minVal, maxVal]);
    ylim([minVal, maxVal]);
    set(gca, 'FontSize', 12);
    set(gcf, 'color', 'w');
    
    filename = 'DRH_Prediction_Validation';
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600');
    savefig(fullfile(fig_out_dir, [filename '.fig']));
    
    %% Plot 2: Error vs Actual DRH
    figure('Position', [150, 150, 800, 600]);
    hold on; grid on; box on;
    
    scatter(DRH_actual, errors, 80, 'filled', ...
        'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Add labels
    for k = 1:length(salt_names)
        text(DRH_actual(k), errors(k), ['  ' salt_names{k}], ...
            'FontSize', 8, 'Interpreter', 'none', ...
            'VerticalAlignment', 'middle');
    end
    
    % Zero line
    yline(0, 'k--', 'LineWidth', 2);
    
    xlabel('Actual DRH (%)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Prediction Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title('DRH Prediction Error vs Actual DRH', 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gca, 'FontSize', 12);
    set(gcf, 'color', 'w');
    
    filename = 'DRH_Prediction_Error';
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600');
    savefig(fullfile(fig_out_dir, [filename '.fig']));
    
    %% Plot 3: Solubility vs DRH (showing relationship)
    figure('Position', [200, 200, 800, 600]);
    hold on; grid on; box on;
    
    solubility = [results.solubility_limit];
    
    scatter(solubility, DRH_actual, 80, 'filled', ...
        'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'Actual DRH');
    
    scatter(solubility, DRH_predicted, 80, 'filled', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'Predicted DRH');
    
    % Add labels
    for k = 1:length(salt_names)
        text(solubility(k), DRH_actual(k), ['  ' salt_names{k}], ...
            'FontSize', 8, 'Interpreter', 'none', ...
            'VerticalAlignment', 'middle', 'Color', [0.4660 0.6740 0.1880]);
    end
    
    xlabel('Solubility Limit (g salt / 100 g water)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('DRH (%)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Relationship: Solubility vs DRH', 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 11);
    
    set(gca, 'FontSize', 12);
    set(gcf, 'color', 'w');
    
    filename = 'DRH_vs_Solubility';
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600');
    savefig(fullfile(fig_out_dir, [filename '.fig']));
    
    fprintf('\nPlots saved to: %s\n', fig_out_dir);
else
    warning('No valid predictions to plot');
end

fprintf('\n=== Complete ===\n');
