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

% Process each salt and collect data for fitting
num_points = 100;
all_salts_data = struct();

% Initialize container arrays for fitting
all_RH_fit = [];
all_gamma_fit = [];
all_weights = [];

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    
    % Get Ion counts
    n_cat = salt_data{s}{12};
    n_an  = salt_data{s}{13};
    nu = n_cat + n_an; % Total number of dissociated ions
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water_mol = zeros(size(RH_vec));
    
    % Calculate for each RH value
    valid_indices = [];
    
    for i = 1:length(RH_vec)
        try
            if func_args == 1
                mf_salt(i) = feval(func_name, RH_vec(i), T);
            else
                mf_salt(i) = feval(func_name, RH_vec(i));
            end
            
            % Check for NaN or Inf
            if isnan(mf_salt(i)) || isinf(mf_salt(i))
                continue;  % Skip this point
            end
            
            mf_water(i) = 1 - mf_salt(i);
            
            % Check for valid mass fractions
            if mf_water(i) <= 0 || mf_water(i) >= 1
                continue;  % Skip this point
            end
            
            % Moles of water and salt (per unit mass basis)
            n_w = mf_water(i) / MWw;
            n_s = mf_salt(i) / MW;
            
            % Molecular Definition: x_w = n_w / (n_w + n_s)
            x_water_mol(i) = n_w / (n_w + n_s);
            
            % Check for valid mole fraction
            if isnan(x_water_mol(i)) || isinf(x_water_mol(i)) || x_water_mol(i) <= 0 || x_water_mol(i) >= 1
                continue;  % Skip this point
            end
            
            valid_indices(end+1) = i;
            
        catch ME
            % Silently skip failed points
            continue;
        end
    end
    
    % Only proceed if we have valid data points
    if length(valid_indices) >= 10  % Need at least 10 valid points
        % Extract valid data
        RH_valid = RH_vec(valid_indices);
        x_water_valid = x_water_mol(valid_indices);
        
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        gamma_w_mol = RH_valid ./ x_water_valid;
        
        % Check for NaN/Inf in gamma
        valid_gamma = ~isnan(gamma_w_mol) & ~isinf(gamma_w_mol) & gamma_w_mol > 0 & gamma_w_mol < 10;
        if sum(valid_gamma) >= 10
            RH_valid = RH_valid(valid_gamma);
            gamma_w_mol = gamma_w_mol(valid_gamma);
            x_water_valid = x_water_valid(valid_gamma);
            
            % Store data
            all_salts_data.(salt_name) = struct();
            all_salts_data.(salt_name).RH = RH_valid;
            all_salts_data.(salt_name).x_water_mol = x_water_valid;
            all_salts_data.(salt_name).gamma_w_mol = gamma_w_mol;
            all_salts_data.(salt_name).display_name = salt_name;
            
            % Calculate the Range of RH for this salt to use as weight
            rh_range = max(RH_valid) - min(RH_valid);
            
            % Assign weight (range) to every point in this dataset
            w_vec = ones(size(RH_valid)) * rh_range;
            
            % Append to master lists (Convert RH to %)
            all_RH_fit = [all_RH_fit, RH_valid * 100]; 
            all_gamma_fit = [all_gamma_fit, gamma_w_mol];
            all_weights = [all_weights, w_vec];
        end
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

%% Calculate Constrained Weighted Polynomial Fit
% Quadratic Model: y = a*x^2 + b*x + c
% Constraint: Pass through (100%, 1) -> 1 = a(10000) + b(100) + c
% Solve for c: c = 1 - 10000a - 100b
% Substitute c back into model: 
% y = a*x^2 + b*x + (1 - 10000a - 100b)
% Rearrange to isolate a and b:
% y - 1 = a(x^2 - 10000) + b(x - 100)

% Filter out any remaining NaN/Inf values before fitting
valid_mask = ~isnan(all_gamma_fit) & ~isinf(all_gamma_fit) & ...
             ~isnan(all_RH_fit) & ~isinf(all_RH_fit) & ...
             all_gamma_fit > 0 & all_gamma_fit < 10 & ...
             all_RH_fit >= 0 & all_RH_fit <= 100;

if sum(valid_mask) < 50
    error('Not enough valid data points for fitting. Only %d valid points found.', sum(valid_mask));
end

all_RH_fit_clean = all_RH_fit(valid_mask);
all_gamma_fit_clean = all_gamma_fit(valid_mask);
all_weights_clean = all_weights(valid_mask);

fprintf('Using %d valid data points for fitting (out of %d total)\n', ...
    length(all_RH_fit_clean), length(all_RH_fit));

% Prepare matrices for linear regression (Y = X*Beta)
% For constrained quadratic: y - 1 = a(x^2 - 10000) + b(x - 100)
Y_vec = all_gamma_fit_clean' - 1;
X_mat = [(all_RH_fit_clean.^2 - 10000)', (all_RH_fit_clean - 100)'];

% Solve using lscov (Weighted Least Squares)
% lscov minimizes (B - A*x)'*diag(w)*(B - A*x)
coeffs = lscov(X_mat, Y_vec, all_weights_clean');

% Extract parameters
a_fit = coeffs(1);
b_fit = coeffs(2);
c_fit = 1 - 10000*a_fit - 100*b_fit;

% Check for valid coefficients
if isnan(a_fit) || isnan(b_fit) || isnan(c_fit) || ...
   isinf(a_fit) || isinf(b_fit) || isinf(c_fit)
    error('Fit failed: coefficients are NaN or Inf');
end

% Generate fit line points
x_fit_line = linspace(0, 100, 200);
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;

%% Goodness of Fit Metrics (R^2 and RMSE)

% Predicted gamma values for the fit data (using clean data)
gamma_pred = a_fit * all_RH_fit_clean.^2 + b_fit * all_RH_fit_clean + c_fit;

% Residuals
residuals = all_gamma_fit_clean - gamma_pred;

% RMSE
RMSE = sqrt(mean(residuals.^2));

% R-squared
SS_res = sum(residuals.^2);
SS_tot = sum((all_gamma_fit - mean(all_gamma_fit)).^2);
R_squared = 1 - SS_res / SS_tot;

% Display fit results
fprintf('\n=== Constrained Fit Results (passes through 100%%, 1) ===\n');
fprintf('Fit equation: gamma_w = %.6f*RH^2 + %.6f*RH + %.6f\n', a_fit, b_fit, c_fit);
fprintf('R^2 = %.4f\n', R_squared);
fprintf('RMSE = %.4f\n', RMSE);
fprintf('Value at 100%% RH: %.6f\n', a_fit*10000 + b_fit*100 + c_fit);

%% Generate distinct colors for all salts
colors = zeros(num_salts, 3);
golden_angle = 0.38196601125; 

for i = 1:num_salts
    hue = mod((i-1) * golden_angle, 1.0);
    sat_cycle = mod(i, 3);
    if sat_cycle == 1, saturation = 0.85 + 0.15 * sin(i * 0.5); 
    elseif sat_cycle == 2, saturation = 0.70 + 0.15 * cos(i * 0.3);
    else, saturation = 0.55 + 0.20 * sin(i * 0.7);
    end
    saturation = max(0.5, min(1.0, saturation));
    
    val_cycle = mod(i, 4);
    if val_cycle == 1, value = 0.90 + 0.10 * cos(i * 0.4);
    elseif val_cycle == 2, value = 0.80 + 0.15 * sin(i * 0.6);
    elseif val_cycle == 3, value = 0.70 + 0.20 * cos(i * 0.5);
    else, value = 0.65 + 0.25 * sin(i * 0.8);
    end
    value = max(0.65, min(1.0, value));
    colors(i,:) = hsv2rgb([hue, saturation, value]);
end

%% FIGURE 1: Activity Coefficient vs Mole Fraction
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.x_water_mol, data.gamma_w_mol, 'LineWidth', 2.5, ...
         'DisplayName', data.display_name, 'color', colors(s,:));
end

plot([0 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Mole Fraction (All Salts, Constrained Fit)', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
xlim([0.3 1.0])
ylim([0 1.2])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save
filename = 'Activity_Coefficient_vs_Mole_Fraction_AllSalts_ConstrainedFit';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, [filename '.fig']))

%% FIGURE 2: Activity Coefficient vs Relative Humidity with Fit
figure('Position', [150, 150, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.gamma_w_mol, 'LineWidth', 2.5, ...
         'DisplayName', data.display_name, 'color', colors(s,:));
end

% Plot the Constrained Fit Line
plot(x_fit_line, y_fit_line, 'k-', 'LineWidth', 3, 'DisplayName', 'Constrained Fit (passes through 100%, 1)')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Relative Humidity (All Salts, Constrained Fit)', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
xlim([0 100])
ylim([0 1.2])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

%% Equation String, Annotate Fit Equation and Statistics
eqn_str = sprintf(['$\\gamma_w = %.6f\\,RH^2 %+ .6f\\,RH %+ .6f$\n' ...
                   '$R^2 = %.4f$,   RMSE = %.4f\n' ...
                   '$\\gamma_w(100\\%%) = %.4f$'], ...
                   a_fit, b_fit, c_fit, R_squared, RMSE, a_fit*10000 + b_fit*100 + c_fit);

annotation('textbox', [0.10 0.70 0.45 0.10], ...
    'String', eqn_str, ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.1);

% Save
filename = 'Activity_Coefficient_vs_RH_AllSalts_ConstrainedFit';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, [filename '.fig']))

disp('Plots generated successfully!')
disp(['  - ' fullfile(fig_out_dir, 'Activity_Coefficient_vs_Mole_Fraction_AllSalts_ConstrainedFit.png')])
disp(['  - ' fullfile(fig_out_dir, 'Activity_Coefficient_vs_RH_AllSalts_ConstrainedFit.png')])