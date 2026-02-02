close all 
clear
clc 

% Suppress fzero warnings (they're handled by our error checking)
warning('off', 'MATLAB:fzero:AbortFcn');

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Define all salts with their molecular weights, valid RH ranges, and ION COUNTS
% Format: {salt_name, MW, RH_min, RH_max, function_name, function_args, num_cations, num_anions}
% function_args: 0 = no args, 1 = needs T

salt_data = {
    % Endothermic salts
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1};
    {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1};
    {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1};
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};
    
    % Exothermic salts
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2};
    {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
    {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1};
    {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2};
    {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
    {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1};
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1};
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2};
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2};
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2};
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2};
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2};
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2};
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2};
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
};

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
    n_cat = salt_data{s}{7};
    n_an  = salt_data{s}{8};
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