close all 
clear
clc 

% Add calculate_mf, util, and data folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', '..', 'util'));
addpath(fullfile(filepath, '..', '..', 'data'));

T = 25; 
MWw = 18.015;

%% Load salt data and compute RH, mole fraction, gamma for each exothermic salt
salt_data = load_salt_data();
% All salts for plots (includes MgNO32); fit excludes MgNO32
plot_salts = {'LiCl', 'CaCl2', 'MgCl2', 'LiBr', 'ZnCl2', 'LiI', 'ZnBr2', 'ZnI2', ...
              'HCl', 'MgNO32', 'LiOH', 'NaOH'};
fit_salts = {'LiCl', 'CaCl2', 'MgCl2', 'LiBr', 'ZnCl2', ...
             'LiI', 'ZnBr2', 'ZnI2', 'HCl', 'LiOH', 'NaOH'};  % Exclude MgNO32 from fit
num_points = 100;

for k = 1:length(plot_salts)
    salt_name = plot_salts{k};
    idx = find(cellfun(@(r) strcmp(r{1}, salt_name), salt_data), 1);
    if isempty(idx), error('Salt %s not found in load_salt_data', salt_name); end
    r = salt_data{idx};
    MW_salt = r{2}; RH_min = r{3}; RH_max = r{4}; func_name = r{5}; func_args = r{6}; T_salt = r{9};
    
    RH_vec = linspace(RH_min, RH_max, num_points);
    for i = 1:num_points
        if func_args == 1
            mf_s = feval(func_name, RH_vec(i), T_salt);
        else
            mf_s = feval(func_name, RH_vec(i));
        end
        mf_w = 1 - mf_s;
        x_w = (mf_w / MWw) / ((mf_w / MWw) + (mf_s / MW_salt));
        if i == 1, mf_salt_vec = zeros(size(RH_vec)); mf_water_vec = mf_salt_vec; x_water_vec = mf_salt_vec; end
        mf_salt_vec(i) = mf_s; mf_water_vec(i) = mf_w; x_water_vec(i) = x_w;
    end
    gamma_vec = RH_vec ./ x_water_vec;
    eval(['RH_' salt_name ' = RH_vec;']);
    eval(['x_water_' salt_name ' = x_water_vec;']);
    eval(['gamma_' salt_name ' = gamma_vec;']);
end

%% Calculate Constrained Weighted Polynomial Fit (Exclude MgNO32)

% Initialize container arrays for fitting
all_RH_fit = [];
all_gamma_fit = [];
all_weights = [];

% Loop through salts to aggregate data and calculate weights
for k = 1:length(fit_salts)
    salt_name = fit_salts{k};
    
    % Dynamically fetch the RH and Gamma vectors from the workspace
    % (Assumes variable naming convention RH_Name and gamma_Name)
    rh_vec = eval(['RH_' salt_name]);
    gamma_vec = eval(['gamma_' salt_name]);
    
    % Calculate the Range of RH for this salt to use as weight
    rh_range = max(rh_vec) - min(rh_vec);
    
    % Assign weight (range) to every point in this dataset
    w_vec = ones(size(rh_vec)) * rh_range;
    
    % Append to master lists (Convert RH to %)
    all_RH_fit = [all_RH_fit, rh_vec * 100]; 
    all_gamma_fit = [all_gamma_fit, gamma_vec];
    all_weights = [all_weights, w_vec];
end

% --- Perform Constrained Least Squares ---
% Quadratic Model: y = a*x^2 + b*x + c
% Constraint: Pass through (100, 1) -> 1 = a(10000) + b(100) + c
% Solve for c: c = 1 - 10000a - 100b
% Substitute c back into model: 
% y = a*x^2 + b*x + (1 - 10000a - 100b)
% Rearrange to isolate a and b:
% y - 1 = a(x^2 - 10000) + b(x - 100)

% Prepare matrices for linear regression (Y = X*Beta)
Y_vec = all_gamma_fit' - 1;
X_mat = [(all_RH_fit.^2 - 10000)', (all_RH_fit - 100)'];

% Solve using lscov (Weighted Least Squares)
% lscov minimizes (B - A*x)'*diag(w)*(B - A*x)
coeffs = lscov(X_mat, Y_vec, all_weights');

% Extract parameters
a_fit = coeffs(1);
b_fit = coeffs(2);
c_fit = 1 - 10000*a_fit - 100*b_fit;

% Generate fit line points
x_fit_line = linspace(0, 100, 200);
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;

%% Goodness of Fit Metrics (R^2 and RMSE)

% Predicted gamma values for the fit data
gamma_pred = a_fit * all_RH_fit.^2 + b_fit * all_RH_fit + c_fit;

% Residuals
residuals = all_gamma_fit - gamma_pred;

% RMSE
RMSE = sqrt(mean(residuals.^2));

% R-squared
SS_res = sum(residuals.^2);
SS_tot = sum((all_gamma_fit - mean(all_gamma_fit)).^2);
R_squared = 1 - SS_res / SS_tot;


%% FIGURE 1: Activity Coefficient vs Mole Fraction
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;
plot(x_water_LiCl, gamma_LiCl, 'LineWidth', 2.5, 'DisplayName', 'LiCl', 'color', [0, 0.5, 0])
plot(x_water_CaCl2, gamma_CaCl2, 'LineWidth', 2.5, 'DisplayName', 'CaCl_2', 'color', [0.9290 0.6940 0.1250])
plot(x_water_MgCl2, gamma_MgCl2, 'LineWidth', 2.5, 'DisplayName', 'MgCl_2', 'color', [0.8500 0.3250 0.0980])
plot(x_water_LiBr, gamma_LiBr, 'LineWidth', 2.5, 'DisplayName', 'LiBr', 'color', [0.3010, 0.7450, 0.9330])
plot(x_water_ZnCl2, gamma_ZnCl2, 'LineWidth', 2.5, 'DisplayName', 'ZnCl_2', 'color', [0.6350 0.0780 0.1840])
plot(x_water_LiI, gamma_LiI, 'LineWidth', 2.5, 'DisplayName', 'LiI', 'color', [0.4940 0.1840 0.5560])
plot(x_water_ZnBr2, gamma_ZnBr2, 'LineWidth', 2.5, 'DisplayName', 'ZnBr_2', 'color', [0.4660 0.6740 0.1880])
plot(x_water_ZnI2, gamma_ZnI2, 'LineWidth', 2.5, 'DisplayName', 'ZnI_2', 'LineStyle', '--', 'color', [0.9290 0.6940 0.1250])
plot(x_water_HCl, gamma_HCl, 'LineWidth', 2.5, 'DisplayName', 'HCl', 'color', [0, 0.4470, 0.7410])
plot(x_water_MgNO32, gamma_MgNO32, 'LineWidth', 2.5, 'DisplayName', 'Mg(NO_3)_2', 'color', [0 1 1])
plot(x_water_LiOH, gamma_LiOH, 'LineWidth', 2.5, 'DisplayName', 'LiOH', 'color', [1 0 1])
plot(x_water_NaOH, gamma_NaOH, 'LineWidth', 2.5, 'DisplayName', 'NaOH', 'color', [0.75 0 0.75])
plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northeast', 'FontSize', 11)
xlim([0.5 1.0])
ylim([0 1.2])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Exothermic_fitted'), '-dpng', '-r600')
savefig(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Exothermic_fitted.fig'))


%% FIGURE 2: Activity Coefficient vs Relative Humidity
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;
plot(RH_LiCl*100, gamma_LiCl, 'LineWidth', 2.5, 'DisplayName', 'LiCl', 'color', [0, 0.5, 0])
plot(RH_CaCl2*100, gamma_CaCl2, 'LineWidth', 2.5, 'DisplayName', 'CaCl_2', 'color', [0.9290 0.6940 0.1250])
plot(RH_MgCl2*100, gamma_MgCl2, 'LineWidth', 2.5, 'DisplayName', 'MgCl_2', 'color', [0.8500 0.3250 0.0980])
plot(RH_LiBr*100, gamma_LiBr, 'LineWidth', 2.5, 'DisplayName', 'LiBr', 'color', [0.3010, 0.7450, 0.9330])
plot(RH_ZnCl2*100, gamma_ZnCl2, 'LineWidth', 2.5, 'DisplayName', 'ZnCl_2', 'color', [0.6350 0.0780 0.1840])
plot(RH_LiI*100, gamma_LiI, 'LineWidth', 2.5, 'DisplayName', 'LiI', 'color', [0.4940 0.1840 0.5560])
plot(RH_ZnBr2*100, gamma_ZnBr2, 'LineWidth', 2.5, 'DisplayName', 'ZnBr_2', 'color', [0.4660 0.6740 0.1880])
plot(RH_ZnI2*100, gamma_ZnI2, 'LineWidth', 2.5, 'DisplayName', 'ZnI_2', 'LineStyle', '--', 'color', [0.9290 0.6940 0.1250])
plot(RH_HCl*100, gamma_HCl, 'LineWidth', 2.5, 'DisplayName', 'HCl', 'color', [0, 0.4470, 0.7410])
plot(RH_MgNO32*100, gamma_MgNO32, 'LineWidth', 2.5, 'DisplayName', 'Mg(NO_3)_2', 'color', [0 1 1])
plot(RH_LiOH*100, gamma_LiOH, 'LineWidth', 2.5, 'DisplayName', 'LiOH', 'color', [1 0 1])
plot(RH_NaOH*100, gamma_NaOH, 'LineWidth', 2.5, 'DisplayName', 'NaOH', 'color', [0.75 0 0.75])

% Plot the Average Fit Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Relative Humidity', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northeast', 'FontSize', 11)
xlim([0 100])
ylim([0 1.2])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 

%% Equation String, Annotate Fit Equation and Statistics

eqn_str = sprintf(['$\\gamma_w = %.6f\\,RH^2 %+ .6f\\,RH %+ .6f$\n' ...
                   '$R^2 = %.4f$,   RMSE = %.4f'], ...
                   a_fit, b_fit, c_fit, R_squared, RMSE);

annotation('textbox', [0.10 0.70 0.45 0.08], ...
    'String', eqn_str, ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.1);

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Exothermic_fitted'), '-dpng', '-r600')
savefig(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Exothermic_fitted.fig'))

disp('Plots generated successfully!')
disp('  - figures/activity_coefficient/activity_coefficient_subsets/Activity_Coefficient_vs_Mole_Fraction_Exothermic_fitted.png')
disp('  - figures/activity_coefficient/activity_coefficient_subsets/Activity_Coefficient_vs_RH_Exothermic_fitted.png')