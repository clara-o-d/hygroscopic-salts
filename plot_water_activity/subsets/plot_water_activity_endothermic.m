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

%% Load salt data and compute RH, mole fraction, gamma for each endothermic salt
salt_data = load_salt_data();
fit_salts = {'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', 'LiNO3', 'KNO3', ...
             'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI'};
num_points = 100;

for k = 1:length(fit_salts)
    salt_name = fit_salts{k};
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

%% Calculate Constrained Weighted Polynomial Fit
% Initialize container arrays for fitting
all_RH_fit = [];
all_gamma_fit = [];
all_weights = [];

% Loop through salts to aggregate data and calculate weights
for k = 1:length(fit_salts)
    salt_name = fit_salts{k};
    
    % Dynamically fetch the RH and Gamma vectors from the workspace
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

% Generate distinct colors for all salts
colors = lines(17);

%% FIGURE 1: Activity Coefficient vs Mole Fraction
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% plot(x_water_NaCl,  gamma_NaCl,  'LineWidth', 2.5, 'DisplayName', 'NaCl',  'color', colors(1,:))
plot(x_water_KCl,   gamma_KCl,   'LineWidth', 2.5, 'DisplayName', 'KCl',   'color', colors(2,:))
plot(x_water_NH4Cl, gamma_NH4Cl, 'LineWidth', 2.5, 'DisplayName', 'NH_4Cl','color', colors(3,:))
plot(x_water_CsCl,  gamma_CsCl,  'LineWidth', 2.5, 'DisplayName', 'CsCl',  'color', colors(4,:))
plot(x_water_NaNO3, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3','color', colors(5,:))
plot(x_water_AgNO3, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3','color', colors(6,:))
plot(x_water_KI,    gamma_KI,    'LineWidth', 2.5, 'DisplayName', 'KI',    'color', colors(7,:))
plot(x_water_LiNO3, gamma_LiNO3, 'LineWidth', 2.5, 'DisplayName', 'LiNO_3','color', colors(8,:))
plot(x_water_KNO3,  gamma_KNO3,  'LineWidth', 2.5, 'DisplayName', 'KNO_3', 'color', colors(9,:))
plot(x_water_NaClO4, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4','color', colors(10,:))
plot(x_water_KClO3, gamma_KClO3, 'LineWidth', 2.5, 'DisplayName', 'KClO_3','color', colors(11,:))
plot(x_water_NaBr,  gamma_NaBr,  'LineWidth', 2.5, 'DisplayName', 'NaBr',  'color', colors(12,:))
plot(x_water_NaI,   gamma_NaI,   'LineWidth', 2.5, 'DisplayName', 'NaI',   'color', colors(13,:))
plot(x_water_KBr,   gamma_KBr,   'LineWidth', 2.5, 'DisplayName', 'KBr',   'color', colors(14,:))
plot(x_water_RbCl,  gamma_RbCl,  'LineWidth', 2.5, 'DisplayName', 'RbCl',  'color', colors(15,:))
plot(x_water_CsBr,  gamma_CsBr,  'LineWidth', 2.5, 'DisplayName', 'CsBr',  'color', colors(16,:))
plot(x_water_CsI,   gamma_CsI,   'LineWidth', 2.5, 'DisplayName', 'CsI',   'color', colors(17,:))

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10)
xlim([0.2 1.0]) % Adjusted xlim since these salts are mostly soluble/high aw
ylim([0.6 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Endothermic'), '-dpng', '-r600')

%% FIGURE 2: Activity Coefficient vs Relative Humidity
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% Use the same colors as Figure 1
% plot(RH_NaCl*100,  gamma_NaCl,  'LineWidth', 2.5, 'DisplayName', 'NaCl',  'color', colors(1,:))
plot(RH_KCl*100,   gamma_KCl,   'LineWidth', 2.5, 'DisplayName', 'KCl',   'color', colors(2,:))
plot(RH_NH4Cl*100, gamma_NH4Cl, 'LineWidth', 2.5, 'DisplayName', 'NH_4Cl','color', colors(3,:))
plot(RH_CsCl*100,  gamma_CsCl,  'LineWidth', 2.5, 'DisplayName', 'CsCl',  'color', colors(4,:))
plot(RH_NaNO3*100, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3','color', colors(5,:))
plot(RH_AgNO3*100, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3','color', colors(6,:))
plot(RH_KI*100,    gamma_KI,    'LineWidth', 2.5, 'DisplayName', 'KI',    'color', colors(7,:))
plot(RH_LiNO3*100, gamma_LiNO3, 'LineWidth', 2.5, 'DisplayName', 'LiNO_3','color', colors(8,:))
plot(RH_KNO3*100,  gamma_KNO3,  'LineWidth', 2.5, 'DisplayName', 'KNO_3', 'color', colors(9,:))
plot(RH_NaClO4*100, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4','color', colors(10,:))
plot(RH_KClO3*100, gamma_KClO3, 'LineWidth', 2.5, 'DisplayName', 'KClO_3','color', colors(11,:))
plot(RH_NaBr*100,  gamma_NaBr,  'LineWidth', 2.5, 'DisplayName', 'NaBr',  'color', colors(12,:))
plot(RH_NaI*100,   gamma_NaI,   'LineWidth', 2.5, 'DisplayName', 'NaI',   'color', colors(13,:))
plot(RH_KBr*100,   gamma_KBr,   'LineWidth', 2.5, 'DisplayName', 'KBr',   'color', colors(14,:))
plot(RH_RbCl*100,  gamma_RbCl,  'LineWidth', 2.5, 'DisplayName', 'RbCl',  'color', colors(15,:))
plot(RH_CsBr*100,  gamma_CsBr,  'LineWidth', 2.5, 'DisplayName', 'CsBr',  'color', colors(16,:))
plot(RH_CsI*100,   gamma_CsI,   'LineWidth', 2.5, 'DisplayName', 'CsI',   'color', colors(17,:))

% Plot the Average Fit Line
% plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Relative Humidity', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10)
xlim([50 100]) % Adjusted xlim to include all data ranges
ylim([0.6 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Endothermic'), '-dpng', '-r600')

disp('Plots generated successfully!')
disp('  - figures/activity_coefficient/activity_coefficient_subsets/Activity_Coefficient_vs_Mole_Fraction_Endothermic.png')
disp('  - figures/activity_coefficient/activity_coefficient_subsets/Activity_Coefficient_vs_RH_Endothermic.png')