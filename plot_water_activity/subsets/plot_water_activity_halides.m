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

%% Load salt data and compute RH, mole fraction, gamma for each halide
salt_data = load_salt_data();
fit_salts = {'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI', 'CaBr2', 'CaI2', ...
             'SrCl2', 'SrBr2', 'SrI2', 'BaCl2', 'BaBr2'};
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


%% Calculate Average Halide Fit
all_RH_fit = [];
all_gamma_fit = [];
all_weights = [];

for k = 1:length(fit_salts)
    salt_name = fit_salts{k};
    rh_vec = eval(['RH_' salt_name]);
    gamma_vec = eval(['gamma_' salt_name]);
    
    % Weighting by range of data coverage
    rh_range = max(rh_vec) - min(rh_vec);
    w_vec = ones(size(rh_vec)) * rh_range;
    
    all_RH_fit = [all_RH_fit, rh_vec * 100]; 
    all_gamma_fit = [all_gamma_fit, gamma_vec];
    all_weights = [all_weights, w_vec];
end

% Constrained Fit: Pass through (100, 1)
Y_vec = all_gamma_fit' - 1;
X_mat = [(all_RH_fit.^2 - 10000)', (all_RH_fit - 100)'];
coeffs = lscov(X_mat, Y_vec, all_weights');

a_fit = coeffs(1);
b_fit = coeffs(2);
c_fit = 1 - 10000*a_fit - 100*b_fit;

x_fit_line = linspace(55, 100, 200); % Plot range
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;


%% FIGURE 1: Activity Coefficient vs Mole Fraction (Halides)
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

plot(x_water_NaBr, gamma_NaBr, 'LineWidth', 2.5, 'DisplayName', 'NaBr', 'color', [0.75, 0, 0])
plot(x_water_NaI, gamma_NaI, 'LineWidth', 2.5, 'DisplayName', 'NaI', 'color', [0.5, 0.5, 0.5])
plot(x_water_KBr, gamma_KBr, 'LineWidth', 2.5, 'DisplayName', 'KBr', 'color', [0, 0.75, 0.75])
plot(x_water_RbCl, gamma_RbCl, 'LineWidth', 2.5, 'DisplayName', 'RbCl', 'color', [0.5, 0, 0.5])
plot(x_water_CsBr, gamma_CsBr, 'LineWidth', 2.5, 'DisplayName', 'CsBr', 'color', [0.8, 0.4, 0])
plot(x_water_CsI, gamma_CsI, 'LineWidth', 2.5, 'DisplayName', 'CsI', 'color', [0.2, 0.8, 0.2])
plot(x_water_CaBr2, gamma_CaBr2, 'LineWidth', 2.5, 'DisplayName', 'CaBr_2', 'color', [0, 0, 1])
plot(x_water_CaI2, gamma_CaI2, 'LineWidth', 2.5, 'DisplayName', 'CaI_2', 'color', [1, 0, 1])
plot(x_water_SrCl2, gamma_SrCl2, 'LineWidth', 2.5, 'DisplayName', 'SrCl_2', 'color', [0.6, 0.3, 0.1])
plot(x_water_SrBr2, gamma_SrBr2, 'LineWidth', 2.5, 'DisplayName', 'SrBr_2', 'color', [0.3, 0.7, 0.5])
plot(x_water_SrI2, gamma_SrI2, 'LineWidth', 2.5, 'DisplayName', 'SrI_2', 'color', [0.7, 0.2, 0.8])
plot(x_water_BaCl2, gamma_BaCl2, 'LineWidth', 2.5, 'DisplayName', 'BaCl_2', 'color', [0.2, 0.6, 0.8])
plot(x_water_BaBr2, gamma_BaBr2, 'LineWidth', 2.5, 'DisplayName', 'BaBr_2', 'color', [0.9, 0.6, 0.2])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Halide Solutions: Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([0.6 1.0]) % Halides have wide range
ylim([0.6 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Halides'), '-dpng', '-r600')


%% FIGURE 2: Activity Coefficient vs Relative Humidity (Halides)
figure('Position', [150, 150, 900, 700]);
hold on; grid on; box on;

plot(RH_NaBr*100, gamma_NaBr, 'LineWidth', 2.5, 'DisplayName', 'NaBr', 'color', [0.75, 0, 0])
plot(RH_NaI*100, gamma_NaI, 'LineWidth', 2.5, 'DisplayName', 'NaI', 'color', [0.5, 0.5, 0.5])
plot(RH_KBr*100, gamma_KBr, 'LineWidth', 2.5, 'DisplayName', 'KBr', 'color', [0, 0.75, 0.75])
plot(RH_RbCl*100, gamma_RbCl, 'LineWidth', 2.5, 'DisplayName', 'RbCl', 'color', [0.5, 0, 0.5])
plot(RH_CsBr*100, gamma_CsBr, 'LineWidth', 2.5, 'DisplayName', 'CsBr', 'color', [0.8, 0.4, 0])
plot(RH_CsI*100, gamma_CsI, 'LineWidth', 2.5, 'DisplayName', 'CsI', 'color', [0.2, 0.8, 0.2])
plot(RH_CaBr2*100, gamma_CaBr2, 'LineWidth', 2.5, 'DisplayName', 'CaBr_2', 'color', [0, 0, 1])
plot(RH_CaI2*100, gamma_CaI2, 'LineWidth', 2.5, 'DisplayName', 'CaI_2', 'color', [1, 0, 1])
plot(RH_SrCl2*100, gamma_SrCl2, 'LineWidth', 2.5, 'DisplayName', 'SrCl_2', 'color', [0.6, 0.3, 0.1])
plot(RH_SrBr2*100, gamma_SrBr2, 'LineWidth', 2.5, 'DisplayName', 'SrBr_2', 'color', [0.3, 0.7, 0.5])
plot(RH_SrI2*100, gamma_SrI2, 'LineWidth', 2.5, 'DisplayName', 'SrI_2', 'color', [0.7, 0.2, 0.8])
plot(RH_BaCl2*100, gamma_BaCl2, 'LineWidth', 2.5, 'DisplayName', 'BaCl_2', 'color', [0.2, 0.6, 0.8])
plot(RH_BaBr2*100, gamma_BaBr2, 'LineWidth', 2.5, 'DisplayName', 'BaBr_2', 'color', [0.9, 0.6, 0.2])

% Average Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Halide Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Halide Solutions: Activity Coefficient vs RH', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([55 100]) % Halides have wide RH range
ylim([0.6 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Halides'), '-dpng', '-r600')

disp('Halide plots generated successfully!')
