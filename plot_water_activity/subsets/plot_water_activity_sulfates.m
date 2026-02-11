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

%% Load salt data and compute RH, mole fraction, gamma for each sulfate
salt_data = load_salt_data();
fit_salts = {'Na2SO4', 'K2SO4', 'NH42SO4', 'MgSO4', 'MnSO4', ...
             'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4'};
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


%% Calculate Average Sulfate Fit
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

x_fit_line = linspace(70, 100, 200); % Plot range
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;


%% FIGURE 1: Activity Coefficient vs Mole Fraction (Sulfates)
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

plot(x_water_Na2SO4,  gamma_Na2SO4,  'LineWidth', 2.5, 'DisplayName', 'Na_2SO_4', 'color', [0.75, 0, 0])
plot(x_water_K2SO4,   gamma_K2SO4,   'LineWidth', 2.5, 'DisplayName', 'K_2SO_4',  'color', [0, 0.75, 0.75])
plot(x_water_NH42SO4, gamma_NH42SO4, 'LineWidth', 2.5, 'DisplayName', '(NH_4)_2SO_4', 'color', [0.5, 0, 0.5])
plot(x_water_MgSO4,   gamma_MgSO4,   'LineWidth', 2.5, 'DisplayName', 'MgSO_4',   'color', [0.6, 0.6, 0.6])
plot(x_water_MnSO4,   gamma_MnSO4,   'LineWidth', 2.5, 'DisplayName', 'MnSO_4',   'color', [0.8, 0.4, 0])
plot(x_water_Li2SO4,  gamma_Li2SO4,  'LineWidth', 2.5, 'DisplayName', 'Li_2SO_4', 'color', [0, 0.5, 0])
plot(x_water_NiSO4,   gamma_NiSO4,   'LineWidth', 2.5, 'DisplayName', 'NiSO_4',   'color', [0.2, 0.8, 0.2])
plot(x_water_CuSO4,   gamma_CuSO4,   'LineWidth', 2.5, 'DisplayName', 'CuSO_4',   'color', [0, 0, 1])
plot(x_water_ZnSO4,   gamma_ZnSO4,   'LineWidth', 2.5, 'DisplayName', 'ZnSO_4',   'color', [0.5, 0.2, 0.2])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Sulfate Solutions: Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([0.9 1.0]) % Sulfates are mostly dilute, so zoom in near 1
ylim([0.7 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Sulfates'), '-dpng', '-r600')


%% FIGURE 2: Activity Coefficient vs Relative Humidity (Sulfates)
figure('Position', [150, 150, 900, 700]);
hold on; grid on; box on;

plot(RH_Na2SO4*100,  gamma_Na2SO4,  'LineWidth', 2.5, 'DisplayName', 'Na_2SO_4', 'color', [0.75, 0, 0])
plot(RH_K2SO4*100,   gamma_K2SO4,   'LineWidth', 2.5, 'DisplayName', 'K_2SO_4',  'color', [0, 0.75, 0.75])
plot(RH_NH42SO4*100, gamma_NH42SO4, 'LineWidth', 2.5, 'DisplayName', '(NH_4)_2SO_4', 'color', [0.5, 0, 0.5])
plot(RH_MgSO4*100,   gamma_MgSO4,   'LineWidth', 2.5, 'DisplayName', 'MgSO_4',   'color', [0.6, 0.6, 0.6])
plot(RH_MnSO4*100,   gamma_MnSO4,   'LineWidth', 2.5, 'DisplayName', 'MnSO_4',   'color', [0.8, 0.4, 0])
plot(RH_Li2SO4*100,  gamma_Li2SO4,  'LineWidth', 2.5, 'DisplayName', 'Li_2SO_4', 'color', [0, 0.5, 0])
plot(RH_NiSO4*100,   gamma_NiSO4,   'LineWidth', 2.5, 'DisplayName', 'NiSO_4',   'color', [0.2, 0.8, 0.2])
plot(RH_CuSO4*100,   gamma_CuSO4,   'LineWidth', 2.5, 'DisplayName', 'CuSO_4',   'color', [0, 0, 1])
plot(RH_ZnSO4*100,   gamma_ZnSO4,   'LineWidth', 2.5, 'DisplayName', 'ZnSO_4',   'color', [0.5, 0.2, 0.2])

% Average Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Sulfate Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Sulfate Solutions: Activity Coefficient vs RH', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([80 100]) % Sulfates data is mostly high RH
ylim([0.7 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Sulfates'), '-dpng', '-r600')

disp('Sulfate plots generated successfully!')