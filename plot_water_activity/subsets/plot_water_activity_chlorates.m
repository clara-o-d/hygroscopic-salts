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

%% Load salt data and compute RH, mole fraction, gamma for each chlorate
salt_data = load_salt_data();
fit_salts = {'NaClO4', 'LiClO4', 'KClO3'};
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


%% Calculate Average Chlorate Fit
% Define the salts to include in the average
fit_salts = {'NaClO4', 'LiClO4', 'KClO3'};

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

x_fit_line = linspace(75, 100, 200); % Plot range
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;


%% FIGURE 1: Activity Coefficient vs Mole Fraction (Chlorates)
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

plot(x_water_NaClO4, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4', 'color', [0.75, 0, 0])
plot(x_water_LiClO4, gamma_LiClO4, 'LineWidth', 2.5, 'DisplayName', 'LiClO_4', 'color', [0, 0.5, 0])
plot(x_water_KClO3,  gamma_KClO3,  'LineWidth', 2.5, 'DisplayName', 'KClO_3',  'color', [0, 0, 0.75])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Chlorate Solutions: Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 1)
xlim([0.7 1.0]) % Chlorates have wider range
ylim([0.8 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Chlorates'), '-dpng', '-r600')

%% FIGURE 2: Activity Coefficient vs Relative Humidity (Chlorates)
figure('Position', [150, 150, 900, 700]);
hold on; grid on; box on;

plot(RH_NaClO4*100, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4', 'color', [0.75, 0, 0])
plot(RH_LiClO4*100, gamma_LiClO4, 'LineWidth', 2.5, 'DisplayName', 'LiClO_4', 'color', [0, 0.5, 0])
plot(RH_KClO3*100,  gamma_KClO3,  'LineWidth', 2.5, 'DisplayName', 'KClO_3',  'color', [0, 0, 0.75])

% Average Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Chlorate Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Chlorate Solutions: Activity Coefficient vs RH', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 1)
xlim([75 100]) % Chlorates have wider RH range
ylim([0.8 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Chlorates'), '-dpng', '-r600')

disp('Chlorate plots generated successfully!')
