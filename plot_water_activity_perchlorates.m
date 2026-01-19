close all 
clear
clc 

% Add calculate_mf folder to path (ensure your functions are here)
addpath('calculate_mf');

T = 25; 
MWw = 18.015;

%% 1. NaClO4
MW_NaClO4 = 122.44;
% Fit valid: [0.7775, 0.9869]
RH_NaClO4 = linspace(0.778, 0.986, 100); 

for i = 1:length(RH_NaClO4)
    mf_salt_NaClO4(i) = calculate_mf_NaClO4(RH_NaClO4(i));
    mf_water_NaClO4(i) = 1 - mf_salt_NaClO4(i);
    x_water_NaClO4(i) = (mf_water_NaClO4(i) / MWw) / ...
        ((mf_water_NaClO4(i) / MWw) + (mf_salt_NaClO4(i) / MW_NaClO4));
end

%% 2. LiClO4
MW_LiClO4 = 106.39;
% Fit valid: [0.7775, 0.9931]
RH_LiClO4 = linspace(0.778, 0.993, 100);

for i = 1:length(RH_LiClO4)
    mf_salt_LiClO4(i) = calculate_mf_LiClO4(RH_LiClO4(i));
    mf_water_LiClO4(i) = 1 - mf_salt_LiClO4(i);
    x_water_LiClO4(i) = (mf_water_LiClO4(i) / MWw) / ...
        ((mf_water_LiClO4(i) / MWw) + (mf_salt_LiClO4(i) / MW_LiClO4));
end


%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_NaClO4 = RH_NaClO4 ./ x_water_NaClO4;
gamma_LiClO4 = RH_LiClO4 ./ x_water_LiClO4;


%% Calculate Average Perchlorate Fit
% Define the salts to include in the average
fit_salts = {'NaClO4', 'LiClO4'};

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


%% FIGURE 1: Activity Coefficient vs Mole Fraction (Perchlorates)
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

plot(x_water_NaClO4, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4', 'color', [0.75, 0, 0])
plot(x_water_LiClO4, gamma_LiClO4, 'LineWidth', 2.5, 'DisplayName', 'LiClO_4', 'color', [0, 0.5, 0])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Perchlorate Solutions: Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 1)
xlim([0.7 1.0]) % Perchlorates have wider range
ylim([0.9 1.5])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print('figures/Activity_Coefficient_vs_Mole_Fraction_Perchlorates', '-dpng', '-r600')


%% FIGURE 2: Activity Coefficient vs Relative Humidity (Perchlorates)
figure('Position', [150, 150, 900, 700]);
hold on; grid on; box on;

plot(RH_NaClO4*100, gamma_NaClO4, 'LineWidth', 2.5, 'DisplayName', 'NaClO_4', 'color', [0.75, 0, 0])
plot(RH_LiClO4*100, gamma_LiClO4, 'LineWidth', 2.5, 'DisplayName', 'LiClO_4', 'color', [0, 0.5, 0])

% Average Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Perchlorate Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Perchlorate Solutions: Activity Coefficient vs RH', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 1)
xlim([75 100]) % Perchlorates have wider RH range
ylim([0.9 1.5])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print('figures/Activity_Coefficient_vs_RH_Perchlorates', '-dpng', '-r600')

disp('Perchlorate plots generated successfully!')
