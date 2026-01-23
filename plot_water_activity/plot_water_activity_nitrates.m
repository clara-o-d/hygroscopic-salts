close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T = 25; 
MWw = 18.015;

%% 1. NaNO3
MW_NaNO3 = 85.00;
% Fit valid: [0.9701, 0.9996] (from calculate_mf_NaNO3_)
RH_NaNO3 = linspace(0.9701, 0.999, 100); 

for i = 1:length(RH_NaNO3)
    mf_salt_NaNO3(i) = calculate_mf_NaNO3(RH_NaNO3(i));
    mf_water_NaNO3(i) = 1 - mf_salt_NaNO3(i);
    x_water_NaNO3(i) = (mf_water_NaNO3(i) / MWw) / ...
        ((mf_water_NaNO3(i) / MWw) + (mf_salt_NaNO3(i) / MW_NaNO3));
end

%% 2. AgNO3
MW_AgNO3 = 169.87;
% Fit valid: approximately [0.896, 0.986] (from data range)
RH_AgNO3 = linspace(0.897, 0.985, 100);

for i = 1:length(RH_AgNO3)
    mf_salt_AgNO3(i) = calculate_mf_AgNO3(RH_AgNO3(i));
    mf_water_AgNO3(i) = 1 - mf_salt_AgNO3(i);
    x_water_AgNO3(i) = (mf_water_AgNO3(i) / MWw) / ...
        ((mf_water_AgNO3(i) / MWw) + (mf_salt_AgNO3(i) / MW_AgNO3));
end

%% 3. LiNO3
MW_LiNO3 = 68.95;
% Fit valid: [0.7353, 0.9967]
RH_LiNO3 = linspace(0.736, 0.995, 100);

for i = 1:length(RH_LiNO3)
    mf_salt_LiNO3(i) = calculate_mf_LiNO3(RH_LiNO3(i));
    mf_water_LiNO3(i) = 1 - mf_salt_LiNO3(i);
    x_water_LiNO3(i) = (mf_water_LiNO3(i) / MWw) / ...
        ((mf_water_LiNO3(i) / MWw) + (mf_salt_LiNO3(i) / MW_LiNO3));
end

%% 3b. NH4NO3
MW_NH4NO3 = 80.043;
% Fit valid: [0.118, 0.732]
RH_NH4NO3 = linspace(0.118, 0.732, 100);

for i = 1:length(RH_NH4NO3)
    mf_salt_NH4NO3(i) = calculate_mf_NH4NO3(RH_NH4NO3(i));
    mf_water_NH4NO3(i) = 1 - mf_salt_NH4NO3(i);
    x_water_NH4NO3(i) = (mf_water_NH4NO3(i) / MWw) / ...
        ((mf_water_NH4NO3(i) / MWw) + (mf_salt_NH4NO3(i) / MW_NH4NO3));
end

%% 4. KNO3
MW_KNO3 = 101.10;
% Fit valid: [0.9315, 0.9967]
RH_KNO3 = linspace(0.932, 0.995, 100);

for i = 1:length(RH_KNO3)
    mf_salt_KNO3(i) = calculate_mf_KNO3(RH_KNO3(i));
    mf_water_KNO3(i) = 1 - mf_salt_KNO3(i);
    x_water_KNO3(i) = (mf_water_KNO3(i) / MWw) / ...
        ((mf_water_KNO3(i) / MWw) + (mf_salt_KNO3(i) / MW_KNO3));
end

%% 5. Ba(NO3)2
MW_BaNO3 = 261.34;
% Fit valid: [0.9859, 0.9958]
RH_BaNO3 = linspace(0.986, 0.995, 100);

for i = 1:length(RH_BaNO3)
    mf_salt_BaNO3(i) = calculate_mf_BaNO32(RH_BaNO3(i));
    mf_water_BaNO3(i) = 1 - mf_salt_BaNO3(i);
    x_water_BaNO3(i) = (mf_water_BaNO3(i) / MWw) / ...
        ((mf_water_BaNO3(i) / MWw) + (mf_salt_BaNO3(i) / MW_BaNO3));
end

%% 6. Ca(NO3)2
MW_CaNO3 = 164.09;
% Fit valid: [0.6464, 0.9955]
RH_CaNO3 = linspace(0.647, 0.995, 100);

for i = 1:length(RH_CaNO3)
    mf_salt_CaNO3(i) = calculate_mf_CaNO32(RH_CaNO3(i));
    mf_water_CaNO3(i) = 1 - mf_salt_CaNO3(i);
    x_water_CaNO3(i) = (mf_water_CaNO3(i) / MWw) / ...
        ((mf_water_CaNO3(i) / MWw) + (mf_salt_CaNO3(i) / MW_CaNO3));
end


%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_NaNO3 = RH_NaNO3 ./ x_water_NaNO3;
gamma_AgNO3 = RH_AgNO3 ./ x_water_AgNO3;
gamma_LiNO3 = RH_LiNO3 ./ x_water_LiNO3;
gamma_NH4NO3 = RH_NH4NO3 ./ x_water_NH4NO3;
gamma_KNO3  = RH_KNO3  ./ x_water_KNO3;
gamma_BaNO3 = RH_BaNO3 ./ x_water_BaNO3;
gamma_CaNO3 = RH_CaNO3 ./ x_water_CaNO3;


%% Calculate Average Nitrate Fit
% Define the salts to include in the average
fit_salts = {'NaNO3', 'AgNO3', 'LiNO3', 'NH4NO3', 'KNO3', 'BaNO3', 'CaNO3'};

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

x_fit_line = linspace(60, 100, 200); % Plot range
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;


%% FIGURE 1: Activity Coefficient vs Mole Fraction (Nitrates)
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

plot(x_water_NaNO3, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3', 'color', [0.75, 0, 0])
plot(x_water_AgNO3, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3', 'color', [0.5, 0.5, 0.5])
plot(x_water_LiNO3, gamma_LiNO3, 'LineWidth', 2.5, 'DisplayName', 'LiNO_3', 'color', [0, 0.5, 0])
plot(x_water_NH4NO3, gamma_NH4NO3, 'LineWidth', 2.5, 'DisplayName', 'NH_4NO_3', 'color', [0.75, 0.5, 0])
plot(x_water_KNO3,  gamma_KNO3,  'LineWidth', 2.5, 'DisplayName', 'KNO_3',  'color', [0, 0.75, 0.75])
plot(x_water_BaNO3, gamma_BaNO3, 'LineWidth', 2.5, 'DisplayName', 'Ba(NO_3)_2', 'color', [0.5, 0, 0.5])
plot(x_water_CaNO3, gamma_CaNO3, 'LineWidth', 2.5, 'DisplayName', 'Ca(NO_3)_2', 'color', [0.8, 0.4, 0])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Nitrate Solutions: Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([0.7 1.0]) % Nitrates have wider range than sulfates
ylim([0.7 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Nitrates'), '-dpng', '-r600')


%% FIGURE 2: Activity Coefficient vs Relative Humidity (Nitrates)
figure('Position', [150, 150, 900, 700]);
hold on; grid on; box on;

plot(RH_NaNO3*100, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3', 'color', [0.75, 0, 0])
plot(RH_AgNO3*100, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3', 'color', [0.5, 0.5, 0.5])
plot(RH_LiNO3*100, gamma_LiNO3, 'LineWidth', 2.5, 'DisplayName', 'LiNO_3', 'color', [0, 0.5, 0])
plot(RH_NH4NO3*100, gamma_NH4NO3, 'LineWidth', 2.5, 'DisplayName', 'NH_4NO_3', 'color', [0.75, 0.5, 0])
plot(RH_KNO3*100,  gamma_KNO3,  'LineWidth', 2.5, 'DisplayName', 'KNO_3',  'color', [0, 0.75, 0.75])
plot(RH_BaNO3*100, gamma_BaNO3, 'LineWidth', 2.5, 'DisplayName', 'Ba(NO_3)_2', 'color', [0.5, 0, 0.5])
plot(RH_CaNO3*100, gamma_CaNO3, 'LineWidth', 2.5, 'DisplayName', 'Ca(NO_3)_2', 'color', [0.8, 0.4, 0])

% Average Line
plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Nitrate Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Nitrate Solutions: Activity Coefficient vs RH', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', 2)
xlim([60 100]) % Nitrates have wider RH range
ylim([0.7 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Nitrates'), '-dpng', '-r600')

disp('Nitrate plots generated successfully!')
