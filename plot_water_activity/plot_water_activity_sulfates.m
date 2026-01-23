close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T = 25; 
MWw = 18.015;

%% 1. Na2SO4
MW_Na2SO4 = 142.04;
% Fit valid: [0.8990, 0.9957]
RH_Na2SO4 = linspace(0.900, 0.995, 100); 

for i = 1:length(RH_Na2SO4)
    mf_salt_Na2SO4(i) = calculate_mf_Na2SO4(RH_Na2SO4(i));
    mf_water_Na2SO4(i) = 1 - mf_salt_Na2SO4(i);
    x_water_Na2SO4(i) = (mf_water_Na2SO4(i) / MWw) / ...
        ((mf_water_Na2SO4(i) / MWw) + (mf_salt_Na2SO4(i) / MW_Na2SO4));
end

%% 2. K2SO4
MW_K2SO4 = 174.26;
% Fit valid: [0.9720, 0.9958]
RH_K2SO4 = linspace(0.973, 0.995, 100);

for i = 1:length(RH_K2SO4)
    mf_salt_K2SO4(i) = calculate_mf_K2SO4(RH_K2SO4(i));
    mf_water_K2SO4(i) = 1 - mf_salt_K2SO4(i);
    x_water_K2SO4(i) = (mf_water_K2SO4(i) / MWw) / ...
        ((mf_water_K2SO4(i) / MWw) + (mf_salt_K2SO4(i) / MW_K2SO4));
end

%% 3. (NH4)2SO4
MW_NH42SO4 = 132.14;
% Fit valid: [0.8310, 0.9959]
RH_NH42SO4 = linspace(0.832, 0.995, 100);

for i = 1:length(RH_NH42SO4)
    mf_salt_NH42SO4(i) = calculate_mf_NH42SO4(RH_NH42SO4(i));
    mf_water_NH42SO4(i) = 1 - mf_salt_NH42SO4(i);
    x_water_NH42SO4(i) = (mf_water_NH42SO4(i) / MWw) / ...
        ((mf_water_NH42SO4(i) / MWw) + (mf_salt_NH42SO4(i) / MW_NH42SO4));
end

%% 4. MgSO4
MW_MgSO4 = 120.37;
% Fit valid: [0.9050, 0.9960]
RH_MgSO4 = linspace(0.906, 0.995, 100);

for i = 1:length(RH_MgSO4)
    mf_salt_MgSO4(i) = calculate_mf_MgSO4(RH_MgSO4(i));
    mf_water_MgSO4(i) = 1 - mf_salt_MgSO4(i);
    x_water_MgSO4(i) = (mf_water_MgSO4(i) / MWw) / ...
        ((mf_water_MgSO4(i) / MWw) + (mf_salt_MgSO4(i) / MW_MgSO4));
end

%% 5. MnSO4
MW_MnSO4 = 151.00;
% Fit valid: [0.8620 or 0.919, 0.9961]
RH_MnSO4 = linspace(0.92, 0.995, 100);

for i = 1:length(RH_MnSO4)
    mf_salt_MnSO4(i) = calculate_mf_MnSO4(RH_MnSO4(i));
    mf_water_MnSO4(i) = 1 - mf_salt_MnSO4(i);
    x_water_MnSO4(i) = (mf_water_MnSO4(i) / MWw) / ...
        ((mf_water_MnSO4(i) / MWw) + (mf_salt_MnSO4(i) / MW_MnSO4));
end

%% 6. Li2SO4
MW_Li2SO4 = 109.94;
% Fit valid: [0.8530, 0.9956]
RH_Li2SO4 = linspace(0.854, 0.995, 100);

for i = 1:length(RH_Li2SO4)
    mf_salt_Li2SO4(i) = calculate_mf_Li2SO4(RH_Li2SO4(i));
    mf_water_Li2SO4(i) = 1 - mf_salt_Li2SO4(i);
    x_water_Li2SO4(i) = (mf_water_Li2SO4(i) / MWw) / ...
        ((mf_water_Li2SO4(i) / MWw) + (mf_salt_Li2SO4(i) / MW_Li2SO4));
end

%% 7. NiSO4
MW_NiSO4 = 154.75;
% Fit valid: [0.9390, 0.9962]
RH_NiSO4 = linspace(0.940, 0.995, 100);

for i = 1:length(RH_NiSO4)
    mf_salt_NiSO4(i) = calculate_mf_NiSO4(RH_NiSO4(i));
    mf_water_NiSO4(i) = 1 - mf_salt_NiSO4(i);
    x_water_NiSO4(i) = (mf_water_NiSO4(i) / MWw) / ...
        ((mf_water_NiSO4(i) / MWw) + (mf_salt_NiSO4(i) / MW_NiSO4));
end

%% 8. CuSO4
MW_CuSO4 = 159.61;
% Fit valid: [0.9750, 0.9963]
RH_CuSO4 = linspace(0.976, 0.995, 100);

for i = 1:length(RH_CuSO4)
    mf_salt_CuSO4(i) = calculate_mf_CuSO4(RH_CuSO4(i));
    mf_water_CuSO4(i) = 1 - mf_salt_CuSO4(i);
    x_water_CuSO4(i) = (mf_water_CuSO4(i) / MWw) / ...
        ((mf_water_CuSO4(i) / MWw) + (mf_salt_CuSO4(i) / MW_CuSO4));
end

%% 9. ZnSO4
MW_ZnSO4 = 161.44;
% Fit valid: [0.9130, 0.9962]
RH_ZnSO4 = linspace(0.914, 0.995, 100);

for i = 1:length(RH_ZnSO4)
    mf_salt_ZnSO4(i) = calculate_mf_ZnSO4(RH_ZnSO4(i));
    mf_water_ZnSO4(i) = 1 - mf_salt_ZnSO4(i);
    x_water_ZnSO4(i) = (mf_water_ZnSO4(i) / MWw) / ...
        ((mf_water_ZnSO4(i) / MWw) + (mf_salt_ZnSO4(i) / MW_ZnSO4));
end


%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_Na2SO4  = RH_Na2SO4  ./ x_water_Na2SO4;
gamma_K2SO4   = RH_K2SO4   ./ x_water_K2SO4;
gamma_NH42SO4 = RH_NH42SO4 ./ x_water_NH42SO4;
gamma_MgSO4   = RH_MgSO4   ./ x_water_MgSO4;
gamma_MnSO4   = RH_MnSO4   ./ x_water_MnSO4;
gamma_Li2SO4  = RH_Li2SO4  ./ x_water_Li2SO4;
gamma_NiSO4   = RH_NiSO4   ./ x_water_NiSO4;
gamma_CuSO4   = RH_CuSO4   ./ x_water_CuSO4;
gamma_ZnSO4   = RH_ZnSO4   ./ x_water_ZnSO4;


%% Calculate Average Sulfate Fit
% Define the salts to include in the average
fit_salts = {'Na2SO4', 'K2SO4', 'NH42SO4', 'MgSO4', 'MnSO4', ...
             'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4'};

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

print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction_Sulfates'), '-dpng', '-r600')


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

print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH_Sulfates'), '-dpng', '-r600')

disp('Sulfate plots generated successfully!')