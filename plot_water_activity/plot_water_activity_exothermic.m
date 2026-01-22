close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T = 25; 
MWw = 18;

%% LiCl
MW_LiCl = 42.4;
RH_LiCl = linspace(0.12, 0.9, 100);

for i = 1:length(RH_LiCl)
    mf_salt_LiCl(i) = calculate_mf_LiCl(RH_LiCl(i), T);
    mf_water_LiCl(i) = 1 - mf_salt_LiCl(i);
    x_water_LiCl(i) = (mf_water_LiCl(i) / MWw) / ((mf_water_LiCl(i) / MWw) + (mf_salt_LiCl(i) / MW_LiCl));
end

%% LiOH
MW_LiOH = 24;
RH_LiOH = linspace(0.85, 0.9, 100);

for i = 1:length(RH_LiOH)
    mf_salt_LiOH(i) = calculate_mf_LiOH(RH_LiOH(i));
    mf_water_LiOH(i) = 1 - mf_salt_LiOH(i);
    x_water_LiOH(i) = (mf_water_LiOH(i) / MWw) / ((mf_water_LiOH(i) / MWw) + (mf_salt_LiOH(i) / MW_LiOH));
end

%% NaOH
MW_NaOH = 40;
RH_NaOH = linspace(0.23, 0.9, 100);

for i = 1:length(RH_NaOH)
    mf_salt_NaOH(i) = calculate_mf_NaOH(RH_NaOH(i));
    mf_water_NaOH(i) = 1 - mf_salt_NaOH(i);
    x_water_NaOH(i) = (mf_water_NaOH(i) / MWw) / ((mf_water_NaOH(i) / MWw) + (mf_salt_NaOH(i) / MW_NaOH));
end

%% HCl
MW_HCl = 36.5;
RH_HCl = linspace(0.17, 0.9, 100);

for i = 1:length(RH_HCl)
    mf_salt_HCl(i) = calculate_mf_HCl(RH_HCl(i));
    mf_water_HCl(i) = 1 - mf_salt_HCl(i);
    x_water_HCl(i) = (mf_water_HCl(i) / MWw) / ((mf_water_HCl(i) / MWw) + (mf_salt_HCl(i) / MW_HCl));
end

%% CaCl2
MW_CaCl2 = 111;
RH_CaCl2 = linspace(0.31, 0.9, 100);

for i = 1:length(RH_CaCl2)
    mf_salt_CaCl2(i) = calculate_mf_CaCl(RH_CaCl2(i), T);
    mf_water_CaCl2(i) = 1 - mf_salt_CaCl2(i);
    x_water_CaCl2(i) = (mf_water_CaCl2(i) / MWw) / ((mf_water_CaCl2(i) / MWw) + (mf_salt_CaCl2(i) / MW_CaCl2));
end

%% MgCl2
MW_MgCl2 = 95.2;
RH_MgCl2 = linspace(0.33, 0.9, 100);

for i = 1:length(RH_MgCl2)
    mf_salt_MgCl2(i) = calculate_mf_MgCl(RH_MgCl2(i));
    mf_water_MgCl2(i) = 1 - mf_salt_MgCl2(i);
    x_water_MgCl2(i) = (mf_water_MgCl2(i) / MWw) / ((mf_water_MgCl2(i) / MWw) + (mf_salt_MgCl2(i) / MW_MgCl2));
end

%% MgNO32
MW_MgNO32 = 148.3;
RH_MgNO32 = linspace(0.55, 0.9, 100);

for i = 1:length(RH_MgNO32)
    mf_salt_MgNO32(i) = calculate_mf_MgNO3(RH_MgNO32(i));
    mf_water_MgNO32(i) = 1 - mf_salt_MgNO32(i);
    x_water_MgNO32(i) = (mf_water_MgNO32(i) / MWw) / ((mf_water_MgNO32(i) / MWw) + (mf_salt_MgNO32(i) / MW_MgNO32));
end

%% LiBr
MW_LiBr = 86.85;
RH_LiBr = linspace(0.07, 0.9, 100);

for i = 1:length(RH_LiBr)
    mf_salt_LiBr(i) = calculate_mf_LiBr(RH_LiBr(i));
    mf_water_LiBr(i) = 1 - mf_salt_LiBr(i);
    x_water_LiBr(i) = (mf_water_LiBr(i) / MWw) / ((mf_water_LiBr(i) / MWw) + (mf_salt_LiBr(i) / MW_LiBr));
end

%% ZnCl2
MW_ZnCl2 = 136.3;
RH_ZnCl2 = linspace(0.07, 0.8, 100);

for i = 1:length(RH_ZnCl2)
    mf_salt_ZnCl2(i) = calculate_mf_ZnCl(RH_ZnCl2(i));
    mf_water_ZnCl2(i) = 1 - mf_salt_ZnCl2(i);
    x_water_ZnCl2(i) = (mf_water_ZnCl2(i) / MWw) / ((mf_water_ZnCl2(i) / MWw) + (mf_salt_ZnCl2(i) / MW_ZnCl2));
end

%% ZnI2
MW_ZnI2 = 319.18;
RH_ZnI2 = linspace(0.25, 0.9, 100);

for i = 1:length(RH_ZnI2)
    mf_salt_ZnI2(i) = calculate_mf_ZnI(RH_ZnI2(i));
    mf_water_ZnI2(i) = 1 - mf_salt_ZnI2(i);
    x_water_ZnI2(i) = (mf_water_ZnI2(i) / MWw) / ((mf_water_ZnI2(i) / MWw) + (mf_salt_ZnI2(i) / MW_ZnI2));
end

%% ZnBr2
MW_ZnBr2 = 225.2;
RH_ZnBr2 = linspace(0.08, 0.85, 100);

for i = 1:length(RH_ZnBr2)
    mf_salt_ZnBr2(i) = calculate_mf_ZnBr(RH_ZnBr2(i));
    mf_water_ZnBr2(i) = 1 - mf_salt_ZnBr2(i);
    x_water_ZnBr2(i) = (mf_water_ZnBr2(i) / MWw) / ((mf_water_ZnBr2(i) / MWw) + (mf_salt_ZnBr2(i) / MW_ZnBr2));
end

%% LiI
MW_LiI = 133.85;
RH_LiI = linspace(0.18, 0.9, 100);

for i = 1:length(RH_LiI)
    mf_salt_LiI(i) = calculate_mf_LiI(RH_LiI(i));
    mf_water_LiI(i) = 1 - mf_salt_LiI(i);
    x_water_LiI(i) = (mf_water_LiI(i) / MWw) / ((mf_water_LiI(i) / MWw) + (mf_salt_LiI(i) / MW_LiI));
end

%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_LiCl = RH_LiCl ./ x_water_LiCl;
gamma_CaCl2 = RH_CaCl2 ./ x_water_CaCl2;
gamma_MgCl2 = RH_MgCl2 ./ x_water_MgCl2;
gamma_LiBr = RH_LiBr ./ x_water_LiBr;
gamma_ZnCl2 = RH_ZnCl2 ./ x_water_ZnCl2;
gamma_LiI = RH_LiI ./ x_water_LiI;
gamma_ZnBr2 = RH_ZnBr2 ./ x_water_ZnBr2;
gamma_ZnI2 = RH_ZnI2 ./ x_water_ZnI2;
gamma_HCl = RH_HCl ./ x_water_HCl;
gamma_MgNO32 = RH_MgNO32 ./ x_water_MgNO32;
gamma_LiOH = RH_LiOH ./ x_water_LiOH;
gamma_NaOH = RH_NaOH ./ x_water_NaOH;

x_water_LiCl
plot(x_water_LiCl, gamma_LiCl, 'LineWidth', 2.5, 'DisplayName', 'LiCl', 'color', [0, 0.5, 0])

%% Calculate Constrained Weighted Polynomial Fit
% Define the salts to include in the fit (Exclude MgNO32)
fit_salts = {'LiCl', 'CaCl2', 'MgCl2', 'LiBr', 'ZnCl2', ...
             'LiI', 'ZnBr2', 'ZnI2', 'HCl', 'LiOH', 'NaOH'};

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
print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction'), '-dpng', '-r600')
savefig(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_Mole_Fraction.fig'))

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
print(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH'), '-dpng', '-r600')
savefig(fullfile(filepath, '..', 'figures', 'activity_coefficient_subsets', 'Activity_Coefficient_vs_RH.fig'))

disp('Plots generated successfully!')
disp('  - figures/activity_coefficient_subsets/Activity_Coefficient_vs_Mole_Fraction.png')
disp('  - figures/activity_coefficient_subsets/Activity_Coefficient_vs_RH.png')