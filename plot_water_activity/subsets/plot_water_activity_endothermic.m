close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T = 25; 
MWw = 18.015;

%% NaCl
MW_NaCl = 58.443;
RH_NaCl = linspace(0.763, 0.9934, 100);

for i = 1:length(RH_NaCl)
    mf_salt_NaCl(i) = calculate_mf_NaCl(RH_NaCl(i));
    mf_water_NaCl(i) = 1 - mf_salt_NaCl(i);
    x_water_NaCl(i) = (mf_water_NaCl(i) / MWw) / ((mf_water_NaCl(i) / MWw) + (mf_salt_NaCl(i) / MW_NaCl));
end

%% KCl
MW_KCl = 74.551;
RH_KCl = linspace(0.853, 0.9935, 100);

for i = 1:length(RH_KCl)
    mf_salt_KCl(i) = calculate_mf_KCl(RH_KCl(i));
    mf_water_KCl(i) = 1 - mf_salt_KCl(i);
    x_water_KCl(i) = (mf_water_KCl(i) / MWw) / ((mf_water_KCl(i) / MWw) + (mf_salt_KCl(i) / MW_KCl));
end

%% NH4Cl
MW_NH4Cl = 53.491;
RH_NH4Cl = linspace(0.813, 0.993, 100);

for i = 1:length(RH_NH4Cl)
    mf_salt_NH4Cl(i) = calculate_mf_NH4Cl(RH_NH4Cl(i));
    mf_water_NH4Cl(i) = 1 - mf_salt_NH4Cl(i);
    x_water_NH4Cl(i) = (mf_water_NH4Cl(i) / MWw) / ((mf_water_NH4Cl(i) / MWw) + (mf_salt_NH4Cl(i) / MW_NH4Cl));
end

%% CsCl
MW_CsCl = 168.363;
RH_CsCl = linspace(0.817, 0.993, 100);

for i = 1:length(RH_CsCl)
    mf_salt_CsCl(i) = calculate_mf_CsCl(RH_CsCl(i));
    mf_water_CsCl(i) = 1 - mf_salt_CsCl(i);
    x_water_CsCl(i) = (mf_water_CsCl(i) / MWw) / ((mf_water_CsCl(i) / MWw) + (mf_salt_CsCl(i) / MW_CsCl));
end

%% NaNO3
MW_NaNO3 = 84.994;
RH_NaNO3 = linspace(0.9701, 0.9996, 100);

for i = 1:length(RH_NaNO3)
    mf_salt_NaNO3(i) = calculate_mf_NaNO3(RH_NaNO3(i));
    mf_water_NaNO3(i) = 1 - mf_salt_NaNO3(i);
    x_water_NaNO3(i) = (mf_water_NaNO3(i) / MWw) / ((mf_water_NaNO3(i) / MWw) + (mf_salt_NaNO3(i) / MW_NaNO3));
end

%% AgNO3
MW_AgNO3 = 169.87;
RH_AgNO3 = linspace(0.862, 0.986, 100);

for i = 1:length(RH_AgNO3)
    mf_salt_AgNO3(i) = calculate_mf_AgNO3(RH_AgNO3(i));
    mf_water_AgNO3(i) = 1 - mf_salt_AgNO3(i);
    x_water_AgNO3(i) = (mf_water_AgNO3(i) / MWw) / ((mf_water_AgNO3(i) / MWw) + (mf_salt_AgNO3(i) / MW_AgNO3));
end

%% KI
MW_KI = 165.998;
RH_KI = linspace(0.9671, 0.999, 100);

for i = 1:length(RH_KI)
    mf_salt_KI(i) = calculate_mf_KI(RH_KI(i));
    mf_water_KI(i) = 1 - mf_salt_KI(i);
    x_water_KI(i) = (mf_water_KI(i) / MWw) / ((mf_water_KI(i) / MWw) + (mf_salt_KI(i) / MW_KI));
end

%% LiNO3
MW_LiNO3 = 68.95;
RH_LiNO3 = linspace(0.7353, 0.9967, 100);

for i = 1:length(RH_LiNO3)
    mf_salt_LiNO3(i) = calculate_mf_LiNO3(RH_LiNO3(i));
    mf_water_LiNO3(i) = 1 - mf_salt_LiNO3(i);
    x_water_LiNO3(i) = (mf_water_LiNO3(i) / MWw) / ((mf_water_LiNO3(i) / MWw) + (mf_salt_LiNO3(i) / MW_LiNO3));
end

%% KNO3
MW_KNO3 = 101.10;
RH_KNO3 = linspace(0.9315, 0.9967, 100);

for i = 1:length(RH_KNO3)
    mf_salt_KNO3(i) = calculate_mf_KNO3(RH_KNO3(i));
    mf_water_KNO3(i) = 1 - mf_salt_KNO3(i);
    x_water_KNO3(i) = (mf_water_KNO3(i) / MWw) / ((mf_water_KNO3(i) / MWw) + (mf_salt_KNO3(i) / MW_KNO3));
end

%% NaClO4
MW_NaClO4 = 122.44;
RH_NaClO4 = linspace(0.7775, 0.9869, 100);

for i = 1:length(RH_NaClO4)
    mf_salt_NaClO4(i) = calculate_mf_NaClO4(RH_NaClO4(i));
    mf_water_NaClO4(i) = 1 - mf_salt_NaClO4(i);
    x_water_NaClO4(i) = (mf_water_NaClO4(i) / MWw) / ((mf_water_NaClO4(i) / MWw) + (mf_salt_NaClO4(i) / MW_NaClO4));
end

%% KClO3
MW_KClO3 = 122.55;
RH_KClO3 = linspace(0.9800, 0.9936, 100);

for i = 1:length(RH_KClO3)
    mf_salt_KClO3(i) = calculate_mf_KClO3(RH_KClO3(i));
    mf_water_KClO3(i) = 1 - mf_salt_KClO3(i);
    x_water_KClO3(i) = (mf_water_KClO3(i) / MWw) / ((mf_water_KClO3(i) / MWw) + (mf_salt_KClO3(i) / MW_KClO3));
end

%% NaBr
MW_NaBr = 102.89;
RH_NaBr = linspace(0.6133, 0.9290, 100);

for i = 1:length(RH_NaBr)
    mf_salt_NaBr(i) = calculate_mf_NaBr(RH_NaBr(i));
    mf_water_NaBr(i) = 1 - mf_salt_NaBr(i);
    x_water_NaBr(i) = (mf_water_NaBr(i) / MWw) / ((mf_water_NaBr(i) / MWw) + (mf_salt_NaBr(i) / MW_NaBr));
end

%% NaI
MW_NaI = 149.89;
RH_NaI = linspace(0.5801, 0.9669, 100);

for i = 1:length(RH_NaI)
    mf_salt_NaI(i) = calculate_mf_NaI(RH_NaI(i));
    mf_water_NaI(i) = 1 - mf_salt_NaI(i);
    x_water_NaI(i) = (mf_water_NaI(i) / MWw) / ((mf_water_NaI(i) / MWw) + (mf_salt_NaI(i) / MW_NaI));
end

%% KBr
MW_KBr = 119.00;
RH_KBr = linspace(0.8325, 0.9528, 100);

for i = 1:length(RH_KBr)
    mf_salt_KBr(i) = calculate_mf_KBr(RH_KBr(i));
    mf_water_KBr(i) = 1 - mf_salt_KBr(i);
    x_water_KBr(i) = (mf_water_KBr(i) / MWw) / ((mf_water_KBr(i) / MWw) + (mf_salt_KBr(i) / MW_KBr));
end

%% RbCl
MW_RbCl = 120.92;
RH_RbCl = linspace(0.7423, 0.9527, 100);

for i = 1:length(RH_RbCl)
    mf_salt_RbCl(i) = calculate_mf_RbCl(RH_RbCl(i));
    mf_water_RbCl(i) = 1 - mf_salt_RbCl(i);
    x_water_RbCl(i) = (mf_water_RbCl(i) / MWw) / ((mf_water_RbCl(i) / MWw) + (mf_salt_RbCl(i) / MW_RbCl));
end

%% CsBr
MW_CsBr = 212.81;
RH_CsBr = linspace(0.8475, 0.9482, 100);

for i = 1:length(RH_CsBr)
    mf_salt_CsBr(i) = calculate_mf_CsBr(RH_CsBr(i));
    mf_water_CsBr(i) = 1 - mf_salt_CsBr(i);
    x_water_CsBr(i) = (mf_water_CsBr(i) / MWw) / ((mf_water_CsBr(i) / MWw) + (mf_salt_CsBr(i) / MW_CsBr));
end

%% CsI
MW_CsI = 259.81;
RH_CsI = linspace(0.9124, 0.9624, 100);

for i = 1:length(RH_CsI)
    mf_salt_CsI(i) = calculate_mf_CsI(RH_CsI(i));
    mf_water_CsI(i) = 1 - mf_salt_CsI(i);
    x_water_CsI(i) = (mf_water_CsI(i) / MWw) / ((mf_water_CsI(i) / MWw) + (mf_salt_CsI(i) / MW_CsI));
end

%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_NaCl  = RH_NaCl  ./ x_water_NaCl;
gamma_KCl   = RH_KCl   ./ x_water_KCl;
gamma_NH4Cl = RH_NH4Cl ./ x_water_NH4Cl;
gamma_CsCl  = RH_CsCl  ./ x_water_CsCl;
gamma_NaNO3 = RH_NaNO3 ./ x_water_NaNO3;
gamma_AgNO3 = RH_AgNO3 ./ x_water_AgNO3;
gamma_KI    = RH_KI    ./ x_water_KI;
gamma_LiNO3 = RH_LiNO3 ./ x_water_LiNO3;
gamma_KNO3  = RH_KNO3  ./ x_water_KNO3;
gamma_NaClO4 = RH_NaClO4 ./ x_water_NaClO4;
gamma_KClO3 = RH_KClO3 ./ x_water_KClO3;
gamma_NaBr  = RH_NaBr  ./ x_water_NaBr;
gamma_NaI   = RH_NaI   ./ x_water_NaI;
gamma_KBr   = RH_KBr   ./ x_water_KBr;
gamma_RbCl  = RH_RbCl  ./ x_water_RbCl;
gamma_CsBr  = RH_CsBr  ./ x_water_CsBr;
gamma_CsI   = RH_CsI   ./ x_water_CsI;

%% Calculate Constrained Weighted Polynomial Fit
% Define the salts to include in the fit
fit_salts = {'NaCl', 'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', 'LiNO3', 'KNO3', 'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI'};

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