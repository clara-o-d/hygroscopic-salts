close all 
clear
clc 

% Add calculate_mf folder to path
addpath('calculate_mf');

T = 25; 
MWw = 18.015;

%% NaCl
MW_NaCl = 58.443;
% Output Range from fit: [0.7620, 0.9934]
RH_NaCl = linspace(0.765, 0.9934, 100);

for i = 1:length(RH_NaCl)
    % Note: T is implicit in the specific function for this salt (25C)
    mf_salt_NaCl(i) = calculate_mf_NaCl_(RH_NaCl(i));
    mf_water_NaCl(i) = 1 - mf_salt_NaCl(i);
    x_water_NaCl(i) = (mf_water_NaCl(i) / MWw) / ((mf_water_NaCl(i) / MWw) + (mf_salt_NaCl(i) / MW_NaCl));
end

%% KCl
MW_KCl = 74.551;
% Output Range from fit: [0.8520, 0.9935]
RH_KCl = linspace(0.852, 0.9935, 100);

for i = 1:length(RH_KCl)
    mf_salt_KCl(i) = calculate_mf_KCl_(RH_KCl(i));
    mf_water_KCl(i) = 1 - mf_salt_KCl(i);
    x_water_KCl(i) = (mf_water_KCl(i) / MWw) / ((mf_water_KCl(i) / MWw) + (mf_salt_KCl(i) / MW_KCl));
end

%% NH4Cl
MW_NH4Cl = 53.491;
% Output Range from fit: [0.8120, 0.9930]
RH_NH4Cl = linspace(0.812, 0.993, 100);

for i = 1:length(RH_NH4Cl)
    mf_salt_NH4Cl(i) = calculate_mf_NH4Cl_(RH_NH4Cl(i));
    mf_water_NH4Cl(i) = 1 - mf_salt_NH4Cl(i);
    x_water_NH4Cl(i) = (mf_water_NH4Cl(i) / MWw) / ((mf_water_NH4Cl(i) / MWw) + (mf_salt_NH4Cl(i) / MW_NH4Cl));
end

%% CsCl
MW_CsCl = 168.363;
% Output Range from fit: [0.8170, 0.9930]
RH_CsCl = linspace(0.817, 0.993, 100);

for i = 1:length(RH_CsCl)
    mf_salt_CsCl(i) = calculate_mf_CsCl_(RH_CsCl(i));
    mf_water_CsCl(i) = 1 - mf_salt_CsCl(i);
    x_water_CsCl(i) = (mf_water_CsCl(i) / MWw) / ((mf_water_CsCl(i) / MWw) + (mf_salt_CsCl(i) / MW_CsCl));
end

%% NaNO3
MW_NaNO3 = 84.994;
% Output Range from fit: [0.9701, 0.9996]
RH_NaNO3 = linspace(0.9701, 0.9996, 100);

for i = 1:length(RH_NaNO3)
    mf_salt_NaNO3(i) = calculate_mf_NaNO3_(RH_NaNO3(i));
    mf_water_NaNO3(i) = 1 - mf_salt_NaNO3(i);
    x_water_NaNO3(i) = (mf_water_NaNO3(i) / MWw) / ((mf_water_NaNO3(i) / MWw) + (mf_salt_NaNO3(i) / MW_NaNO3));
end

%% AgNO3
MW_AgNO3 = 169.87;
% Output Range from fit: [0.8620, 0.9860]
RH_AgNO3 = linspace(0.862, 0.986, 100);

for i = 1:length(RH_AgNO3)
    mf_salt_AgNO3(i) = calculate_mf_AgNO3_(RH_AgNO3(i));
    mf_water_AgNO3(i) = 1 - mf_salt_AgNO3(i);
    x_water_AgNO3(i) = (mf_water_AgNO3(i) / MWw) / ((mf_water_AgNO3(i) / MWw) + (mf_salt_AgNO3(i) / MW_AgNO3));
end

%% KI
MW_KI = 165.998;
% Output Range from fit: [0.9737, 1.0000] ?
RH_KI = linspace(0.97, 0.999, 100);

for i = 1:length(RH_KI)
    mf_salt_KI(i) = calculate_mf_KI_(RH_KI(i));
    mf_water_KI(i) = 1 - mf_salt_KI(i);
    x_water_KI(i) = (mf_water_KI(i) / MWw) / ((mf_water_KI(i) / MWw) + (mf_salt_KI(i) / MW_KI));
end

%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_NaCl  = RH_NaCl  ./ x_water_NaCl;
gamma_KCl   = RH_KCl   ./ x_water_KCl;
gamma_NH4Cl = RH_NH4Cl ./ x_water_NH4Cl;
% gamma_CsCl  = RH_CsCl  ./ x_water_CsCl;
gamma_NaNO3 = RH_NaNO3 ./ x_water_NaNO3;
gamma_AgNO3 = RH_AgNO3 ./ x_water_AgNO3;
gamma_KI    = RH_KI    ./ x_water_KI;

%% Calculate Constrained Weighted Polynomial Fit
% Define the salts to include in the fit
fit_salts = {'NaCl', 'KCl', 'NH4Cl', 'NaNO3', 'AgNO3', 'KI'};

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

%% FIGURE 1: Activity Coefficient vs Mole Fraction
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;
plot(x_water_NaCl,  gamma_NaCl,  'LineWidth', 2.5, 'DisplayName', 'NaCl',  'color', [0, 0.4470, 0.7410])
plot(x_water_KCl,   gamma_KCl,   'LineWidth', 2.5, 'DisplayName', 'KCl',   'color', [0.8500 0.3250 0.0980])
plot(x_water_NH4Cl, gamma_NH4Cl, 'LineWidth', 2.5, 'DisplayName', 'NH_4Cl','color', [0.9290 0.6940 0.1250])
% plot(x_water_CsCl,  gamma_CsCl,  'LineWidth', 2.5, 'DisplayName', 'CsCl',  'color', [0.4940 0.1840 0.5560])
plot(x_water_NaNO3, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3','color', [0.4660 0.6740 0.1880])
plot(x_water_AgNO3, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3','color', [0.3010 0.7450 0.9330])
plot(x_water_KI,    gamma_KI,    'LineWidth', 2.5, 'DisplayName', 'KI',    'color', [0.6350 0.0780 0.1840])

plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Mole Fraction', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 11)
xlim([0.8 1.0]) % Adjusted xlim since these salts are mostly soluble/high aw
ylim([0.9 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print('figures/Activity_Coefficient_vs_Mole_Fraction_New', '-dpng', '-r600')

%% FIGURE 2: Activity Coefficient vs Relative Humidity
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;
plot(RH_NaCl*100,  gamma_NaCl,  'LineWidth', 2.5, 'DisplayName', 'NaCl',  'color', [0, 0.4470, 0.7410])
plot(RH_KCl*100,   gamma_KCl,   'LineWidth', 2.5, 'DisplayName', 'KCl',   'color', [0.8500 0.3250 0.0980])
plot(RH_NH4Cl*100, gamma_NH4Cl, 'LineWidth', 2.5, 'DisplayName', 'NH_4Cl','color', [0.9290 0.6940 0.1250])
% plot(RH_CsCl*100,  gamma_CsCl,  'LineWidth', 2.5, 'DisplayName', 'CsCl',  'color', [0.4940 0.1840 0.5560])
plot(RH_NaNO3*100, gamma_NaNO3, 'LineWidth', 2.5, 'DisplayName', 'NaNO_3','color', [0.4660 0.6740 0.1880])
plot(RH_AgNO3*100, gamma_AgNO3, 'LineWidth', 2.5, 'DisplayName', 'AgNO_3','color', [0.3010 0.7450 0.9330])
plot(RH_KI*100,    gamma_KI,    'LineWidth', 2.5, 'DisplayName', 'KI',    'color', [0.6350 0.0780 0.1840])

% Plot the Average Fit Line
% plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', 'Average')

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Relative Humidity', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northwest', 'FontSize', 11)
xlim([70 100]) % Adjusted xlim to focus on the data range
ylim([0.9 1.1])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print('figures/Activity_Coefficient_vs_RH_New', '-dpng', '-r600')

disp('Plots generated successfully!')
disp('  - figures/Activity_Coefficient_vs_Mole_Fraction_New.png')
disp('  - figures/Activity_Coefficient_vs_RH_New.png')