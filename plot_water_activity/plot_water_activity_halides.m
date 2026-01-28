close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T = 25; 
MWw = 18.015;

%% 1. NaBr
MW_NaBr = 102.89;
% Fit valid: [0.6133, 0.9290]
RH_NaBr = linspace(0.614, 0.928, 100); 

for i = 1:length(RH_NaBr)
    mf_salt_NaBr(i) = calculate_mf_NaBr(RH_NaBr(i));
    mf_water_NaBr(i) = 1 - mf_salt_NaBr(i);
    x_water_NaBr(i) = (mf_water_NaBr(i) / MWw) / ...
        ((mf_water_NaBr(i) / MWw) + (mf_salt_NaBr(i) / MW_NaBr));
end

%% 2. NaI
MW_NaI = 149.89;
% Fit valid: [0.5801, 0.9669]
RH_NaI = linspace(0.581, 0.966, 100);

for i = 1:length(RH_NaI)
    mf_salt_NaI(i) = calculate_mf_NaI(RH_NaI(i));
    mf_water_NaI(i) = 1 - mf_salt_NaI(i);
    x_water_NaI(i) = (mf_water_NaI(i) / MWw) / ...
        ((mf_water_NaI(i) / MWw) + (mf_salt_NaI(i) / MW_NaI));
end

%% 3. KBr
MW_KBr = 119.00;
% Fit valid: [0.8325, 0.9528]
RH_KBr = linspace(0.833, 0.952, 100);

for i = 1:length(RH_KBr)
    mf_salt_KBr(i) = calculate_mf_KBr(RH_KBr(i));
    mf_water_KBr(i) = 1 - mf_salt_KBr(i);
    x_water_KBr(i) = (mf_water_KBr(i) / MWw) / ...
        ((mf_water_KBr(i) / MWw) + (mf_salt_KBr(i) / MW_KBr));
end

%% 4. RbCl
MW_RbCl = 120.92;
% Fit valid: [0.7423, 0.9527]
RH_RbCl = linspace(0.743, 0.952, 100);

for i = 1:length(RH_RbCl)
    mf_salt_RbCl(i) = calculate_mf_RbCl(RH_RbCl(i));
    mf_water_RbCl(i) = 1 - mf_salt_RbCl(i);
    x_water_RbCl(i) = (mf_water_RbCl(i) / MWw) / ...
        ((mf_water_RbCl(i) / MWw) + (mf_salt_RbCl(i) / MW_RbCl));
end

%% 5. CsBr
MW_CsBr = 212.81;
% Fit valid: [0.8475, 0.9482]
RH_CsBr = linspace(0.848, 0.948, 100);

for i = 1:length(RH_CsBr)
    mf_salt_CsBr(i) = calculate_mf_CsBr(RH_CsBr(i));
    mf_water_CsBr(i) = 1 - mf_salt_CsBr(i);
    x_water_CsBr(i) = (mf_water_CsBr(i) / MWw) / ...
        ((mf_water_CsBr(i) / MWw) + (mf_salt_CsBr(i) / MW_CsBr));
end

%% 6. CsI
MW_CsI = 259.81;
% Fit valid: [0.9124, 0.9624]
RH_CsI = linspace(0.913, 0.962, 100);

for i = 1:length(RH_CsI)
    mf_salt_CsI(i) = calculate_mf_CsI(RH_CsI(i));
    mf_water_CsI(i) = 1 - mf_salt_CsI(i);
    x_water_CsI(i) = (mf_water_CsI(i) / MWw) / ...
        ((mf_water_CsI(i) / MWw) + (mf_salt_CsI(i) / MW_CsI));
end

%% 7. CaBr2
MW_CaBr2 = 199.89;
% Fit valid: [0.6395, 0.9540]
RH_CaBr2 = linspace(0.640, 0.953, 100);

for i = 1:length(RH_CaBr2)
    mf_salt_CaBr2(i) = calculate_mf_CaBr2(RH_CaBr2(i));
    mf_water_CaBr2(i) = 1 - mf_salt_CaBr2(i);
    x_water_CaBr2(i) = (mf_water_CaBr2(i) / MWw) / ...
        ((mf_water_CaBr2(i) / MWw) + (mf_salt_CaBr2(i) / MW_CaBr2));
end

%% 8. CaI2
MW_CaI2 = 293.89;
% Fit valid: [0.8321, 0.9524]
RH_CaI2 = linspace(0.833, 0.952, 100);

for i = 1:length(RH_CaI2)
    mf_salt_CaI2(i) = calculate_mf_CaI2(RH_CaI2(i));
    mf_water_CaI2(i) = 1 - mf_salt_CaI2(i);
    x_water_CaI2(i) = (mf_water_CaI2(i) / MWw) / ...
        ((mf_water_CaI2(i) / MWw) + (mf_salt_CaI2(i) / MW_CaI2));
end

%% 9. SrCl2
MW_SrCl2 = 158.53;
% Fit valid: [0.8059, 0.9778]
RH_SrCl2 = linspace(0.806, 0.977, 100);

for i = 1:length(RH_SrCl2)
    mf_salt_SrCl2(i) = calculate_mf_SrCl2(RH_SrCl2(i));
    mf_water_SrCl2(i) = 1 - mf_salt_SrCl2(i);
    x_water_SrCl2(i) = (mf_water_SrCl2(i) / MWw) / ...
        ((mf_water_SrCl2(i) / MWw) + (mf_salt_SrCl2(i) / MW_SrCl2));
end

%% 10. SrBr2
MW_SrBr2 = 247.43;
% Fit valid: [0.7776, 0.9571]
RH_SrBr2 = linspace(0.778, 0.957, 100);

for i = 1:length(RH_SrBr2)
    mf_salt_SrBr2(i) = calculate_mf_SrBr2(RH_SrBr2(i));
    mf_water_SrBr2(i) = 1 - mf_salt_SrBr2(i);
    x_water_SrBr2(i) = (mf_water_SrBr2(i) / MWw) / ...
        ((mf_water_SrBr2(i) / MWw) + (mf_salt_SrBr2(i) / MW_SrBr2));
end

%% 11. SrI2
MW_SrI2 = 341.43;
% Fit valid: [0.6785, 0.9569]
RH_SrI2 = linspace(0.679, 0.956, 100);

for i = 1:length(RH_SrI2)
    mf_salt_SrI2(i) = calculate_mf_SrI2(RH_SrI2(i));
    mf_water_SrI2(i) = 1 - mf_salt_SrI2(i);
    x_water_SrI2(i) = (mf_water_SrI2(i) / MWw) / ...
        ((mf_water_SrI2(i) / MWw) + (mf_salt_SrI2(i) / MW_SrI2));
end

%% 12. BaCl2
MW_BaCl2 = 208.23;
% Fit valid: [0.9375, 0.9731]
RH_BaCl2 = linspace(0.938, 0.973, 100);

for i = 1:length(RH_BaCl2)
    mf_salt_BaCl2(i) = calculate_mf_BaCl2(RH_BaCl2(i));
    mf_water_BaCl2(i) = 1 - mf_salt_BaCl2(i);
    x_water_BaCl2(i) = (mf_water_BaCl2(i) / MWw) / ...
        ((mf_water_BaCl2(i) / MWw) + (mf_salt_BaCl2(i) / MW_BaCl2));
end

%% 13. BaBr2
MW_BaBr2 = 297.14;
% Fit valid: [0.8221, 0.9587]
RH_BaBr2 = linspace(0.823, 0.958, 100);

for i = 1:length(RH_BaBr2)
    mf_salt_BaBr2(i) = calculate_mf_BaBr2(RH_BaBr2(i));
    mf_water_BaBr2(i) = 1 - mf_salt_BaBr2(i);
    x_water_BaBr2(i) = (mf_water_BaBr2(i) / MWw) / ...
        ((mf_water_BaBr2(i) / MWw) + (mf_salt_BaBr2(i) / MW_BaBr2));
end


%% Calculate Activity Coefficients (gamma_w = a_w / x_w)
gamma_NaBr = RH_NaBr ./ x_water_NaBr;
gamma_NaI = RH_NaI ./ x_water_NaI;
gamma_KBr = RH_KBr ./ x_water_KBr;
gamma_RbCl = RH_RbCl ./ x_water_RbCl;
gamma_CsBr = RH_CsBr ./ x_water_CsBr;
gamma_CsI = RH_CsI ./ x_water_CsI;
gamma_CaBr2 = RH_CaBr2 ./ x_water_CaBr2;
gamma_CaI2 = RH_CaI2 ./ x_water_CaI2;
gamma_SrCl2 = RH_SrCl2 ./ x_water_SrCl2;
gamma_SrBr2 = RH_SrBr2 ./ x_water_SrBr2;
gamma_SrI2 = RH_SrI2 ./ x_water_SrI2;
gamma_BaCl2 = RH_BaCl2 ./ x_water_BaCl2;
gamma_BaBr2 = RH_BaBr2 ./ x_water_BaBr2;


%% Calculate Average Halide Fit
% Define the salts to include in the average
fit_salts = {'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI', 'CaBr2', 'CaI2', 'SrCl2', 'SrBr2', 'SrI2', 'BaCl2', 'BaBr2'};

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
