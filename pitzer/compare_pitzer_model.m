close all 
clear
clc 

% Add calculate_mf folder to path
addpath('calculate_mf');

T = 25; % Temperature in Celsius
MWw = 18;

%% Calculate data for all salts

% LiCl
MW_LiCl = 42.4;
RH_LiCl = linspace(0.12, 0.9, 100);
for i = 1:length(RH_LiCl)
    mf_salt_LiCl(i) = calculate_mf_LiCl(RH_LiCl(i), T);
    mf_water_LiCl(i) = 1 - mf_salt_LiCl(i);
    x_water_LiCl(i) = (mf_water_LiCl(i) / MWw) / ((mf_water_LiCl(i) / MWw) + (mf_salt_LiCl(i) / MW_LiCl));
    % Calculate molality (mol salt / kg water)
    molality_LiCl(i) = (mf_salt_LiCl(i) / MW_LiCl) / (mf_water_LiCl(i) / 1000);
end

% CaCl2
MW_CaCl2 = 111;
RH_CaCl2 = linspace(0.31, 0.9, 100);
for i = 1:length(RH_CaCl2)
    mf_salt_CaCl2(i) = calculate_mf_CaCl(RH_CaCl2(i), T);
    mf_water_CaCl2(i) = 1 - mf_salt_CaCl2(i);
    x_water_CaCl2(i) = (mf_water_CaCl2(i) / MWw) / ((mf_water_CaCl2(i) / MWw) + (mf_salt_CaCl2(i) / MW_CaCl2));
    molality_CaCl2(i) = (mf_salt_CaCl2(i) / MW_CaCl2) / (mf_water_CaCl2(i) / 1000);
end

% MgCl2
MW_MgCl2 = 95.2;
RH_MgCl2 = linspace(0.33, 0.9, 100);
for i = 1:length(RH_MgCl2)
    mf_salt_MgCl2(i) = calculate_mf_MgCl(RH_MgCl2(i));
    mf_water_MgCl2(i) = 1 - mf_salt_MgCl2(i);
    x_water_MgCl2(i) = (mf_water_MgCl2(i) / MWw) / ((mf_water_MgCl2(i) / MWw) + (mf_salt_MgCl2(i) / MW_MgCl2));
    molality_MgCl2(i) = (mf_salt_MgCl2(i) / MW_MgCl2) / (mf_water_MgCl2(i) / 1000);
end

% LiBr
MW_LiBr = 86.85;
RH_LiBr = linspace(0.07, 0.9, 100);
for i = 1:length(RH_LiBr)
    mf_salt_LiBr(i) = calculate_mf_LiBr(RH_LiBr(i));
    mf_water_LiBr(i) = 1 - mf_salt_LiBr(i);
    x_water_LiBr(i) = (mf_water_LiBr(i) / MWw) / ((mf_water_LiBr(i) / MWw) + (mf_salt_LiBr(i) / MW_LiBr));
    molality_LiBr(i) = (mf_salt_LiBr(i) / MW_LiBr) / (mf_water_LiBr(i) / 1000);
end

% ZnCl2
MW_ZnCl2 = 136.3;
RH_ZnCl2 = linspace(0.07, 0.8, 100);
for i = 1:length(RH_ZnCl2)
    mf_salt_ZnCl2(i) = calculate_mf_ZnCl(RH_ZnCl2(i));
    mf_water_ZnCl2(i) = 1 - mf_salt_ZnCl2(i);
    x_water_ZnCl2(i) = (mf_water_ZnCl2(i) / MWw) / ((mf_water_ZnCl2(i) / MWw) + (mf_salt_ZnCl2(i) / MW_ZnCl2));
    molality_ZnCl2(i) = (mf_salt_ZnCl2(i) / MW_ZnCl2) / (mf_water_ZnCl2(i) / 1000);
end

%% Calculate Pitzer model predictions for each salt

% Pitzer parameters from literature (Pitzer & Mayorga, 1973; Kim & Frederick, 1988)
% For 1:1 electrolytes (LiCl, LiBr): beta0, beta1, Cphi
% For 2:1 electrolytes (CaCl2, MgCl2, ZnCl2): beta0, beta1, beta2, Cphi, alpha1, alpha2

% LiCl - 1:1 electrolyte
beta0_LiCl = 0.1494;
beta1_LiCl = 0.3074;
Cphi_LiCl = 0.0008;
for i = 1:length(molality_LiCl)
    aw_pitzer_LiCl(i) = pitzer_water_activity(molality_LiCl(i), 2, beta0_LiCl, beta1_LiCl, 0, Cphi_LiCl, 2.0, 0, T);
end

% CaCl2 - 2:1 electrolyte
beta0_CaCl2 = 0.3159;
beta1_CaCl2 = 1.614;
beta2_CaCl2 = -2.84;
Cphi_CaCl2 = -0.00034;
for i = 1:length(molality_CaCl2)
    aw_pitzer_CaCl2(i) = pitzer_water_activity(molality_CaCl2(i), 3, beta0_CaCl2, beta1_CaCl2, beta2_CaCl2, Cphi_CaCl2, 2.0, 12.0, T);
end

% MgCl2 - 2:1 electrolyte
beta0_MgCl2 = 0.3524;
beta1_MgCl2 = 1.6815;
beta2_MgCl2 = -3.10;
Cphi_MgCl2 = 0.0025;
for i = 1:length(molality_MgCl2)
    aw_pitzer_MgCl2(i) = pitzer_water_activity(molality_MgCl2(i), 3, beta0_MgCl2, beta1_MgCl2, beta2_MgCl2, Cphi_MgCl2, 2.0, 12.0, T);
end

% LiBr - 1:1 electrolyte
beta0_LiBr = 0.1748;
beta1_LiBr = 0.2547;
Cphi_LiBr = 0.0053;
for i = 1:length(molality_LiBr)
    aw_pitzer_LiBr(i) = pitzer_water_activity(molality_LiBr(i), 2, beta0_LiBr, beta1_LiBr, 0, Cphi_LiBr, 2.0, 0, T);
end

% ZnCl2 - 2:1 electrolyte
beta0_ZnCl2 = 0.3070;
beta1_ZnCl2 = 1.5790;
beta2_ZnCl2 = -28.5;
Cphi_ZnCl2 = -0.0048;
for i = 1:length(molality_ZnCl2)
    aw_pitzer_ZnCl2(i) = pitzer_water_activity(molality_ZnCl2(i), 3, beta0_ZnCl2, beta1_ZnCl2, beta2_ZnCl2, Cphi_ZnCl2, 2.0, 12.0, T);
end

%% Create comparison plots

% Figure 1: Water Activity Comparison (RH vs Molality)
figure('Position', [100, 100, 1200, 900]);

% LiCl
subplot(3, 2, 1)
hold on; grid on; box on;
plot(molality_LiCl, RH_LiCl, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0, 0.5, 0])
plot(molality_LiCl, aw_pitzer_LiCl, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Water Activity (a_w)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiCl', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% CaCl2
subplot(3, 2, 2)
hold on; grid on; box on;
plot(molality_CaCl2, RH_CaCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.9290 0.6940 0.1250])
plot(molality_CaCl2, aw_pitzer_CaCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Water Activity (a_w)', 'FontSize', 11, 'FontWeight', 'bold')
title('CaCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% MgCl2
subplot(3, 2, 3)
hold on; grid on; box on;
plot(molality_MgCl2, RH_MgCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.8500 0.3250 0.0980])
plot(molality_MgCl2, aw_pitzer_MgCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Water Activity (a_w)', 'FontSize', 11, 'FontWeight', 'bold')
title('MgCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% LiBr
subplot(3, 2, 4)
hold on; grid on; box on;
plot(molality_LiBr, RH_LiBr, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.3010, 0.7450, 0.9330])
plot(molality_LiBr, aw_pitzer_LiBr, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Water Activity (a_w)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiBr', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% ZnCl2
subplot(3, 2, 5)
hold on; grid on; box on;
plot(molality_ZnCl2, RH_ZnCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.6350 0.0780 0.1840])
plot(molality_ZnCl2, aw_pitzer_ZnCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Water Activity (a_w)', 'FontSize', 11, 'FontWeight', 'bold')
title('ZnCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

sgtitle('Water Activity: Data vs Pitzer Model', 'FontSize', 16, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 9]; 
print('figures/Pitzer_Comparison_Water_Activity', '-dpng', '-r600')
savefig('figures/Pitzer_Comparison_Water_Activity.fig')

%% Figure 2: Osmotic Coefficient Comparison
figure('Position', [100, 100, 1200, 900]);

% Calculate osmotic coefficients from your data
phi_data_LiCl = -log(RH_LiCl) ./ (0.018015 * 2 * molality_LiCl);
phi_data_CaCl2 = -log(RH_CaCl2) ./ (0.018015 * 3 * molality_CaCl2);
phi_data_MgCl2 = -log(RH_MgCl2) ./ (0.018015 * 3 * molality_MgCl2);
phi_data_LiBr = -log(RH_LiBr) ./ (0.018015 * 2 * molality_LiBr);
phi_data_ZnCl2 = -log(RH_ZnCl2) ./ (0.018015 * 3 * molality_ZnCl2);

% Calculate osmotic coefficients from Pitzer model
phi_pitzer_LiCl = -log(aw_pitzer_LiCl) ./ (0.018015 * 2 * molality_LiCl);
phi_pitzer_CaCl2 = -log(aw_pitzer_CaCl2) ./ (0.018015 * 3 * molality_CaCl2);
phi_pitzer_MgCl2 = -log(aw_pitzer_MgCl2) ./ (0.018015 * 3 * molality_MgCl2);
phi_pitzer_LiBr = -log(aw_pitzer_LiBr) ./ (0.018015 * 2 * molality_LiBr);
phi_pitzer_ZnCl2 = -log(aw_pitzer_ZnCl2) ./ (0.018015 * 3 * molality_ZnCl2);

% LiCl
subplot(3, 2, 1)
hold on; grid on; box on;
plot(molality_LiCl, phi_data_LiCl, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0, 0.5, 0])
plot(molality_LiCl, phi_pitzer_LiCl, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Osmotic Coefficient (\phi)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiCl', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% CaCl2
subplot(3, 2, 2)
hold on; grid on; box on;
plot(molality_CaCl2, phi_data_CaCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.9290 0.6940 0.1250])
plot(molality_CaCl2, phi_pitzer_CaCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Osmotic Coefficient (\phi)', 'FontSize', 11, 'FontWeight', 'bold')
title('CaCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% MgCl2
subplot(3, 2, 3)
hold on; grid on; box on;
plot(molality_MgCl2, phi_data_MgCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.8500 0.3250 0.0980])
plot(molality_MgCl2, phi_pitzer_MgCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Osmotic Coefficient (\phi)', 'FontSize', 11, 'FontWeight', 'bold')
title('MgCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% LiBr
subplot(3, 2, 4)
hold on; grid on; box on;
plot(molality_LiBr, phi_data_LiBr, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.3010, 0.7450, 0.9330])
plot(molality_LiBr, phi_pitzer_LiBr, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Osmotic Coefficient (\phi)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiBr', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

% ZnCl2
subplot(3, 2, 5)
hold on; grid on; box on;
plot(molality_ZnCl2, phi_data_ZnCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Your Data', 'color', [0.6350 0.0780 0.1840])
plot(molality_ZnCl2, phi_pitzer_ZnCl2, 's--', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Pitzer Model', 'color', [0.8, 0, 0])
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Osmotic Coefficient (\phi)', 'FontSize', 11, 'FontWeight', 'bold')
title('ZnCl_2', 'FontSize', 13, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 9)
set(gca, 'FontSize', 10)

sgtitle('Osmotic Coefficient: Data vs Pitzer Model', 'FontSize', 16, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 9]; 
print('figures/Pitzer_Comparison_Osmotic_Coefficient', '-dpng', '-r600')
savefig('figures/Pitzer_Comparison_Osmotic_Coefficient.fig')

%% Figure 3: Residual Analysis (Relative Error)
figure('Position', [100, 100, 1200, 900]);

% Calculate relative errors
rel_error_LiCl = 100 * (RH_LiCl - aw_pitzer_LiCl) ./ RH_LiCl;
rel_error_CaCl2 = 100 * (RH_CaCl2 - aw_pitzer_CaCl2) ./ RH_CaCl2;
rel_error_MgCl2 = 100 * (RH_MgCl2 - aw_pitzer_MgCl2) ./ RH_MgCl2;
rel_error_LiBr = 100 * (RH_LiBr - aw_pitzer_LiBr) ./ RH_LiBr;
rel_error_ZnCl2 = 100 * (RH_ZnCl2 - aw_pitzer_ZnCl2) ./ RH_ZnCl2;

% LiCl
subplot(3, 2, 1)
hold on; grid on; box on;
plot(molality_LiCl, rel_error_LiCl, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'color', [0, 0.5, 0])
plot([min(molality_LiCl) max(molality_LiCl)], [0 0], 'k--', 'LineWidth', 1.5)
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Relative Error (%)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiCl', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% CaCl2
subplot(3, 2, 2)
hold on; grid on; box on;
plot(molality_CaCl2, rel_error_CaCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'color', [0.9290 0.6940 0.1250])
plot([min(molality_CaCl2) max(molality_CaCl2)], [0 0], 'k--', 'LineWidth', 1.5)
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Relative Error (%)', 'FontSize', 11, 'FontWeight', 'bold')
title('CaCl_2', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% MgCl2
subplot(3, 2, 3)
hold on; grid on; box on;
plot(molality_MgCl2, rel_error_MgCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'color', [0.8500 0.3250 0.0980])
plot([min(molality_MgCl2) max(molality_MgCl2)], [0 0], 'k--', 'LineWidth', 1.5)
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Relative Error (%)', 'FontSize', 11, 'FontWeight', 'bold')
title('MgCl_2', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% LiBr
subplot(3, 2, 4)
hold on; grid on; box on;
plot(molality_LiBr, rel_error_LiBr, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'color', [0.3010, 0.7450, 0.9330])
plot([min(molality_LiBr) max(molality_LiBr)], [0 0], 'k--', 'LineWidth', 1.5)
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Relative Error (%)', 'FontSize', 11, 'FontWeight', 'bold')
title('LiBr', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% ZnCl2
subplot(3, 2, 5)
hold on; grid on; box on;
plot(molality_ZnCl2, rel_error_ZnCl2, 'o-', 'LineWidth', 2, 'MarkerSize', 4, 'color', [0.6350 0.0780 0.1840])
plot([min(molality_ZnCl2) max(molality_ZnCl2)], [0 0], 'k--', 'LineWidth', 1.5)
xlabel('Molality (mol/kg)', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Relative Error (%)', 'FontSize', 11, 'FontWeight', 'bold')
title('ZnCl_2', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

sgtitle('Relative Error: (Data - Pitzer) / Data × 100%', 'FontSize', 16, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 9]; 
print('figures/Pitzer_Comparison_Residuals', '-dpng', '-r600')
savefig('figures/Pitzer_Comparison_Residuals.fig')

disp('Pitzer model comparison completed successfully!')
disp('Generated files:')
disp('  - figures/Pitzer_Comparison_Water_Activity.png')
disp('  - figures/Pitzer_Comparison_Osmotic_Coefficient.png')
disp('  - figures/Pitzer_Comparison_Residuals.png')

%% Print statistical summary
fprintf('\n=== STATISTICAL SUMMARY ===\n\n')
fprintf('Salt       | Mean Error (%%) | RMSE (%%) | Max Error (%%)\n')
fprintf('-----------|----------------|----------|-------------\n')
fprintf('LiCl       | %8.3f       | %7.3f  | %7.3f\n', mean(abs(rel_error_LiCl)), rms(rel_error_LiCl), max(abs(rel_error_LiCl)))
fprintf('CaCl2      | %8.3f       | %7.3f  | %7.3f\n', mean(abs(rel_error_CaCl2)), rms(rel_error_CaCl2), max(abs(rel_error_CaCl2)))
fprintf('MgCl2      | %8.3f       | %7.3f  | %7.3f\n', mean(abs(rel_error_MgCl2)), rms(rel_error_MgCl2), max(abs(rel_error_MgCl2)))
fprintf('LiBr       | %8.3f       | %7.3f  | %7.3f\n', mean(abs(rel_error_LiBr)), rms(rel_error_LiBr), max(abs(rel_error_LiBr)))
fprintf('ZnCl2      | %8.3f       | %7.3f  | %7.3f\n', mean(abs(rel_error_ZnCl2)), rms(rel_error_ZnCl2), max(abs(rel_error_ZnCl2)))
fprintf('\n')
