close all 
clear
clc 

% Add calculate_mf folder to path if needed
addpath('calculate_mf');

T = 25; 
MWw = 18.015;

%% NaCl
MW_NaCl = 58.443;
% Fit valid > 0.762 (User req: 73)
RH_NaCl = linspace(0.765, 0.99, 100);

for i = 1:length(RH_NaCl)
   U_NaCl_gg(i) = 1/calculate_mf_NaCl_(RH_NaCl(i)) - 1;
end
% 1:1 salt -> 2 ions
U_NaCl_molmol = U_NaCl_gg * (2 * MWw / MW_NaCl)^-1;

%% KCl
MW_KCl = 74.551;
% Fit valid > 0.852 (User req: 83)
RH_KCl = linspace(0.855, 0.99, 100);

for i = 1:length(RH_KCl)
   U_KCl_gg(i) = 1/calculate_mf_KCl_(RH_KCl(i)) - 1;
end
U_KCl_molmol = U_KCl_gg * (2 * MWw / MW_KCl)^-1;

%% NH4Cl
MW_NH4Cl = 53.491;
% Fit valid > 0.812 (User req: 76)
RH_NH4Cl = linspace(0.815, 0.99, 100);

for i = 1:length(RH_NH4Cl)
   U_NH4Cl_gg(i) = 1/calculate_mf_NH4Cl_(RH_NH4Cl(i)) - 1;
end
U_NH4Cl_molmol = U_NH4Cl_gg * (2 * MWw / MW_NH4Cl)^-1;

%% CsCl
MW_CsCl = 168.363;
% Fit valid > 0.817 (User req: 60)
RH_CsCl = linspace(0.82, 0.99, 100);

for i = 1:length(RH_CsCl)
   U_CsCl_gg(i) = 1/calculate_mf_CsCl_(RH_CsCl(i)) - 1;
end
U_CsCl_molmol = U_CsCl_gg * (2 * MWw / MW_CsCl)^-1;

%% NaNO3
MW_NaNO3 = 84.994;
% Fit valid > 0.970 (User req: 71) - RESTRICTED TO FIT RANGE
RH_NaNO3 = linspace(0.971, 0.995, 100);

for i = 1:length(RH_NaNO3)
   U_NaNO3_gg(i) = 1/calculate_mf_NaNO3_(RH_NaNO3(i)) - 1;
end
U_NaNO3_molmol = U_NaNO3_gg * (2 * MWw / MW_NaNO3)^-1;

%% AgNO3
MW_AgNO3 = 169.87;
% Fit valid > 0.862
RH_AgNO3 = linspace(0.865, 0.99, 100);

for i = 1:length(RH_AgNO3)
   U_AgNO3_gg(i) = 1/calculate_mf_AgNO3_(RH_AgNO3(i)) - 1;
end
U_AgNO3_molmol = U_AgNO3_gg * (2 * MWw / MW_AgNO3)^-1;

%% KI
MW_KI = 165.998;
% Fit valid > 0.973 (User req: 65) - RESTRICTED TO FIT RANGE
RH_KI = linspace(0.975, 0.995, 100);

for i = 1:length(RH_KI)
   U_KI_gg(i) = 1/calculate_mf_KI_(RH_KI(i)) - 1;
end
U_KI_molmol = U_KI_gg * (2 * MWw / MW_KI)^-1;


%% PLOTTING

% --- g/g Uptake ---
figure('Position', [100, 100, 800, 600]);
hold on
plot(RH_NaCl,  U_NaCl_gg,  'LineWidth', 2, 'DisplayName', 'NaCl')
plot(RH_KCl,   U_KCl_gg,   'LineWidth', 2, 'DisplayName', 'KCl')
plot(RH_NH4Cl, U_NH4Cl_gg, 'LineWidth', 2, 'DisplayName', 'NH_4Cl')
plot(RH_CsCl,  U_CsCl_gg,  'LineWidth', 2, 'DisplayName', 'CsCl')
plot(RH_NaNO3, U_NaNO3_gg, 'LineWidth', 2, 'DisplayName', 'NaNO_3')
plot(RH_AgNO3, U_AgNO3_gg, 'LineWidth', 2, 'DisplayName', 'AgNO_3')
plot(RH_KI,    U_KI_gg,    'LineWidth', 2, 'DisplayName', 'KI')

xlabel('Relative Humidity (RH)')
ylabel('Uptake (g/g)')
title('Endothermic Salts: Water Uptake (Mass Basis)')
legend('Location', 'best', 'Interpreter', 'tex')
grid on
set(gcf,'color','w');

% Save figure
print('Uptake_Endothermic_gg','-dpng','-r600')


% --- mol/mol Uptake ---
figure('Position', [150, 150, 800, 600]);
hold on
plot(RH_NaCl*100,  U_NaCl_molmol,  'LineWidth', 2, 'DisplayName', 'NaCl')
plot(RH_KCl*100,   U_KCl_molmol,   'LineWidth', 2, 'DisplayName', 'KCl')
plot(RH_NH4Cl*100, U_NH4Cl_molmol, 'LineWidth', 2, 'DisplayName', 'NH_4Cl')
plot(RH_CsCl*100,  U_CsCl_molmol,  'LineWidth', 2, 'DisplayName', 'CsCl')
plot(RH_NaNO3*100, U_NaNO3_molmol, 'LineWidth', 2, 'DisplayName', 'NaNO_3')
plot(RH_AgNO3*100, U_AgNO3_molmol, 'LineWidth', 2, 'DisplayName', 'AgNO_3')
plot(RH_KI*100,    U_KI_molmol,    'LineWidth', 2, 'DisplayName', 'KI')

xlabel('Relative Humidity (%)')
ylabel('Uptake (mol water / mol dissociation particle)')
title('Endothermic Salts: Water Uptake (Molar Basis)')
legend('Location', 'best', 'Interpreter', 'tex')
grid on
set(gcf,'color','w');

% Save figure
print('Uptake_Endothermic_molmol','-dpng','-r600')

disp('Endothermic salt plots generated.')