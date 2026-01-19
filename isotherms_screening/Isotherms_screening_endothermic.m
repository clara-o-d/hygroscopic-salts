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

%% LiNO3
MW_LiNO3 = 68.95;
% Fit valid > 0.7353
RH_LiNO3 = linspace(0.736, 0.99, 100);

for i = 1:length(RH_LiNO3)
   U_LiNO3_gg(i) = 1/calculate_mf_LiNO3(RH_LiNO3(i)) - 1;
end
U_LiNO3_molmol = U_LiNO3_gg * (2 * MWw / MW_LiNO3)^-1;

%% KNO3
MW_KNO3 = 101.10;
% Fit valid > 0.9315 - RESTRICTED TO FIT RANGE
RH_KNO3 = linspace(0.932, 0.995, 100);

for i = 1:length(RH_KNO3)
   U_KNO3_gg(i) = 1/calculate_mf_KNO3(RH_KNO3(i)) - 1;
end
U_KNO3_molmol = U_KNO3_gg * (2 * MWw / MW_KNO3)^-1;

%% NaClO4
MW_NaClO4 = 122.44;
% Fit valid > 0.7775
RH_NaClO4 = linspace(0.778, 0.99, 100);

for i = 1:length(RH_NaClO4)
   U_NaClO4_gg(i) = 1/calculate_mf_NaClO4(RH_NaClO4(i)) - 1;
end
U_NaClO4_molmol = U_NaClO4_gg * (2 * MWw / MW_NaClO4)^-1;

%% KClO3
MW_KClO3 = 122.55;
% Fit valid > 0.9800 - RESTRICTED TO FIT RANGE
RH_KClO3 = linspace(0.981, 0.995, 100);

for i = 1:length(RH_KClO3)
   U_KClO3_gg(i) = 1/calculate_mf_KClO3(RH_KClO3(i)) - 1;
end
U_KClO3_molmol = U_KClO3_gg * (2 * MWw / MW_KClO3)^-1;

%% NaBr
MW_NaBr = 102.89;
% Fit valid > 0.6133
RH_NaBr = linspace(0.614, 0.99, 100);

for i = 1:length(RH_NaBr)
   U_NaBr_gg(i) = 1/calculate_mf_NaBr(RH_NaBr(i)) - 1;
end
U_NaBr_molmol = U_NaBr_gg * (2 * MWw / MW_NaBr)^-1;

%% NaI
MW_NaI = 149.89;
% Fit valid > 0.5801
RH_NaI = linspace(0.581, 0.99, 100);

for i = 1:length(RH_NaI)
   U_NaI_gg(i) = 1/calculate_mf_NaI(RH_NaI(i)) - 1;
end
U_NaI_molmol = U_NaI_gg * (2 * MWw / MW_NaI)^-1;

%% KBr
MW_KBr = 119.00;
% Fit valid > 0.8325
RH_KBr = linspace(0.833, 0.99, 100);

for i = 1:length(RH_KBr)
   U_KBr_gg(i) = 1/calculate_mf_KBr(RH_KBr(i)) - 1;
end
U_KBr_molmol = U_KBr_gg * (2 * MWw / MW_KBr)^-1;

%% RbCl
MW_RbCl = 120.92;
% Fit valid > 0.7423
RH_RbCl = linspace(0.743, 0.99, 100);

for i = 1:length(RH_RbCl)
   U_RbCl_gg(i) = 1/calculate_mf_RbCl(RH_RbCl(i)) - 1;
end
U_RbCl_molmol = U_RbCl_gg * (2 * MWw / MW_RbCl)^-1;

%% CsBr
MW_CsBr = 212.81;
% Fit valid > 0.8475
RH_CsBr = linspace(0.848, 0.99, 100);

for i = 1:length(RH_CsBr)
   U_CsBr_gg(i) = 1/calculate_mf_CsBr(RH_CsBr(i)) - 1;
end
U_CsBr_molmol = U_CsBr_gg * (2 * MWw / MW_CsBr)^-1;

%% CsI
MW_CsI = 259.81;
% Fit valid > 0.9124 - RESTRICTED TO FIT RANGE
RH_CsI = linspace(0.913, 0.995, 100);

for i = 1:length(RH_CsI)
   U_CsI_gg(i) = 1/calculate_mf_CsI(RH_CsI(i)) - 1;
end
U_CsI_molmol = U_CsI_gg * (2 * MWw / MW_CsI)^-1;


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
plot(RH_LiNO3, U_LiNO3_gg, 'LineWidth', 2, 'DisplayName', 'LiNO_3')
plot(RH_KNO3,  U_KNO3_gg,  'LineWidth', 2, 'DisplayName', 'KNO_3')
plot(RH_NaClO4, U_NaClO4_gg, 'LineWidth', 2, 'DisplayName', 'NaClO_4')
plot(RH_KClO3, U_KClO3_gg, 'LineWidth', 2, 'DisplayName', 'KClO_3')
plot(RH_NaBr,  U_NaBr_gg,  'LineWidth', 2, 'DisplayName', 'NaBr')
plot(RH_NaI,   U_NaI_gg,   'LineWidth', 2, 'DisplayName', 'NaI')
plot(RH_KBr,   U_KBr_gg,   'LineWidth', 2, 'DisplayName', 'KBr')
plot(RH_RbCl,  U_RbCl_gg,  'LineWidth', 2, 'DisplayName', 'RbCl')
plot(RH_CsBr,  U_CsBr_gg,  'LineWidth', 2, 'DisplayName', 'CsBr')
plot(RH_CsI,   U_CsI_gg,   'LineWidth', 2, 'DisplayName', 'CsI')

xlabel('Relative Humidity (RH)')
ylabel('Uptake (g/g)')
title('Endothermic Salts: Water Uptake (Mass Basis)')
legend('Location', 'best', 'Interpreter', 'tex')
grid on
set(gcf,'color','w');

% Save figure
print(fullfile(filepath, '..', 'figures', 'Uptake_Endothermic_gg'),'-dpng','-r600')


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
plot(RH_LiNO3*100, U_LiNO3_molmol, 'LineWidth', 2, 'DisplayName', 'LiNO_3')
plot(RH_KNO3*100,  U_KNO3_molmol,  'LineWidth', 2, 'DisplayName', 'KNO_3')
plot(RH_NaClO4*100, U_NaClO4_molmol, 'LineWidth', 2, 'DisplayName', 'NaClO_4')
plot(RH_KClO3*100, U_KClO3_molmol, 'LineWidth', 2, 'DisplayName', 'KClO_3')
plot(RH_NaBr*100,  U_NaBr_molmol,  'LineWidth', 2, 'DisplayName', 'NaBr')
plot(RH_NaI*100,   U_NaI_molmol,   'LineWidth', 2, 'DisplayName', 'NaI')
plot(RH_KBr*100,   U_KBr_molmol,   'LineWidth', 2, 'DisplayName', 'KBr')
plot(RH_RbCl*100,  U_RbCl_molmol,  'LineWidth', 2, 'DisplayName', 'RbCl')
plot(RH_CsBr*100,  U_CsBr_molmol,  'LineWidth', 2, 'DisplayName', 'CsBr')
plot(RH_CsI*100,   U_CsI_molmol,   'LineWidth', 2, 'DisplayName', 'CsI')

xlabel('Relative Humidity (%)')
ylabel('Uptake (mol water / mol dissociation particle)')
title('Endothermic Salts: Water Uptake (Molar Basis)')
legend('Location', 'best', 'Interpreter', 'tex')
grid on
set(gcf,'color','w');

% Save figure
print(fullfile(filepath, '..', 'figures', 'Uptake_Endothermic_molmol'),'-dpng','-r600')

disp('Endothermic salt plots generated.')