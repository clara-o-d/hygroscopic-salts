close all 
clear
clc 

% PlotDefaults_Slides()

T=25; 
MWw=18;

%% LiCl
MW_LiCl=42.4;
RH_LiCl=linspace(0.12,0.9);

for i=1:length(RH_LiCl)
   U_LiCl_gg(i)=1/calculate_mf_LiCl(RH_LiCl(i),T)-1;
end

U_LiCl_molmol=U_LiCl_gg*(2*MWw/(MW_LiCl))^-1;

%% LiOH
MW_LiOH=24;
RH_LiOH=linspace(0.85,0.9);

for i=1:length(RH_LiOH)
   U_LiOH_gg(i)=1/calculate_mf_LiOH(RH_LiOH(i))-1;
end

U_LiOH_molmol=U_LiOH_gg*(2*MWw/(MW_LiOH))^-1;

%% NaOH
MW_NaOH=40;
RH_NaOH=linspace(0.23,0.9);

for i=1:length(RH_NaOH)
   U_NaOH_gg(i)=1/calculate_mf_NaOH(RH_NaOH(i))-1;
end

U_NaOH_molmol=U_NaOH_gg*(2*MWw/(MW_NaOH))^-1;

%% HCl
MW_HCl=36.5;
RH_HCl=linspace(0.17,0.9);

for i=1:length(RH_HCl)
   U_HCl_gg(i)=1/calculate_mf_HCl(RH_HCl(i))-1;
end

U_HCl_molmol=U_HCl_gg*(2*MWw/(MW_HCl))^-1;


%% CaCl2
MW_CaCl2=111;
RH_CaCl2=linspace(0.31,0.9);

for i=1:length(RH_CaCl2)
   U_CaCl2_gg(i)=1/calculate_mf_CaCl(RH_CaCl2(i),T)-1;
end

U_CaCl2_molmol=U_CaCl2_gg*(3*MWw/(MW_CaCl2))^-1;

%% MgCl2
MW_MgCl2=95.2;
RH_MgCl2=linspace(0.33,0.9);

for i=1:length(RH_MgCl2)
   U_MgCl2_gg(i)=1/calculate_mf_MgCl(RH_MgCl2(i))-1;
end

U_MgCl2_molmol=U_MgCl2_gg*(3*MWw/(MW_MgCl2))^-1;

%% MgNO32
MW_MgNO32=148.3;
RH_MgNO32=linspace(0.55,0.9);

for i=1:length(RH_MgNO32)
   U_MgNO32_gg(i)=1/calculate_mf_MgNO3(RH_MgNO32(i))-1;
end

U_MgNO32_molmol=U_MgNO32_gg*(3*MWw/(MW_MgNO32))^-1;

%% LiBr
MW_LiBr=86.85;
RH_LiBr=linspace(0.07,0.9);

for i=1:length(RH_LiBr)
   U_LiBr_gg(i)=1/calculate_mf_LiBr(RH_LiBr(i))-1;
end

U_LiBr_molmol=U_LiBr_gg*(2*MWw/(MW_LiBr))^-1;

%% ZnCl

MW_ZnCl2=136.3;
RH_ZnCl2=linspace(0.07,0.8);

for i=1:length(RH_ZnCl2)
   U_ZnCl2_gg(i)=1/calculate_mf_ZnCl(RH_ZnCl2(i))-1;
end

U_ZnCl2_molmol=U_ZnCl2_gg*(3*MWw/(MW_ZnCl2))^-1;

%% ZnI
MW_ZnI2=319.18;
RH_ZnI2=linspace(0.25,0.9);

for i=1:length(RH_ZnI2)
   U_ZnI2_gg(i)=1/calculate_mf_ZnI(RH_ZnI2(i))-1;
end

U_ZnI2_molmol=U_ZnI2_gg*(3*MWw/(MW_ZnI2))^-1;

%% ZnBr2
MW_ZnBr2=225.2;
RH_ZnBr2=linspace(0.08,0.85);

for i=1:length(RH_ZnBr2)
   U_ZnBr2_gg(i)=1/calculate_mf_ZnBr(RH_ZnBr2(i))-1;
end

U_ZnBr2_molmol=U_ZnBr2_gg*(3*MWw/(MW_ZnBr2))^-1;

%% LiI
MW_LiI=133.85;
RH_LiI=linspace(0.18,0.9);

for i=1:length(RH_LiI)
   U_LiI_gg(i)=1/calculate_mf_LiI(RH_LiI(i))-1;
end

U_LiI_molmol=U_LiI_gg*(2*MWw/(MW_LiI))^-1;
%}

%% Ideal solution 
RHideal=linspace(0.01,0.9);
U_ideal_molmol=RHideal./(1-RHideal);

%% g/g Uptake
plot(RH_LiCl,U_LiCl_gg, 'color', [0, 0.5, 0])
hold on
plot(RH_CaCl2,U_CaCl2_gg, 'color',[0.9290 0.6940 0.1250])
plot(RH_MgCl2,U_MgCl2_gg,'color',[0.8500 0.3250 0.0980])
plot(RH_LiBr,U_LiBr_gg,'color',[0.3010, 0.7450, 0.9330])
%plot(RH_ZnCl2,U_ZnCl2_gg,'color',[0.6350 0.0780 0.1840])
%plot(RH_LiI,U_LiI_gg,'color',[0.4940 0.1840 0.5560])
%plot(RH_ZnBr2,U_ZnBr2_gg,'color',[0.4660 0.6740 0.1880])
%plot(RH_ZnI2,U_ZnI2_gg,'--','color',[0.9290 0.6940 0.1250])
plot(RH_HCl,U_HCl_gg)
%plot(RH_MgNO32,U_MgNO32_gg,'color',[0 1 1])
%plot(RH_LiOH,U_LiOH_gg,'color',[1 0 1])
%plot(RH_ZnI2,U_ZnI2_gg,'--','color',[0.9290 0.6940 0.1250])
plot(RH_NaOH,U_NaOH_gg,'color',[1 0 1])

xlabel('Relative Humidity (RH)')
ylabel('Uptake (g/g)')
title('Water Uptake vs Relative Humidity')
legend('LiCl', 'CaCl_2', 'MgCl_2', 'LiBr', 'HCl', 'Mg(NO_3)_2', 'NaOH', 'Location', 'best', 'Interpreter', 'tex')
grid on

set(gcf,'color','w');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4]; 
print('Uptake_screening_gg','-dpng','-r600')
close all 

%% Uptake mol/mol
plot(RH_LiCl*100,U_LiCl_molmol, 'color', [0, 0.5, 0])
hold on
%plot(RHideal*100,U_ideal_molmol)
plot(RH_CaCl2*100,U_CaCl2_molmol, 'color',[0.9290 0.6940 0.1250])
plot(RH_MgCl2*100,U_MgCl2_molmol,'color',[0.8500 0.3250 0.0980])
plot(RH_LiBr*100,U_LiBr_molmol,'color',[0.3010, 0.7450, 0.9330])
plot(RH_ZnCl2*100,U_ZnCl2_molmol,'color',[0.6350 0.0780 0.1840])
plot(RH_LiI*100,U_LiI_molmol,'color',[0.4940 0.1840 0.5560])
plot(RH_ZnBr2*100,U_ZnBr2_molmol,'color',[0.4660 0.6740 0.1880])
plot(RH_HCl*100,U_HCl_molmol)
%plot(RH_MgNO32*100,U_MgNO32_molmol,'color',[0 1 1])
plot(RH_LiOH*100,U_LiOH_molmol,'color',[1 0 1])
plot(RH_ZnI2*100,U_ZnI2_molmol,'--','color',[0.9290 0.6940 0.1250])
plot(RH_NaOH*100,U_NaOH_molmol,'color',[1 0 1])

xlabel('Relative Humidity (%)')
ylabel('Uptake (mol/mol)')
title('Water Uptake vs Relative Humidity')
legend('LiCl', 'CaCl_2', 'MgCl_2', 'LiBr', 'ZnCl_2', 'LiI', 'ZnBr_2', 'HCl', 'Mg(NO_3)_2', 'LiOH', 'ZnI_2', 'NaOH', 'Location', 'best', 'Interpreter', 'tex')
grid on

set(gcf,'color','w');

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5/1.6 4/1.6]; 
print('Uptake_screening_molmol','-dpng','-r600')
