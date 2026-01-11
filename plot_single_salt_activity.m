function plot_single_salt_activity(salt_name)

T = 25;
MWw = 18;

switch lower(salt_name)
    case 'licl'
        MW_salt = 42.4;
        RH_range = linspace(0.12, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_LiCl(RH_range(i), T);
        end
        color = [0, 0.5, 0];
        display_name = 'LiCl';
        
    case 'cacl2'
        MW_salt = 111;
        RH_range = linspace(0.31, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_CaCl(RH_range(i), T);
        end
        color = [0.9290 0.6940 0.1250];
        display_name = 'CaCl_2';
        
    case 'mgcl2'
        MW_salt = 95.2;
        RH_range = linspace(0.33, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_MgCl(RH_range(i));
        end
        color = [0.8500 0.3250 0.0980];
        display_name = 'MgCl_2';
        
    case 'libr'
        MW_salt = 86.85;
        RH_range = linspace(0.07, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_LiBr(RH_range(i));
        end
        color = [0.3010, 0.7450, 0.9330];
        display_name = 'LiBr';
        
    case 'zncl2'
        MW_salt = 136.3;
        RH_range = linspace(0.07, 0.8, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_ZnCl(RH_range(i));
        end
        color = [0.6350 0.0780 0.1840];
        display_name = 'ZnCl_2';
        
    case 'lii'
        MW_salt = 133.85;
        RH_range = linspace(0.18, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_LiI(RH_range(i));
        end
        color = [0.4940 0.1840 0.5560];
        display_name = 'LiI';
        
    case 'znbr2'
        MW_salt = 225.2;
        RH_range = linspace(0.08, 0.85, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_ZnBr(RH_range(i));
        end
        color = [0.4660 0.6740 0.1880];
        display_name = 'ZnBr_2';
        
    case 'zni2'
        MW_salt = 319.18;
        RH_range = linspace(0.25, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_ZnI(RH_range(i));
        end
        color = [0.9290 0.6940 0.1250];
        display_name = 'ZnI_2';
        
    case 'hcl'
        MW_salt = 36.5;
        RH_range = linspace(0.17, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_HCl(RH_range(i));
        end
        color = [0, 0.4470, 0.7410];
        display_name = 'HCl';
        
    case 'mgno32'
        MW_salt = 148.3;
        RH_range = linspace(0.55, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_MgNO3(RH_range(i));
        end
        color = [0 1 1];
        display_name = 'Mg(NO_3)_2';
        
    case 'lioh'
        MW_salt = 24;
        RH_range = linspace(0.85, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_LiOH(RH_range(i));
        end
        color = [1 0 1];
        display_name = 'LiOH';
        
    case 'naoh'
        MW_salt = 40;
        RH_range = linspace(0.23, 0.9, 100);
        for i = 1:length(RH_range)
            mf_salt(i) = calculate_mf_NaOH(RH_range(i));
        end
        color = [0.75 0 0.75];
        display_name = 'NaOH';
        
    otherwise
        error('Unknown salt. Options: LiCl, CaCl2, MgCl2, LiBr, ZnCl2, LiI, ZnBr2, ZnI2, HCl, MgNO32, LiOH, NaOH');
end

mf_water = 1 - mf_salt;
x_water = (mf_water / MWw) ./ ((mf_water / MWw) + (mf_salt / MW_salt));
a_water = RH_range;
gamma_water = a_water ./ x_water;

x_ideal = linspace(min(x_water), 1.0, 100);
a_ideal = x_ideal;

figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1)
hold on; grid on; box on;
plot(x_water, a_water, 'LineWidth', 3, 'Color', color, 'DisplayName', display_name)
plot(x_ideal, a_ideal, 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (Raoult''s Law)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity (a_w = RH)', 'FontSize', 14, 'FontWeight', 'bold')
title(['Water Activity for ' display_name], 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'southeast', 'FontSize', 12)
xlim([min(x_water)-0.05 1.0])
ylim([0 1.0])
set(gca, 'FontSize', 12)

subplot(1, 2, 2)
hold on; grid on; box on;
plot(x_water, gamma_water, 'LineWidth', 3, 'Color', color, 'DisplayName', display_name)
plot([min(x_water)-0.05 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(['Activity Coefficient for ' display_name], 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'northeast', 'FontSize', 12)
xlim([min(x_water)-0.05 1.0])
set(gca, 'FontSize', 12)

set(gcf, 'color', 'w');

filename = ['Water_Activity_' salt_name];
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 5];
print(filename, '-dpng', '-r600')
savefig([filename '.fig'])

disp(['Plot saved as: ' filename '.png'])

fprintf('\nThermodynamic Data for %s:\n', display_name)
fprintf('----------------------------------------\n')
fprintf('  RH Range: %.2f - %.2f\n', min(RH_range), max(RH_range))
fprintf('  x_w Range: %.3f - %.3f\n', min(x_water), max(x_water))
fprintf('  gamma_w Range: %.3f - %.3f\n', min(gamma_water), max(gamma_water))
fprintf('  Avg deviation from ideality: %.2f%%\n', (mean(gamma_water) - 1) * 100)
fprintf('----------------------------------------\n')

end

