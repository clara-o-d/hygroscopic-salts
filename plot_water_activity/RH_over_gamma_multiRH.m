wclose all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'mole_fraction_water');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Define all salts
salt_data = {
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1};
    {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1};
    {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1};
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2};
    {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
    {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1};
    {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2};
    {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
    {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1};
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1};
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1};
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2};
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2};
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2};
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2};
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2};
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2};
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2};
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
};

% Process each salt
num_points = 100;
all_salts_data = struct();

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    n_cat = salt_data{s}{7};
    n_an  = salt_data{s}{8};
    nu = n_cat + n_an;
    
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water_mol = zeros(size(RH_vec));
    x_water_ion = zeros(size(RH_vec));
    
    success = true;
    for i = 1:length(RH_vec)
        try
            if func_args == 1
                mf_salt(i) = feval(func_name, RH_vec(i), T);
            else
                mf_salt(i) = feval(func_name, RH_vec(i));
            end
            mf_water(i) = 1 - mf_salt(i);
            n_w = mf_water(i) / MWw;
            n_s = mf_salt(i) / MW;
            x_water_mol(i) = n_w / (n_w + n_s);
            x_water_ion(i) = n_w / (n_w + (nu * n_s));
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        all_salts_data.(salt_name).RH = RH_vec;
        all_salts_data.(salt_name).x_water_mol = x_water_mol;
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).MW = MW;
        all_salts_data.(salt_name).nu = nu;
        mf_water_from_x_mol = (x_water_mol * MWw) ./ (x_water_mol * MWw + (1 - x_water_mol) * MW);
        r_ion = (nu * x_water_ion) ./ (1 - x_water_ion);
        mf_water_from_x_ion = (r_ion * MWw) ./ (r_ion * MWw + MW);
        all_salts_data.(salt_name).mf_water_from_x_mol = mf_water_from_x_mol;
        all_salts_data.(salt_name).mf_water_from_x_ion = mf_water_from_x_ion;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

%% Generate colors for salts
colors = zeros(num_salts, 3);
golden_angle = 0.38196601125; 

for i = 1:num_salts
    hue = mod((i-1) * golden_angle, 1.0);
    sat_cycle = mod(i, 3);
    if sat_cycle == 1, saturation = 0.85 + 0.15 * sin(i * 0.5); 
    elseif sat_cycle == 2, saturation = 0.70 + 0.15 * cos(i * 0.3);
    else, saturation = 0.55 + 0.20 * sin(i * 0.7);
    end
    saturation = max(0.5, min(1.0, saturation));
    
    val_cycle = mod(i, 4);
    if val_cycle == 1, value = 0.90 + 0.10 * cos(i * 0.4);
    elseif val_cycle == 2, value = 0.80 + 0.15 * sin(i * 0.6);
    elseif val_cycle == 3, value = 0.70 + 0.20 * cos(i * 0.5);
    else, value = 0.65 + 0.25 * sin(i * 0.8);
    end
    value = max(0.65, min(1.0, value));
    colors(i,:) = hsv2rgb([hue, saturation, value]);
end

%% UPTAKE VS RH: MOLAR AND MASS BASIS (All Salts)

figure('Position', [50, 50, 2400, 900]);

% Subplot 1: Molar Basis (Molecular)
subplot(1, 2, 1);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.x_water_mol, 'LineWidth', 2.5, ...
         'DisplayName', salt_name, 'color', colors(s,:));
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar Basis (Molecular)', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 10, 'NumColumns', 3)
xlim([0 100])
ylim([0 1])
set(gca, 'FontSize', 18)

% Subplot 2: Mass Basis (from Molecular)
subplot(1, 2, 2);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.mf_water_from_x_mol, 'LineWidth', 2.5, ...
         'DisplayName', salt_name, 'color', colors(s,:));
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Mass-Based Water Uptake (Mass Fraction)', 'FontSize', 20, 'FontWeight', 'bold')
title('Mass Basis (from Molecular)', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 10, 'NumColumns', 3)
xlim([0 100])
ylim([0 1])
set(gca, 'FontSize', 18)

sgtitle('Water Uptake vs RH: Molar and Mass Basis for All Salts', 'FontSize', 26, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

filename = 'Water_Uptake_vs_RH_Molar_and_Mass_All_Salts';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

%% FACETED PLOTS AT MULTIPLE RH VALUES

target_RH_values = [0.50, 0.60, 0.70, 0.80, 0.90];
num_RH_panels = length(target_RH_values);
data_at_RH = struct();

for rh_idx = 1:num_RH_panels
    target_RH = target_RH_values(rh_idx);
    x_water_mol_at_RH = [];
    x_water_ion_at_RH = [];
    mf_water_mol_at_RH = [];
    mf_water_ion_at_RH = [];
    salt_names_at_RH = {};
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        if target_RH >= min(data.RH) && target_RH <= max(data.RH)
            [~, idx] = min(abs(data.RH - target_RH));
            x_water_mol_at_RH = [x_water_mol_at_RH; data.x_water_mol(idx)];
            x_water_ion_at_RH = [x_water_ion_at_RH; data.x_water_ion(idx)];
            mf_water_mol_at_RH = [mf_water_mol_at_RH; data.mf_water_from_x_mol(idx)];
            mf_water_ion_at_RH = [mf_water_ion_at_RH; data.mf_water_from_x_ion(idx)];
            salt_names_at_RH{end+1} = salt_name;
        end
    end
    
    data_at_RH(rh_idx).target_RH = target_RH;
    data_at_RH(rh_idx).x_water_mol = x_water_mol_at_RH;
    data_at_RH(rh_idx).x_water_ion = x_water_ion_at_RH;
    data_at_RH(rh_idx).mf_water_mol = mf_water_mol_at_RH;
    data_at_RH(rh_idx).mf_water_ion = mf_water_ion_at_RH;
    data_at_RH(rh_idx).salt_names = salt_names_at_RH;
    fprintf('Found %d salts at %.0f%% RH\n', length(salt_names_at_RH), target_RH*100);
end

%% PLOT 1: Molecular Basis - Molar and Mass

figure('Position', [50, 50, 2000, 1200]);

for rh_idx = 1:num_RH_panels
    subplot(2, 3, rh_idx);
    hold on; grid on; box on;
    
    target_RH = data_at_RH(rh_idx).target_RH;
    x_water_mol = data_at_RH(rh_idx).x_water_mol;
    mf_water_mol = data_at_RH(rh_idx).mf_water_mol;
    salt_names_at_RH = data_at_RH(rh_idx).salt_names;
    
    if isempty(salt_names_at_RH)
        continue;
    end
    
    x_positions = 1:length(salt_names_at_RH);
    
    yyaxis left
    scatter(x_positions, x_water_mol, 80, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    ylabel('Molar Uptake', 'FontSize', 14, 'FontWeight', 'bold')
    ylim([0 1.05])
    
    yyaxis right
    scatter(x_positions, mf_water_mol, 80, 'filled', 'MarkerFaceColor', [0.8 0.4 0.2], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'Marker', 's');
    ylabel('Mass Fraction', 'FontSize', 14, 'FontWeight', 'bold')
    ylim([0 1.05])
    
    title(sprintf('RH = %.0f%%', target_RH*100), 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Salt', 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_RH, 'XTickLabelRotation', 45)
    xlim([0 length(salt_names_at_RH)+1])
    set(gca, 'FontSize', 11)
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
end

sgtitle('Molar and Mass-Based Water Uptake (Molecular Basis) - Faceted by RH', 'FontSize', 20, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

subplot(2, 3, 6);
axis off; hold on;
h1 = scatter(1, 1, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
h2 = scatter(1, 0.5, 100, 'filled', 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's');
legend([h1, h2], {'Molar (circles)', 'Mass (squares)'}, 'Location', 'best', 'FontSize', 14)
xlim([0 2]); ylim([0 2]);

filename = 'Molar_and_Mass_Water_Uptake_Molecular_Faceted_by_RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

%% PLOT 2: Ionic Basis - Molar and Mass

figure('Position', [50, 50, 2000, 1200]);

for rh_idx = 1:num_RH_panels
    subplot(2, 3, rh_idx);
    hold on; grid on; box on;
    
    target_RH = data_at_RH(rh_idx).target_RH;
    x_water_ion = data_at_RH(rh_idx).x_water_ion;
    mf_water_ion = data_at_RH(rh_idx).mf_water_ion;
    salt_names_at_RH = data_at_RH(rh_idx).salt_names;
    
    if isempty(salt_names_at_RH)
        continue;
    end
    
    x_positions = 1:length(salt_names_at_RH);
    
    yyaxis left
    scatter(x_positions, x_water_ion, 80, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    ylabel('Molar Uptake', 'FontSize', 14, 'FontWeight', 'bold')
    ylim([0 1.05])
    
    yyaxis right
    scatter(x_positions, mf_water_ion, 80, 'filled', 'MarkerFaceColor', [0.4 0.8 0.2], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'Marker', 's');
    ylabel('Mass Fraction', 'FontSize', 14, 'FontWeight', 'bold')
    ylim([0 1.05])
    
    title(sprintf('RH = %.0f%%', target_RH*100), 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Salt', 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_RH, 'XTickLabelRotation', 45)
    xlim([0 length(salt_names_at_RH)+1])
    set(gca, 'FontSize', 11)
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
end

sgtitle('Molar and Mass-Based Water Uptake (Ionic Basis) - Faceted by RH', 'FontSize', 20, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

subplot(2, 3, 6);
axis off; hold on;
h1 = scatter(1, 1, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
h2 = scatter(1, 0.5, 100, 'filled', 'MarkerFaceColor', [0.4 0.8 0.2], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's');
legend([h1, h2], {'Molar (circles)', 'Mass (squares)'}, 'Location', 'best', 'FontSize', 14)
xlim([0 2]); ylim([0 2]);

filename = 'Molar_and_Mass_Water_Uptake_Ionic_Faceted_by_RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

%% PLOT 3: Molar Uptake - Molecular vs Ionic

figure('Position', [50, 50, 2000, 1200]);

for rh_idx = 1:num_RH_panels
    subplot(2, 3, rh_idx);
    hold on; grid on; box on;
    
    target_RH = data_at_RH(rh_idx).target_RH;
    x_water_mol = data_at_RH(rh_idx).x_water_mol;
    x_water_ion = data_at_RH(rh_idx).x_water_ion;
    salt_names_at_RH = data_at_RH(rh_idx).salt_names;
    
    if isempty(salt_names_at_RH)
        continue;
    end
    
    x_positions = 1:length(salt_names_at_RH);
    offset = 0.15;
    
    scatter(x_positions - offset, x_water_mol, 80, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    scatter(x_positions + offset, x_water_ion, 80, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    plot([0 length(salt_names_at_RH)+1], [1 1], 'k--', 'LineWidth', 1.0);
    
    title(sprintf('RH = %.0f%%', target_RH*100), 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Salt', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel('Molar Water Uptake', 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_RH, 'XTickLabelRotation', 45)
    xlim([0 length(salt_names_at_RH)+1])
    ylim([0 1.05])
    set(gca, 'FontSize', 11)
end

sgtitle('Molar Water Uptake: Molecular vs Ionic Basis - Faceted by RH', 'FontSize', 20, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

subplot(2, 3, 6);
axis off; hold on;
h1 = scatter(1, 1.2, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
h2 = scatter(1, 0.8, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
h3 = plot([0.5 1.5], [0.4 0.4], 'k--', 'LineWidth', 1.5);
legend([h1, h2, h3], {'Molecular', 'Ionic', 'Ideal (x_w=1)'}, 'Location', 'best', 'FontSize', 14)
xlim([0 2]); ylim([0 2]);

filename = 'Molar_Water_Uptake_Comparison_Faceted_by_RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

%% PLOT 4: Mass Uptake - Molecular vs Ionic

figure('Position', [50, 50, 2000, 1200]);

for rh_idx = 1:num_RH_panels
    subplot(2, 3, rh_idx);
    hold on; grid on; box on;
    
    target_RH = data_at_RH(rh_idx).target_RH;
    mf_water_mol = data_at_RH(rh_idx).mf_water_mol;
    mf_water_ion = data_at_RH(rh_idx).mf_water_ion;
    salt_names_at_RH = data_at_RH(rh_idx).salt_names;
    
    if isempty(salt_names_at_RH)
        continue;
    end
    
    x_positions = 1:length(salt_names_at_RH);
    offset = 0.15;
    
    scatter(x_positions - offset, mf_water_mol, 80, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'Marker', 's');
    scatter(x_positions + offset, mf_water_ion, 80, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'Marker', 's');
    
    title(sprintf('RH = %.0f%%', target_RH*100), 'FontSize', 16, 'FontWeight', 'bold')
    xlabel('Salt', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel('Mass-Based Uptake', 'FontSize', 14, 'FontWeight', 'bold')
    set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_RH, 'XTickLabelRotation', 45)
    xlim([0 length(salt_names_at_RH)+1])
    ylim([0 1.05])
    set(gca, 'FontSize', 11)
end

sgtitle('Mass-Based Water Uptake: Molecular vs Ionic Basis - Faceted by RH', 'FontSize', 20, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

subplot(2, 3, 6);
axis off; hold on;
h1 = scatter(1, 1, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's');
h2 = scatter(1, 0.5, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's');
legend([h1, h2], {'Molecular', 'Ionic'}, 'Location', 'best', 'FontSize', 14)
xlim([0 2]); ylim([0 2]);

filename = 'Mass_Water_Uptake_Comparison_Faceted_by_RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

disp(['All faceted plots generated in: ' fig_out_dir])
