close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'mole_fraction_water');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Load canonical salt data from data/load_salt_data.m
salt_data = load_salt_data();
exclude = {'NH4NO3', 'MgNO32'};
keep = cellfun(@(r) ~any(strcmp(r{1}, exclude)), salt_data);
salt_data = salt_data(keep);

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
    n_cat = salt_data{s}{12};
    n_an  = salt_data{s}{13};
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

%% UPTAKE VS RH: MOLAR AND MASS BASIS (All Salts) - WITH LABELS

figure('Position', [50, 50, 2400, 900]);

% Subplot 1: Molar Basis (Molecular) - with labels
subplot(1, 2, 1);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.x_water_mol, 'LineWidth', 2.5, 'color', colors(s,:));
    
    % Add text label at the beginning of each curve
    x_coords = data.RH*100;
    y_coords = data.x_water_mol;
    text(x_coords(1), y_coords(1), ['  ' salt_name], ...
         'Color', colors(s,:), 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left');
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar Basis (Molecular)', 'FontSize', 24, 'FontWeight', 'bold')
xlim([0 100])
ylim([0 1])
set(gca, 'FontSize', 18)

% Subplot 2: Mass Basis (from Molecular) - with labels
subplot(1, 2, 2);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.mf_water_from_x_mol, 'LineWidth', 2.5, 'color', colors(s,:));
    
    % Add text label at the beginning of each curve
    x_coords = data.RH*100;
    y_coords = data.mf_water_from_x_mol;
    text(x_coords(1), y_coords(1), ['  ' salt_name], ...
         'Color', colors(s,:), 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left');
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Mass-Based Water Uptake (Mass Fraction)', 'FontSize', 20, 'FontWeight', 'bold')
title('Mass Basis (from Molecular)', 'FontSize', 24, 'FontWeight', 'bold')
xlim([0 100])
ylim([0 1])
set(gca, 'FontSize', 18)

sgtitle('Water Uptake vs RH: Molar and Mass Basis for All Salts (Labeled)', 'FontSize', 26, 'FontWeight', 'bold')
set(gcf, 'color', 'w');

filename = 'Water_Uptake_vs_RH_Molar_and_Mass_All_Salts_Labeled';
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
