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

if isempty(salt_data)
    error('salt_data is empty!');
end

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    
    % Get Ion counts
    n_cat = salt_data{s}{12};
    n_an  = salt_data{s}{13};
    nu = n_cat + n_an; % Total number of dissociated ions
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    
    x_water_mol = zeros(size(RH_vec)); % Molecular basis
    x_water_ion = zeros(size(RH_vec)); % Ionic basis
    
    % Calculate for each RH value
    success = true;
    for i = 1:length(RH_vec)
        try
            if func_args == 1
                mf_salt(i) = feval(func_name, RH_vec(i), T);
            else
                mf_salt(i) = feval(func_name, RH_vec(i));
            end
            mf_water(i) = 1 - mf_salt(i);
            
            % Moles of water and salt (per unit mass basis)
            n_w = mf_water(i) / MWw;
            n_s = mf_salt(i) / MW;
            
            % 1. Molecular Definition: x_w = n_w / (n_w + n_s)
            x_water_mol(i) = n_w / (n_w + n_s);
            
            % 2. Ionic Definition: x_w = n_w / (n_w + nu * n_s)
            x_water_ion(i) = n_w / (n_w + (nu * n_s));
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        gamma_w_mol = RH_vec ./ x_water_mol;
        gamma_w_ion = RH_vec ./ x_water_ion;
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        
        all_salts_data.(salt_name).x_water_mol = x_water_mol;
        all_salts_data.(salt_name).gamma_w_mol = gamma_w_mol;
        
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).gamma_w_ion = gamma_w_ion;
        
        all_salts_data.(salt_name).mf_water = mf_water;
        all_salts_data.(salt_name).mf_salt = mf_salt;
        all_salts_data.(salt_name).MW = MW;
        all_salts_data.(salt_name).nu = nu;
        all_salts_data.(salt_name).display_name = salt_name;
        
        % Calculate mass-based water uptake from molar water uptake
        % Mass fraction of water: mf_w = (x_w * MWw) / (x_w * MWw + (1-x_w) * MW)
        % For molecular basis
        mf_water_from_x_mol = (x_water_mol * MWw) ./ (x_water_mol * MWw + (1 - x_water_mol) * MW);
        % For ionic basis - need to account for ion count
        % x_water_ion = n_w / (n_w + nu * n_s)
        % We can derive mf_water from x_water_ion:
        % Let r = n_w / n_s, then x_water_ion = r / (r + nu)
        % So r = nu * x_water_ion / (1 - x_water_ion)
        % mf_water = (r * MWw) / (r * MWw + MW)
        r_ion = (nu * x_water_ion) ./ (1 - x_water_ion);
        mf_water_from_x_ion = (r_ion * MWw) ./ (r_ion * MWw + MW);
        
        all_salts_data.(salt_name).mf_water_from_x_mol = mf_water_from_x_mol;
        all_salts_data.(salt_name).mf_water_from_x_ion = mf_water_from_x_ion;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

% --- Color Generation ---
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

%% --- PLOTTING LOGIC ---

% Definitions mapping
definitions = {'mol', 'ion'};

for d = 1:length(definitions)
    mode = definitions{d};
    
    % Determine Titles and Filenames based on mode
    if strcmp(mode, 'mol')
        % Molecular basis
        title_suffix = ' (Molecular Basis)'; 
        file_suffix  = '_Molecular';
        xlabel_str   = 'Relative Humidity (%)';
        ylabel_str   = 'Molar Water Uptake (x_w = RH/\gamma_w)';
    else
        % Ionic basis
        title_suffix = ' (Ionic Basis)';
        file_suffix  = '_Ionic';
        xlabel_str   = 'Relative Humidity (%)';
        ylabel_str   = 'Molar Water Uptake (x_w = RH/\gamma_w, Ionic)';
    end
    
    x_field = ['x_water_' mode];
    
    %% Plot Mole Fraction vs RH
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.RH*100, data.(x_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));

        x_coords = data.RH*100;
        y_coords = data.(x_field);
        text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    xlabel(xlabel_str, 'FontSize', 20, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 20, 'FontWeight', 'bold')
    title(['Molar Water Uptake' title_suffix ' vs Relative Humidity'], 'FontSize', 24, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 14, 'NumColumns', 3)
    xlim([0 100])
    ylim([0 1])
    set(gca, 'FontSize', 18)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('MoleFraction_Water%s_vs_RH', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    savefig(fullfile(fig_out_dir, filename))
    
end

%% --- DUAL Y-AXIS PLOTS: MOLAR AND MASS-BASED WATER UPTAKE ---

for d = 1:length(definitions)
    mode = definitions{d};
    
    % Determine Titles and Filenames based on mode
    if strcmp(mode, 'mol')
        % Molecular basis
        title_suffix = ' (Molecular Basis)'; 
        file_suffix  = '_Molecular';
        xlabel_str   = 'Relative Humidity (%)';
        mf_field = 'mf_water_from_x_mol';
    else
        % Ionic basis
        title_suffix = ' (Ionic Basis)';
        file_suffix  = '_Ionic';
        xlabel_str   = 'Relative Humidity (%)';
        mf_field = 'mf_water_from_x_ion';
    end
    
    x_field = ['x_water_' mode];
    
    %% Plot Molar and Mass-based Water Uptake vs RH (Dual Y-axis)
    figure('Position', [100, 100, 1400, 800]);
    
    % Left y-axis: Molar water uptake
    yyaxis left
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.RH*100, data.(x_field), 'LineWidth', 2.5, ...
             'DisplayName', [data.display_name ' (Molar)'], 'color', colors(s,:));
    end
    
    ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
    ylim([0 1])
    
    % Right y-axis: Mass-based water uptake
    yyaxis right
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.RH*100, data.(mf_field), 'LineWidth', 2.5, 'LineStyle', '--', ...
             'DisplayName', [data.display_name ' (Mass)'], 'color', colors(s,:));
    end
    
    ylabel('Mass-Based Water Uptake (Mass Fraction)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
    ylim([0 1])
    
    xlabel(xlabel_str, 'FontSize', 20, 'FontWeight', 'bold')
    title(['Molar and Mass-Based Water Uptake' title_suffix ' vs Relative Humidity'], 'FontSize', 24, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 12, 'NumColumns', 3)
    xlim([0 100])
    set(gca, 'FontSize', 18)
    set(gcf, 'color', 'w');
    
    % Set colors for y-axes
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    % Save
    filename = sprintf('Molar_and_Mass_Water_Uptake%s_vs_RH', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    savefig(fullfile(fig_out_dir, filename))
    
end

%% --- CATEGORICAL PLOT AT 80% RH ---

% Target RH for comparison
target_RH = 0.80;

% Extract data at 80% RH for each salt
x_water_mol_at_80 = [];
x_water_ion_at_80 = [];
mf_water_mol_at_80 = [];
mf_water_ion_at_80 = [];
salt_names_at_80 = {};

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    
    % Check if 80% RH is within the valid range for this salt
    if target_RH >= min(data.RH) && target_RH <= max(data.RH)
        % Find closest RH value to 0.80
        [~, idx] = min(abs(data.RH - target_RH));
        
        x_water_mol_at_80 = [x_water_mol_at_80; data.x_water_mol(idx)];
        x_water_ion_at_80 = [x_water_ion_at_80; data.x_water_ion(idx)];
        mf_water_mol_at_80 = [mf_water_mol_at_80; data.mf_water_from_x_mol(idx)];
        mf_water_ion_at_80 = [mf_water_ion_at_80; data.mf_water_from_x_ion(idx)];
        salt_names_at_80{end+1} = salt_name;
    end
end

fprintf('Found %d salts with data at 80%% RH\n', length(salt_names_at_80));

% Plot Molecular Basis at 80% RH
figure('Position', [100, 100, 1400, 700]);
hold on; grid on; box on;

x_positions = 1:length(salt_names_at_80);
scatter(x_positions, x_water_mol_at_80, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Add horizontal line at x=1 for reference (ideal solution)
plot([0 length(salt_names_at_80)+1], [1 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal (x_w = 1)');

xlabel('Salt', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar Water Uptake (Molecular Basis) at 80% RH', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_80, 'XTickLabelRotation', 45)
xlim([0 length(salt_names_at_80)+1])
ylim([0 1.05])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w');

% Save
filename = 'MoleFraction_Water_Molecular_at_80RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

% Plot Ionic Basis at 80% RH
figure('Position', [100, 100, 1400, 700]);
hold on; grid on; box on;

scatter(x_positions, x_water_ion_at_80, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Add horizontal line at x=1 for reference (ideal solution)
plot([0 length(salt_names_at_80)+1], [1 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal (x_w = 1)');

xlabel('Salt', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Molar Water Uptake (x_w = RH/\gamma_w, Ionic)', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar Water Uptake (Ionic Basis) at 80% RH', 'FontSize', 24, 'FontWeight', 'bold')
set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_80, 'XTickLabelRotation', 45)
xlim([0 length(salt_names_at_80)+1])
ylim([0 1.05])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w');

% Save
filename = 'MoleFraction_Water_Ionic_at_80RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

% Combined comparison plot (both bases side-by-side)
figure('Position', [100, 100, 1400, 700]);
hold on; grid on; box on;

% Offset positions for side-by-side comparison
offset = 0.15;
scatter(x_positions - offset, x_water_mol_at_80, 100, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Molecular Basis');
scatter(x_positions + offset, x_water_ion_at_80, 100, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Ionic Basis');

% Add horizontal line at x=1 for reference (ideal solution)
plot([0 length(salt_names_at_80)+1], [1 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal (x_w = 1)');

xlabel('Salt', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar Water Uptake at 80% RH (Comparison)', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 18)
set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_80, 'XTickLabelRotation', 45)
xlim([0 length(salt_names_at_80)+1])
ylim([0 1.05])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w');

% Save
filename = 'MoleFraction_Water_Comparison_at_80RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

%% --- DUAL Y-AXIS PLOT AT 80% RH: MOLAR AND MASS-BASED ---

% Plot Molecular Basis at 80% RH (Dual Y-axis)
figure('Position', [100, 100, 1400, 700]);
hold on; grid on; box on;

x_positions = 1:length(salt_names_at_80);

% Left y-axis: Molar water uptake
yyaxis left
scatter(x_positions, x_water_mol_at_80, 120, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Molar Water Uptake');
ylabel('Molar Water Uptake (x_w = RH/\gamma_w)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
ylim([0 1.05])

% Right y-axis: Mass-based water uptake
yyaxis right
scatter(x_positions, mf_water_mol_at_80, 120, 'filled', 'MarkerFaceColor', [0.8 0.4 0.2], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's', 'DisplayName', 'Mass-Based Water Uptake');
ylabel('Mass-Based Water Uptake (Mass Fraction)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
ylim([0 1.05])

xlabel('Salt', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar and Mass-Based Water Uptake (Molecular Basis) at 80% RH', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 18)
set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_80, 'XTickLabelRotation', 45)
xlim([0 length(salt_names_at_80)+1])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w');

% Set colors for y-axes
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% Save
filename = 'Molar_and_Mass_Water_Uptake_Molecular_at_80RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

% Plot Ionic Basis at 80% RH (Dual Y-axis)
figure('Position', [100, 100, 1400, 700]);
hold on; grid on; box on;

% Left y-axis: Molar water uptake
yyaxis left
scatter(x_positions, x_water_ion_at_80, 120, 'filled', 'MarkerFaceColor', [0.8 0.2 0.4], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Molar Water Uptake');
ylabel('Molar Water Uptake (x_w = RH/\gamma_w, Ionic)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
ylim([0 1.05])

% Right y-axis: Mass-based water uptake
yyaxis right
scatter(x_positions, mf_water_ion_at_80, 120, 'filled', 'MarkerFaceColor', [0.4 0.8 0.2], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'Marker', 's', 'DisplayName', 'Mass-Based Water Uptake');
ylabel('Mass-Based Water Uptake (Mass Fraction)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k')
ylim([0 1.05])

xlabel('Salt', 'FontSize', 20, 'FontWeight', 'bold')
title('Molar and Mass-Based Water Uptake (Ionic Basis) at 80% RH', 'FontSize', 24, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 18)
set(gca, 'XTick', x_positions, 'XTickLabel', salt_names_at_80, 'XTickLabelRotation', 45)
xlim([0 length(salt_names_at_80)+1])
set(gca, 'FontSize', 18)
set(gcf, 'color', 'w');

% Set colors for y-axes
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% Save
filename = 'Molar_and_Mass_Water_Uptake_Ionic_at_80RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
savefig(fullfile(fig_out_dir, filename))

disp(['All plots generated successfully in: ' fig_out_dir])
