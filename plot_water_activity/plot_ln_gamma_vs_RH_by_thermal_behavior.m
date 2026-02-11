close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient');
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
        
        % Calculate water activity (aw = RH)
        aw = RH_vec;
        
        % Calculate logarithmic quantities
        ln_gamma_w_mol = log(gamma_w_mol);
        ln_gamma_w_ion = log(gamma_w_ion);
        ln_aw = log(aw);
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        
        all_salts_data.(salt_name).x_water_mol = x_water_mol;
        all_salts_data.(salt_name).gamma_w_mol = gamma_w_mol;
        all_salts_data.(salt_name).ln_gamma_w_mol = ln_gamma_w_mol;
        
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).gamma_w_ion = gamma_w_ion;
        all_salts_data.(salt_name).ln_gamma_w_ion = ln_gamma_w_ion;
        
        all_salts_data.(salt_name).aw = aw;
        all_salts_data.(salt_name).ln_aw = ln_aw;
        
        all_salts_data.(salt_name).mf_water = mf_water;
        all_salts_data.(salt_name).mf_salt = mf_salt;
        all_salts_data.(salt_name).MW = MW;
        all_salts_data.(salt_name).nu = nu;
        all_salts_data.(salt_name).display_name = salt_name;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

% --- Classification: Endothermic vs Exothermic ---
% Define endothermic salts (based on comments in salt_data and enthalpy data)
endothermic_salts = {'NaCl', 'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', ...
                     'LiNO3', 'KNO3', 'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', ...
                     'RbCl', 'CsBr', 'CsI', 'Na2SO4', 'K2SO4', 'NH42SO4', ...
                     'MgSO4', 'MnSO4', 'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4', ...
                     'BaNO3', 'BaCl2'};
% All other salts are exothermic

is_endothermic = false(num_salts, 1);
for s = 1:num_salts
    salt_name = salt_names{s};
    is_endothermic(s) = any(strcmp(salt_name, endothermic_salts));
end

%% --- Plot ln(Gamma) vs RH colored by Endothermic/Exothermic ---
% Define colors for endothermic and exothermic with transparency
% Using alpha = 0.5 for transparency (allows seeing through overlapping lines)
alpha = 0.5;
base_color_endothermic = [0.8500, 0.3250, 0.0980]; % Orange/red
base_color_exothermic = [0, 0.4470, 0.7410]; % Blue

% Count salts in each category
num_endothermic = sum(is_endothermic);
num_exothermic = sum(~is_endothermic);

% Plot for Molecular Basis
ln_gamma_field = 'ln_gamma_w_mol';

figure('Position', [300, 300, 1200, 800]);
hold on; grid on; box on;

% Plot endothermic salts with transparency
% Try to use Color property with alpha (works in MATLAB R2014b+)
for s = 1:num_salts
    if is_endothermic(s)
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2, ...
                 'color', base_color_endothermic, 'HandleVisibility', 'off');
        % Set transparency if supported (MATLAB R2014b+)
        try
            h.Color(4) = alpha;
        catch
            % Fallback: use lighter color if transparency not supported
            h.Color = base_color_endothermic * alpha + [1, 1, 1] * (1-alpha);
        end
    end
end

% Plot exothermic salts with transparency
for s = 1:num_salts
    if ~is_endothermic(s)
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2, ...
                 'color', base_color_exothermic, 'HandleVisibility', 'off');
        % Set transparency if supported (MATLAB R2014b+)
        try
            h.Color(4) = alpha;
        catch
            % Fallback: use lighter color if transparency not supported
            h.Color = base_color_exothermic * alpha + [1, 1, 1] * (1-alpha);
        end
    end
end

% Add reference line
plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')

% Add legend entries for endothermic/exothermic with counts
h_endothermic = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', base_color_endothermic, ...
                     'DisplayName', sprintf('Endothermic (n=%d)', num_endothermic));
h_exothermic = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', base_color_exothermic, ...
                     'DisplayName', sprintf('Exothermic (n=%d)', num_exothermic));

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('ln(\gamma_w) vs Relative Humidity - Molecular Basis (Colored by Endothermic/Exothermic)', ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 12)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save
filename = 'ln_Activity_Coefficient_vs_RH_by_Thermal_Behavior_Molecular';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')

disp('Molecular basis plot generated successfully!')
disp(['  - ' fullfile(fig_out_dir, [filename '.png'])])

% Plot for Ionic Basis
ln_gamma_field = 'ln_gamma_w_ion';

figure('Position', [300, 300, 1200, 800]);
hold on; grid on; box on;

% Plot endothermic salts with transparency
% Try to use Color property with alpha (works in MATLAB R2014b+)
for s = 1:num_salts
    if is_endothermic(s)
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2, ...
                 'color', base_color_endothermic, 'HandleVisibility', 'off');
        % Set transparency if supported (MATLAB R2014b+)
        try
            h.Color(4) = alpha;
        catch
            % Fallback: use lighter color if transparency not supported
            h.Color = base_color_endothermic * alpha + [1, 1, 1] * (1-alpha);
        end
    end
end

% Plot exothermic salts with transparency
for s = 1:num_salts
    if ~is_endothermic(s)
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2, ...
                 'color', base_color_exothermic, 'HandleVisibility', 'off');
        % Set transparency if supported (MATLAB R2014b+)
        try
            h.Color(4) = alpha;
        catch
            % Fallback: use lighter color if transparency not supported
            h.Color = base_color_exothermic * alpha + [1, 1, 1] * (1-alpha);
        end
    end
end

% Add reference line
plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')

% Add legend entries for endothermic/exothermic with counts
h_endothermic = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', base_color_endothermic, ...
                     'DisplayName', sprintf('Endothermic (n=%d)', num_endothermic));
h_exothermic = plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', base_color_exothermic, ...
                     'DisplayName', sprintf('Exothermic (n=%d)', num_exothermic));

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('ln(\gamma_w) vs Relative Humidity - Ionic Basis (Colored by Endothermic/Exothermic)', ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 12)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save
filename = 'ln_Activity_Coefficient_vs_RH_by_Thermal_Behavior_Ionic';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')

disp('Ionic basis plot generated successfully!')
disp(['  - ' fullfile(fig_out_dir, [filename '.png'])])
