close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Define all salts with their molecular weights, valid RH ranges, and ION COUNTS
% Format: {salt_name, MW, RH_min, RH_max, function_name, function_args, num_cations, num_anions}
% function_args: 0 = no args, 1 = needs T

salt_data = {
    % Endothermic salts
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
    
    % Exothermic salts
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
    
    % Sulfates
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1};
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1};
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2};
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2};
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2};
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2};
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2};
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2};
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2};
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
};

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
    n_cat = salt_data{s}{7};
    n_an  = salt_data{s}{8};
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
