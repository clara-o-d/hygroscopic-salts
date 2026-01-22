close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures');
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
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl_', 0, 1, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl_', 0, 1, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl_', 0, 1, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl_', 0, 1, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3_', 0, 1, 1};
    {'AgNO3', 169.87, 0.865, 0.99, 'calculate_mf_AgNO3_', 0, 1, 1};
    {'KI', 165.998, 0.975, 0.995, 'calculate_mf_KI_', 0, 1, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
    % {'NH4NO3', 80.043, 0.118, 0.732, 'calculate_mf_NH4NO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.995, 'calculate_mf_KClO3', 0, 1, 1};
    {'NaBr', 102.89, 0.614, 0.99, 'calculate_mf_NaBr', 0, 1, 1};
    {'NaI', 149.89, 0.581, 0.99, 'calculate_mf_NaI', 0, 1, 1};
    {'KBr', 119.00, 0.833, 0.99, 'calculate_mf_KBr', 0, 1, 1};
    {'RbCl', 120.92, 0.743, 0.99, 'calculate_mf_RbCl', 0, 1, 1};
    {'CsBr', 212.81, 0.848, 0.99, 'calculate_mf_CsBr', 0, 1, 1};
    {'CsI', 259.81, 0.913, 0.995, 'calculate_mf_CsI', 0, 1, 1};
    
    % Exothermic salts
    {'LiCl', 42.4, 0.12, 0.9, 'calculate_mf_LiCl', 1, 1, 1};
    {'LiOH', 24, 0.85, 0.9, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.9, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.9, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.9, 'calculate_mf_CaCl', 1, 1, 2};
    {'MgCl2', 95.2, 0.33, 0.9, 'calculate_mf_MgCl', 0, 1, 2};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
    {'LiBr', 86.85, 0.07, 0.9, 'calculate_mf_LiBr', 0, 1, 1};
    {'ZnCl2', 136.3, 0.07, 0.8, 'calculate_mf_ZnCl', 0, 1, 2};
    {'ZnI2', 319.18, 0.25, 0.9, 'calculate_mf_ZnI', 0, 1, 2};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
    {'LiI', 133.85, 0.18, 0.9, 'calculate_mf_LiI', 0, 1, 1};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.8990, 0.9957, 'calculate_mf_Na2SO4_', 0, 2, 1};
    {'K2SO4', 174.26, 0.9720, 0.9958, 'calculate_mf_K2SO4_', 0, 2, 1};
    {'NH42SO4', 132.14, 0.8310, 0.9959, 'calculate_mf_NH42SO4_', 0, 2, 1};
    {'MgSO4', 120.37, 0.9050, 0.9960, 'calculate_mf_MgSO4_', 0, 1, 1};
    {'MnSO4', 151.00, 0.8620, 0.9961, 'calculate_mf_MnSO4_', 0, 1, 1};
    {'Li2SO4', 109.94, 0.8530, 0.9956, 'calculate_mf_Li2SO4_', 0, 2, 1};
    {'NiSO4', 154.75, 0.9390, 0.9962, 'calculate_mf_NiSO4_', 0, 1, 1};
    {'CuSO4', 159.61, 0.9750, 0.9963, 'calculate_mf_CuSO4_', 0, 1, 1};
    {'ZnSO4', 161.44, 0.9130, 0.9962, 'calculate_mf_ZnSO4_', 0, 1, 1};
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9859, 0.9958, 'calculate_mf_BaNO3', 0, 1, 2};
    {'CaNO3', 164.09, 0.6464, 0.9955, 'calculate_mf_CaNO3', 0, 1, 2};
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6395, 0.9540, 'calculate_mf_CaBr2', 0, 1, 2};
    {'CaI2', 293.89, 0.8321, 0.9524, 'calculate_mf_CaI2', 0, 1, 2};
    {'SrCl2', 158.53, 0.8059, 0.9778, 'calculate_mf_SrCl2', 0, 1, 2};
    {'SrBr2', 247.43, 0.7776, 0.9571, 'calculate_mf_SrBr2', 0, 1, 2};
    {'SrI2', 341.43, 0.6785, 0.9569, 'calculate_mf_SrI2', 0, 1, 2};
    {'BaCl2', 208.23, 0.9375, 0.9731, 'calculate_mf_BaCl2', 0, 1, 2};
    {'BaBr2', 297.14, 0.8221, 0.9587, 'calculate_mf_BaBr2', 0, 1, 2};
    
    % Chlorates
    {'LiClO4', 106.39, 0.7775, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
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
        % ORIGINAL naming for the standard molecular definition
        title_suffix = ''; 
        file_suffix  = '';
        xlabel_str   = 'Mole Fraction of Water (x_w)';
        ylabel_str   = 'Water Activity Coefficient (\gamma_w)';
    else
        % Explicit naming for the ionic definition
        title_suffix = ' (Ionic Basis)';
        file_suffix  = '_Ionic';
        xlabel_str   = 'Mole Fraction of Water (x_w, Ionic Basis)';
        ylabel_str   = 'Water Activity Coefficient (\gamma_w, Ionic)';
    end
    
    x_field = ['x_water_' mode];
    gamma_field = ['gamma_w_' mode];
    
    %% Plot Gamma vs Mole Fraction
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.(x_field), data.(gamma_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));
    end
    
    plot([0 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
    xlabel(xlabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    title(['Water Activity Coefficient' title_suffix ' vs Mole Fraction'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    xlim([0.3 1.0])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('Activity_Coefficient%s_vs_MoleFraction', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
    %% Plot Gamma vs RH
    figure('Position', [150, 150, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.RH*100, data.(gamma_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));
    end
    
    plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
    xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    title(['Water Activity Coefficient' title_suffix ' vs Relative Humidity'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    xlim([0 100])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('Activity_Coefficient%s_vs_RH', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
end

disp('All plots generated successfully in: ../figures/')