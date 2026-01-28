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
    % {'NH4NO3', 80.043, 0.118, 0.732, 'calculate_mf_NH4NO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1}; % Adjusted Max (Limit 0.9936)
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};   % Adjusted Max (Limit 0.9290)
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};     % Adjusted Max (Limit 0.9669)
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};     % Adjusted Max (Limit 0.9528)
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};   % Adjusted Max (Limit 0.9527)
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};   % Adjusted Max (Limit 0.9482)
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};     % Adjusted Max (Limit 0.9624)
    
    % Exothermic salts (No source data provided to verify)
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
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1}; % Adjusted Boundaries
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};   % Adjusted Boundaries
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};% Adjusted Boundaries
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};   % Adjusted Boundaries
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};   % Adjusted Min (Limit 0.9190)
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1}; % Adjusted Boundaries
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};   % Adjusted Min (Limit 0.9710)
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};   % Adjusted Boundaries
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};   % Adjusted Min (Limit 0.9380)
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2}; % Adjusted Boundaries
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2}; % Adjusted Boundaries
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2}; % Adjusted Boundaries
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};   % Adjusted Boundaries
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2}; % Adjusted Boundaries
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2}; % Adjusted Boundaries
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};   % Adjusted Boundaries
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2}; % Adjusted Boundaries
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2}; % Adjusted Boundaries
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1}; % Adjusted Min (Limit 0.7775)
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
        ylabel_str   = 'ln(\gamma_w)';
    else
        % Explicit naming for the ionic definition
        title_suffix = ' (Ionic Basis)';
        file_suffix  = '_Ionic';
        xlabel_str   = 'Mole Fraction of Water (x_w, Ionic Basis)';
        ylabel_str   = 'ln(\gamma_w, Ionic)';
    end
    
    x_field = ['x_water_' mode];
    ln_gamma_field = ['ln_gamma_w_' mode];
    
    %% Plot ln(Gamma) vs Mole Fraction
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.(x_field), data.(ln_gamma_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));

        x_coords = data.(x_field);
        y_coords = data.(ln_gamma_field);
        text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    
    plot([0 1], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')
    xlabel(xlabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    title(['ln(\gamma_w)' title_suffix ' vs Mole Fraction'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    xlim([0.3 1.0])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('ln_Activity_Coefficient%s_vs_MoleFraction', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
    %% Plot ln(Gamma) vs RH
    figure('Position', [150, 150, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));

        x_coords = data.(x_field);
        y_coords = data.(ln_gamma_field);
        text(data.RH(1)*100, y_coords(1), ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    
    plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')
    xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    title(['ln(\gamma_w)' title_suffix ' vs Relative Humidity'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    xlim([0 100])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('ln_Activity_Coefficient%s_vs_RH', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
    %% Plot ln(Gamma) vs ln(aw)
    figure('Position', [200, 200, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.ln_aw, data.(ln_gamma_field), 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));

        x_coords = data.ln_aw;
        y_coords = data.(ln_gamma_field);
        text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    
    plot([-1 0], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')
    xlabel('ln(a_w)', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel(ylabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    title(['ln(\gamma_w)' title_suffix ' vs ln(a_w)'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('ln_Activity_Coefficient%s_vs_ln_aw', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
    %% Plot ln(aw) vs Mole Fraction
    figure('Position', [250, 250, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        plot(data.(x_field), data.ln_aw, 'LineWidth', 2.5, ...
             'DisplayName', data.display_name, 'color', colors(s,:));

        x_coords = data.(x_field);
        y_coords = data.ln_aw;
        text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    
    xlabel(xlabel_str, 'FontSize', 14, 'FontWeight', 'bold')
    ylabel('ln(a_w)', 'FontSize', 14, 'FontWeight', 'bold')
    title(['ln(a_w)' title_suffix ' vs Mole Fraction'], 'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
    xlim([0.3 1.0])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    % Save
    filename = sprintf('ln_aw%s_vs_MoleFraction', file_suffix);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    
end

disp('All plots generated successfully in: ../figures/activity_coefficient')
