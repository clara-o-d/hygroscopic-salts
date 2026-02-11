close all 
clear
clc 

% Add calculate_mf, util, and data folders to path
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
