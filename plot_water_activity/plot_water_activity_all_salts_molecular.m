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
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water_mol = zeros(size(RH_vec)); % Molecular basis
    
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
            
            % Molecular Definition: x_w = n_w / (n_w + n_s)
            % Each mole of salt counted as one mole regardless of dissociation
            x_water_mol(i) = n_w / (n_w + n_s);
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficient (gamma = a_w / x_w)
        % Using molecular basis: a_w = RH (water activity)
        gamma_w_mol = RH_vec ./ x_water_mol;
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        all_salts_data.(salt_name).x_water_mol = x_water_mol;
        all_salts_data.(salt_name).gamma_w_mol = gamma_w_mol;
        all_salts_data.(salt_name).mf_water = mf_water;
        all_salts_data.(salt_name).mf_salt = mf_salt;
        all_salts_data.(salt_name).MW = MW;
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

%% Plot Gamma vs Mole Fraction (Molecular Basis)
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.x_water_mol, data.gamma_w_mol, 'LineWidth', 2.5, ...
         'DisplayName', data.display_name, 'color', colors(s,:));

    x_coords = data.x_water_mol;
    y_coords = data.gamma_w_mol;
    text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
         'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
end

plot([0 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Mole Fraction of Water (x_w, Molecular Basis)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Mole Fraction (Molecular Basis)', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
xlim([0.3 1.0])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save
filename = 'Activity_Coefficient_Molecular_vs_MoleFraction';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')

%% Plot Gamma vs RH (Molecular Basis)
figure('Position', [150, 150, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    plot(data.RH*100, data.gamma_w_mol, 'LineWidth', 2.5, ...
         'DisplayName', data.display_name, 'color', colors(s,:));

    x_coords = data.RH*100;
    y_coords = data.gamma_w_mol;
    text(x_coords(1), y_coords(1), ['  ' data.display_name], ...
         'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
end

plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)')
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title('Water Activity Coefficient vs Relative Humidity (Molecular Basis)', 'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 8, 'NumColumns', 3)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save
filename = 'Activity_Coefficient_Molecular_vs_RH';
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')

%% --- Extra Plot: Gamma @ 90% vs Lower Bound (DRH Estimate) ---
% Target RH for comparison
target_RH = 0.90;

figure('Position', [200, 200, 900, 700]);
hold on; grid on; box on;

% Arrays to store valid points for correlation line (optional)
x_vals = [];
y_vals = [];

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    
    % 1. Get Lower Bound RH (The minimum RH used for this salt)
    rh_lower_bound = min(data.RH);
    
    % 2. Interpolate Gamma_w (Molecular) at exactly 90% RH
    % We use 'linear' interpolation. Extrapolation is turned off (returns NaN)
    % to avoid plotting salts that aren't dissolved at 90% (e.g. K2SO4 starts at 97%)
    gamma_at_90 = interp1(data.RH, data.gamma_w_mol, target_RH, 'linear', NaN);
    
    % 3. Plot only if valid (Salt is dissolved at 90% RH)
    if ~isnan(gamma_at_90)
        % Scatter point
        scatter(rh_lower_bound * 100, gamma_at_90, 80, colors(s,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'DisplayName', data.display_name);
        
        % Text Label (offset slightly for readability)
        text(rh_lower_bound * 100, gamma_at_90 + 0.002, ['  ' data.display_name], ...
             'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
             
        % Collect for simple stats/reference
        x_vals = [x_vals, rh_lower_bound * 100];
        y_vals = [y_vals, gamma_at_90];
    end
end

% Formatting
xlabel('Lower Bound RH / Approx. DRH (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(['Water Activity Coeff. (\gamma_{w,mol}) at ' num2str(target_RH*100) '% RH'], ...
    'FontSize', 12, 'FontWeight', 'bold');
title(['Screening: Non-Ideality at ' num2str(target_RH*100) '% RH vs. Deliquescence limit (Molecular Basis)'], ...
    'FontSize', 14, 'FontWeight', 'bold');

% Add a reference line for Ideal Behavior
yline(1.0, 'k--', 'Ideal (\gamma = 1)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

% Cosmetic adjustments
xlim([0 100]);
set(gca, 'FontSize', 12);
set(gcf, 'color', 'w');

% Save the extra figure
print(fullfile(fig_out_dir, 'Gamma90_vs_LowerBoundRH_Molecular'), '-dpng', '-r600');
disp('Extra scatter plot generated: Gamma90_vs_LowerBoundRH_Molecular');

disp('All plots generated successfully in: ../figures/activity_coefficient')
