close all 
clear
clc 

% Add necessary folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(filepath); % For PlotDefaults_Slides

% Apply plot defaults
PlotDefaults_Slides();

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015; % g/mol

% Load canonical salt data
salt_data = load_salt_data();
exclude = {'NH4NO3', 'MgNO32'}; % Exclude problematic salts
keep = cellfun(@(r) ~any(strcmp(r{1}, exclude)), salt_data);
salt_data = salt_data(keep);

% Process each salt
num_points = 100;
all_salts_data = struct();

fprintf('Processing %d salts...\n', length(salt_data));

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    
    % Get ion counts
    n_cat = salt_data{s}{12};
    n_an  = salt_data{s}{13};
    nu = n_cat + n_an; % Total number of dissociated ions
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water_ion = zeros(size(RH_vec));
    
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
            
            % Ionic mole fraction: x_w = n_w / (n_w + nu * n_s)
            x_water_ion(i) = n_w / (n_w + (nu * n_s));
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        % Using ionic basis: a_w = RH
        gamma_w_ion = RH_vec ./ x_water_ion;
        
        % Calculate logarithmic quantities
        ln_gamma_w_ion = log(gamma_w_ion);
        ln_x_water_ion = log(x_water_ion);
        ln_RH = log(RH_vec); % Note: RH is fractional (0-1)
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        all_salts_data.(salt_name).ln_RH = ln_RH;
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).gamma_w_ion = gamma_w_ion;
        all_salts_data.(salt_name).ln_gamma_w_ion = ln_gamma_w_ion;
        all_salts_data.(salt_name).ln_x_water_ion = ln_x_water_ion;
        all_salts_data.(salt_name).display_name = salt_name;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

% --- Color Generation ---
colors = generate_colors(num_salts);

%% --- PLOT: ln(gamma) and ln(x_w) vs ln(RH) (NO LABELS) ---
% This plot shows both ln(gamma_w) and ln(x_w) on the same axes
% The relationship: ln(a_w) = ln(gamma_w) + ln(x_w)
% So we expect: ln(gamma_w) + ln(x_w) = ln(RH)

figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    
    % Plot ln(gamma) vs ln(RH) - solid lines (NO text labels)
    h1 = plot(data.ln_RH, data.ln_gamma_w_ion, '-', 'LineWidth', 2, ...
         'color', colors(s,:), 'HandleVisibility', 'off');
    
    % Plot ln(x_w) vs ln(RH) - dashed lines with same color (NO text labels)
    h2 = plot(data.ln_RH, data.ln_x_water_ion, '--', 'LineWidth', 2, ...
         'color', colors(s,:), 'HandleVisibility', 'off');
end

% Add ideal reference line (ln(gamma) = 0, ln(x_w) = 0)
xlim_vals = xlim;
plot(xlim_vals, [0 0], 'k-', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0 or ln(x_w) = 0)');

% Add identity line ln(x_w) = ln(RH) for ideal solution
plot(xlim_vals, xlim_vals, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ideal (ln(x_w) = ln(RH))');

% Add legend entries for line styles
plot(NaN, NaN, 'k-', 'LineWidth', 2, 'DisplayName', 'ln(\gamma_w)');
plot(NaN, NaN, 'k--', 'LineWidth', 2, 'DisplayName', 'ln(x_w)');

xlabel('ln(RH) = ln(a_w)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('ln(\gamma_w) or ln(x_w) (Ionic Basis)', 'FontSize', 20, 'FontWeight', 'bold');
title('ln(\gamma_w) and ln(x_w) vs ln(RH): Chemical Potential Components', 'FontSize', 22, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 16);
set(gcf, 'color', 'w');

% Set scaling factor for print
alpha = 0.7;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5*alpha 4*alpha];

% Save figure
print(fullfile(fig_out_dir, 'ln_gamma_and_ln_xw_vs_ln_RH_No_Labels'), '-dtiff', '-r600');
fprintf('Saved: %s\n', fullfile(fig_out_dir, 'ln_gamma_and_ln_xw_vs_ln_RH_No_Labels.tif'));

% Additional annotation explaining the physics
annotation('textbox', [0.15, 0.75, 0.3, 0.15], ...
    'String', {'\mu_w - \mu_w^0 = RT ln(a_w)', ...
               '             = RT [ln(\gamma_w) + ln(x_w)]', ...
               '', ...
               'Solid: ln(\gamma_w) [non-ideality]', ...
               'Dashed: ln(x_w) [composition]'}, ...
    'FontSize', 14, 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
    'FitBoxToText', 'on');

%% Helper function for color generation
function colors = generate_colors(num_salts)
    colors = zeros(num_salts, 3);
    golden_angle = 0.38196601125; 
    for i = 1:num_salts
        hue = mod((i-1) * golden_angle, 1.0);
        sat = 0.85 + 0.15 * sin(i * 0.5); 
        val = 0.85 + 0.15 * cos(i * 0.4);
        colors(i,:) = hsv2rgb([hue, sat, val]);
    end
end
