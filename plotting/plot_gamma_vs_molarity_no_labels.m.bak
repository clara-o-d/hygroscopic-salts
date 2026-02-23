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
MWw_kg = 0.018015; % kg/mol

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
    molality = zeros(size(RH_vec));
    molarity = zeros(size(RH_vec));
    
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
            
            % Molality: m = moles_salt / kg_water
            molality(i) = (mf_salt(i) / MW) / (mf_water(i) / 1000);
            
            % Approximate molarity (assuming solution density ~ 1 g/mL)
            % More accurate: M = m * rho / (1 + m*MW/1000)
            % For simplicity, approximate density as 1 + 0.7*mf_salt (empirical)
            rho_approx = 1 + 0.7 * mf_salt(i); % g/mL
            molarity(i) = molality(i) * rho_approx / (1 + molality(i) * MW / 1000);
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        % Using ionic basis
        gamma_w_ion = RH_vec ./ x_water_ion;
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).gamma_w_ion = gamma_w_ion;
        all_salts_data.(salt_name).molality = molality;
        all_salts_data.(salt_name).molarity = molarity;
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
colors = generate_colors(num_salts);

%% --- PLOT: Activity Coefficient vs Molarity (NO LABELS) ---
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    
    % Plot gamma vs molarity (NO text labels)
    plot(data.molarity, data.gamma_w_ion, 'LineWidth', 2, ...
         'color', colors(s,:), 'HandleVisibility', 'off');
end

% Add ideal reference line
plot([0 max(cellfun(@(x) max(all_salts_data.(x).molarity), salt_names))], [1 1], ...
     'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');

xlabel('Molarity (M)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('\gamma_w (Water Activity Coefficient, Ionic Basis)', 'FontSize', 20, 'FontWeight', 'bold');
title('Water Activity Coefficient vs Molarity', 'FontSize', 22, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 16);
set(gcf, 'color', 'w');

% Set scaling factor for print
alpha = 0.7;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5*alpha 4*alpha];

% Save figure
print(fullfile(fig_out_dir, 'Activity_Coefficient_vs_Molarity_No_Labels'), '-dtiff', '-r600');
fprintf('Saved: %s\n', fullfile(fig_out_dir, 'Activity_Coefficient_vs_Molarity_No_Labels.tif'));

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
