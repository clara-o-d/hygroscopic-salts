close all
clear
clc

% Script: licl_temp_RH.m
% 3D plot of LiCl mass-based water uptake as a function of RH and temperature
% Uses temperature-dependent Pitzer parameters from pitzer_binary.csv

% Add necessary paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% 1. Load Pitzer Parameters from CSV
fprintf('Loading Pitzer parameters for LiCl...\n');
csv_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_binary.csv');
pitzer_data = readtable(csv_file);

% Find LiCl parameters (Cl- and Li+)
licl_idx = find(strcmp(pitzer_data.species1, 'Cl-') & strcmp(pitzer_data.species2, 'Li+'));
if isempty(licl_idx)
    % Try reverse order
    licl_idx = find(strcmp(pitzer_data.species1, 'Li+') & strcmp(pitzer_data.species2, 'Cl-'));
end

if isempty(licl_idx)
    error('LiCl parameters not found in pitzer_binary.csv');
end

% Extract temperature-dependent coefficients
params.beta0_a1 = pitzer_data.beta0_a1(licl_idx);
params.beta0_a2 = pitzer_data.beta0_a2(licl_idx);
params.beta0_a3 = pitzer_data.beta0_a3(licl_idx);
params.beta0_a4 = pitzer_data.beta0_a4(licl_idx);

params.beta1_a1 = pitzer_data.beta1_a1(licl_idx);
params.beta1_a2 = pitzer_data.beta1_a2(licl_idx);
params.beta1_a3 = pitzer_data.beta1_a3(licl_idx);
params.beta1_a4 = pitzer_data.beta1_a4(licl_idx);

params.beta2_a1 = pitzer_data.beta2_a1(licl_idx);
params.beta2_a2 = pitzer_data.beta2_a2(licl_idx);
params.beta2_a3 = pitzer_data.beta2_a3(licl_idx);
params.beta2_a4 = pitzer_data.beta2_a4(licl_idx);

params.cphi_a1 = pitzer_data.cphi_a1(licl_idx);
params.cphi_a2 = pitzer_data.cphi_a2(licl_idx);
params.cphi_a3 = pitzer_data.cphi_a3(licl_idx);
params.cphi_a4 = pitzer_data.cphi_a4(licl_idx);

params.alpha1 = pitzer_data.alpha1(licl_idx);
params.alpha2 = pitzer_data.alpha2(licl_idx);

fprintf('LiCl Pitzer parameters loaded successfully.\n');

%% 2. Define Ranges
% Temperature range (Celsius) - extended to cover cold desert nights
T_C_vec = linspace(0, 50, 60);  % 0°C to 50°C (60 points for better resolution)
% RH range - extended to cover very dry conditions
RH_vec = linspace(0.10, 0.95, 60);  % 10% to 95% RH (60 points for better resolution)

% Verify the ranges
fprintf('\nTemperature vector: min=%.2f°C, max=%.2f°C, n=%d points\n', ...
    min(T_C_vec), max(T_C_vec), length(T_C_vec));
fprintf('RH vector: min=%.2f%%, max=%.2f%%, n=%d points\n', ...
    min(RH_vec)*100, max(RH_vec)*100, length(RH_vec));

% Constants
MW_LiCl = 42.4;      % g/mol
MW_water = 18.015;   % g/mol
nu = 2;              % Number of ions (Li+ and Cl-)
Tr = 298.15;         % Reference temperature (K)

%% 3. Calculate Water Uptake for each (T, RH) point
fprintf('Calculating water uptake for %d temperature points and %d RH points...\n', ...
    length(T_C_vec), length(RH_vec));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(T_C_vec), max(T_C_vec));
fprintf('RH range: %.1f%% to %.1f%%\n', min(RH_vec)*100, max(RH_vec)*100);

% Initialize arrays
water_uptake = zeros(length(T_C_vec), length(RH_vec));  % kg water / (kg water + kg salt)
water_uptake_ratio = zeros(length(T_C_vec), length(RH_vec));  % kg water / kg salt
mole_fraction_water = zeros(length(T_C_vec), length(RH_vec));
molality = zeros(length(T_C_vec), length(RH_vec));

for i = 1:length(T_C_vec)
    T_C = T_C_vec(i);
    T_K = T_C + 273.15;  % Convert to Kelvin
    
    if mod(i, 12) == 0 || i == 1
        fprintf('  Processing temperature %.1f°C (%d/%d)...\n', T_C, i, length(T_C_vec));
    end
    
    % Calculate temperature-dependent Pitzer parameters
    beta0 = calc_temp_param(params.beta0_a1, params.beta0_a2, params.beta0_a3, params.beta0_a4, T_K, Tr);
    beta1 = calc_temp_param(params.beta1_a1, params.beta1_a2, params.beta1_a3, params.beta1_a4, T_K, Tr);
    beta2 = calc_temp_param(params.beta2_a1, params.beta2_a2, params.beta2_a3, params.beta2_a4, T_K, Tr);
    cphi  = calc_temp_param(params.cphi_a1, params.cphi_a2, params.cphi_a3, params.cphi_a4, T_K, Tr);
    
    for j = 1:length(RH_vec)
        RH = RH_vec(j);
        
        % For each RH and T, we need to find the molality (m) such that
        % the water activity calculated from Pitzer model equals RH
        
        % Better initial guess based on empirical relationship
        % At low RH, LiCl solutions are very concentrated (high molality)
        % At high RH, solutions are dilute (low molality)
        if RH < 0.2
            m_guess = 15;  % Very concentrated
        elseif RH < 0.4
            m_guess = 8;
        elseif RH < 0.6
            m_guess = 4;
        elseif RH < 0.8
            m_guess = 2;
        else
            m_guess = 0.5;  % Dilute
        end
        
        % Solve for molality using fzero with robust bracketing
        m_solution = NaN;
        
        % Try bracket search first (most robust)
        try
            % Wide range to ensure we bracket the solution
            m_solution = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, ...
                params.alpha1, params.alpha2, nu, T_C), [0.01, 25], ...
                optimset('Display', 'off', 'TolX', 1e-5));
            
            if m_solution < 0
                m_solution = NaN;
            end
        catch
            % If bracketing fails, try from initial guess
            try
                m_solution = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, ...
                    params.alpha1, params.alpha2, nu, T_C), m_guess, ...
                    optimset('Display', 'off', 'TolX', 1e-5));
                
                if m_solution < 0 || m_solution > 30
                    m_solution = NaN;
                end
            catch
                % Last resort: try multiple starting points
                for m_try = [0.1, 1, 5, 10, 15, 20]
                    try
                        m_test = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, ...
                            params.alpha1, params.alpha2, nu, T_C), m_try, ...
                            optimset('Display', 'off', 'TolX', 1e-5));
                        if m_test > 0 && m_test < 30
                            m_solution = m_test;
                            break;
                        end
                    catch
                        continue;
                    end
                end
            end
        end
        
        molality(i, j) = m_solution;
        
        if ~isnan(m_solution) && m_solution > 0
            % Calculate mass fraction of salt
            % molality = moles_salt / kg_water
            % For 1 kg water: moles_salt = m_solution
            % mass_salt = m_solution * MW_LiCl (grams)
            mass_salt_g = m_solution * MW_LiCl;  % grams
            mass_water_g = 1000;  % grams (1 kg)
            
            mf_salt = mass_salt_g / (mass_salt_g + mass_water_g);
            mf_water = 1 - mf_salt;
            
            % Calculate mole fraction of water (ionic basis)
            % n_water = mass_water / MW_water
            % n_ions = nu * moles_salt = nu * m_solution (for 1 kg water)
            n_water = mass_water_g / MW_water;
            n_ions = nu * m_solution;
            
            x_water = n_water / (n_water + n_ions);
            mole_fraction_water(i, j) = x_water;
            
            % Water uptake = kg_water / (kg_water + kg_salt) = mass fraction of water
            water_uptake(i, j) = mf_water;  % dimensionless
            
            % Water uptake ratio = kg_water / kg_salt
            % For 1 kg water: mass_salt = m_solution * MW_LiCl / 1000 (kg)
            water_uptake_ratio(i, j) = 1.0 / (m_solution * MW_LiCl / 1000);  % kg H2O / kg LiCl
        else
            water_uptake(i, j) = NaN;
            water_uptake_ratio(i, j) = NaN;
            mole_fraction_water(i, j) = NaN;
        end
    end
end

fprintf('Calculation complete!\n');

% Report statistics on successful calculations
n_total = length(T_C_vec) * length(RH_vec);
n_valid = sum(~isnan(water_uptake(:)));
n_failed = sum(isnan(water_uptake(:)));
fprintf('Solution statistics: %d/%d points solved successfully (%.1f%%)\n', ...
    n_valid, n_total, n_valid/n_total*100);
if n_failed > 100
    fprintf('  Warning: %d points failed to converge (solver issues or below deliquescence)\n', n_failed);
end

%% 4. Load Atacama Desert Climate Data
fprintf('Loading Atacama desert climate data...\n');

% Load Atacama RH and Temperature data
atacama_rh_file = fullfile(filepath, '..', 'data', 'Atacama_RH.csv');
atacama_temp_file = fullfile(filepath, '..', 'data', 'Atacama_Temp.csv');

if exist(atacama_rh_file, 'file') && exist(atacama_temp_file, 'file')
    % Read the CSV files
    atacama_rh_data = readmatrix(atacama_rh_file);
    atacama_temp_data = readmatrix(atacama_temp_file);
    
    % Extract RH and Temperature values
    atacama_RH_full = atacama_rh_data(:, 2);  % Second column is RH (0-1)
    atacama_T_full = atacama_temp_data(:, 2);  % Second column is Temperature (°C)
    
    % Remove any NaN values from each array independently first
    atacama_RH_full = atacama_RH_full(~isnan(atacama_RH_full));
    atacama_T_full = atacama_T_full(~isnan(atacama_T_full));
    
    % Take the minimum length to ensure arrays are compatible
    min_length = min(length(atacama_RH_full), length(atacama_T_full));
    atacama_RH = atacama_RH_full(1:min_length);
    atacama_T = atacama_T_full(1:min_length);
    
    fprintf('  Loaded %d Atacama data points\n', length(atacama_RH));
    fprintf('  RH range: %.1f%% - %.1f%%\n', min(atacama_RH)*100, max(atacama_RH)*100);
    fprintf('  Temp range: %.1f°C - %.1f°C\n', min(atacama_T), max(atacama_T));
    
    % Interpolate water uptake for Atacama conditions
    atacama_water_uptake = zeros(size(atacama_RH));
    atacama_water_uptake_ratio = zeros(size(atacama_RH));
    atacama_molality = zeros(size(atacama_RH));
    atacama_mole_fraction = zeros(size(atacama_RH));
    
    fprintf('  Interpolating water uptake for Atacama conditions...\n');
    fprintf('  Grid coverage: RH [%.1f%%-%.1f%%], T [%.1f-%.1f°C]\n', ...
        min(RH_vec)*100, max(RH_vec)*100, min(T_C_vec), max(T_C_vec));
    
    for i = 1:length(atacama_RH)
        % Interpolate from our calculated grid
        % Use interp2 for 2D interpolation
        atacama_water_uptake(i) = interp2(RH_vec, T_C_vec, water_uptake, ...
            atacama_RH(i), atacama_T(i), 'linear', NaN);
        atacama_water_uptake_ratio(i) = interp2(RH_vec, T_C_vec, water_uptake_ratio, ...
            atacama_RH(i), atacama_T(i), 'linear', NaN);
        atacama_molality(i) = interp2(RH_vec, T_C_vec, molality, ...
            atacama_RH(i), atacama_T(i), 'linear', NaN);
        atacama_mole_fraction(i) = interp2(RH_vec, T_C_vec, mole_fraction_water, ...
            atacama_RH(i), atacama_T(i), 'linear', NaN);
    end
    
    % Remove any NaN values from interpolation
    valid_atacama = ~isnan(atacama_water_uptake);
    n_valid = sum(valid_atacama);
    n_invalid = sum(~valid_atacama);
    
    if n_invalid > 0
        fprintf('  Warning: %d points outside grid range (RH or T out of bounds)\n', n_invalid);
        % Show some examples of out-of-range points
        invalid_idx = find(~valid_atacama);
        if length(invalid_idx) <= 5
            for k = 1:length(invalid_idx)
                idx = invalid_idx(k);
                fprintf('    Point %d: RH=%.1f%%, T=%.1f°C (outside grid)\n', ...
                    idx, atacama_RH(idx)*100, atacama_T(idx));
            end
        else
            fprintf('    Example: RH=%.1f%%, T=%.1f°C\n', ...
                atacama_RH(invalid_idx(1))*100, atacama_T(invalid_idx(1)));
        end
    end
    
    atacama_RH = atacama_RH(valid_atacama);
    atacama_T = atacama_T(valid_atacama);
    atacama_water_uptake = atacama_water_uptake(valid_atacama);
    atacama_water_uptake_ratio = atacama_water_uptake_ratio(valid_atacama);
    atacama_molality = atacama_molality(valid_atacama);
    atacama_mole_fraction = atacama_mole_fraction(valid_atacama);
    
    fprintf('  Successfully interpolated %d / %d Atacama points (%.1f%%)\n', ...
        n_valid, length(atacama_RH_full), n_valid/length(atacama_RH_full)*100);
    atacama_data_available = true;
else
    fprintf('  Warning: Atacama data files not found, skipping overlay\n');
    atacama_data_available = false;
end

%% 5. Create 3D Surface Plot - Water Uptake
fprintf('Generating 3D plots...\n');

% Create meshgrid for plotting
[RH_grid, T_grid] = meshgrid(RH_vec * 100, T_C_vec);  % RH in percent

% Figure 1: Water Uptake (kg water / kg salt)
figure('Position', [100, 100, 1000, 700]);
surf(RH_grid, T_grid, water_uptake, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

% Add Atacama desert data line
if atacama_data_available
    plot3(atacama_RH * 100, atacama_T, atacama_water_uptake, ...
        'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot3(atacama_RH * 100, atacama_T, atacama_water_uptake, ...
        'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    legend('Location', 'best', 'FontSize', 10);
end

colormap(jet);
colorbar;
xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Water Uptake: kg H_2O / (kg H_2O + kg LiCl)', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl Water Uptake vs RH and Temperature (Pitzer Model)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
view(-30, 30);
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_WaterUptake_3D'), '-dpng', '-r300');

% Figure 2: Contour Plot of Water Uptake
figure('Position', [150, 150, 900, 700]);
contourf(RH_grid, T_grid, water_uptake, 20, 'LineColor', 'none');
hold on;
colormap(jet);
colorbar;
% Add contour lines
contour(RH_grid, T_grid, water_uptake, 10, 'LineColor', 'k', 'LineWidth', 0.5);

% Add Atacama desert data line
if atacama_data_available
    plot(atacama_RH * 100, atacama_T, 'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot(atacama_RH * 100, atacama_T, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    % Add arrow to show direction of daily cycle
    n_pts = length(atacama_RH);
    if n_pts > 10
        arrow_idx = round(n_pts/4);  % Arrow at 1/4 of the cycle
        quiver(atacama_RH(arrow_idx)*100, atacama_T(arrow_idx), ...
            (atacama_RH(arrow_idx+1) - atacama_RH(arrow_idx))*100, ...
            (atacama_T(arrow_idx+1) - atacama_T(arrow_idx)), ...
            0.5, 'r', 'LineWidth', 2, 'MaxHeadSize', 1.5, 'HandleVisibility', 'off');
    end
    legend('Location', 'best', 'FontSize', 10);
end

xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl Water Uptake: kg H_2O / (kg H_2O + kg LiCl) - Contour Plot', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_WaterUptake_Contour'), '-dpng', '-r300');

% Figure 3: Water Uptake Ratio (kg H2O / kg LiCl)
figure('Position', [175, 175, 1000, 700]);
surf(RH_grid, T_grid, water_uptake_ratio, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

% Add Atacama desert data line
if atacama_data_available
    plot3(atacama_RH * 100, atacama_T, atacama_water_uptake_ratio, ...
        'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot3(atacama_RH * 100, atacama_T, atacama_water_uptake_ratio, ...
        'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    legend('Location', 'best', 'FontSize', 10);
end

colormap(jet);
colorbar;
xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Water Uptake: kg H_2O / kg LiCl', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl Water Uptake Ratio vs RH and Temperature (Pitzer Model)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
view(-30, 30);
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_WaterUptakeRatio_3D'), '-dpng', '-r300');

% Figure 4: Contour Plot of Water Uptake Ratio
figure('Position', [225, 225, 900, 700]);
contourf(RH_grid, T_grid, water_uptake_ratio, 20, 'LineColor', 'none');
hold on;
colormap(jet);
colorbar;
% Add contour lines
contour(RH_grid, T_grid, water_uptake_ratio, 10, 'LineColor', 'k', 'LineWidth', 0.5);

% Add Atacama desert data line
if atacama_data_available
    plot(atacama_RH * 100, atacama_T, 'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot(atacama_RH * 100, atacama_T, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    % Add arrow to show direction of daily cycle
    n_pts = length(atacama_RH);
    if n_pts > 10
        arrow_idx = round(n_pts/4);  % Arrow at 1/4 of the cycle
        quiver(atacama_RH(arrow_idx)*100, atacama_T(arrow_idx), ...
            (atacama_RH(arrow_idx+1) - atacama_RH(arrow_idx))*100, ...
            (atacama_T(arrow_idx+1) - atacama_T(arrow_idx)), ...
            0.5, 'r', 'LineWidth', 2, 'MaxHeadSize', 1.5, 'HandleVisibility', 'off');
    end
    legend('Location', 'best', 'FontSize', 10);
end

xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl Water Uptake Ratio: kg H_2O / kg LiCl - Contour Plot', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_WaterUptakeRatio_Contour'), '-dpng', '-r300');

% Figure 5: Mole Fraction of Water
figure('Position', [275, 275, 1000, 700]);
surf(RH_grid, T_grid, mole_fraction_water, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

% Add Atacama desert data line
if atacama_data_available
    plot3(atacama_RH * 100, atacama_T, atacama_mole_fraction, ...
        'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot3(atacama_RH * 100, atacama_T, atacama_mole_fraction, ...
        'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    legend('Location', 'best', 'FontSize', 10);
end

colormap(jet);
colorbar;
xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Mole Fraction of Water (Ionic Basis)', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl: Water Mole Fraction vs RH and Temperature', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
view(-30, 30);
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_MoleFraction_3D'), '-dpng', '-r300');

% Figure 6: Molality
figure('Position', [325, 325, 1000, 700]);
surf(RH_grid, T_grid, molality, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

% Add Atacama desert data line
if atacama_data_available
    plot3(atacama_RH * 100, atacama_T, atacama_molality, ...
        'r-', 'LineWidth', 3, 'DisplayName', 'Atacama Desert Daily Cycle');
    plot3(atacama_RH * 100, atacama_T, atacama_molality, ...
        'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    legend('Location', 'best', 'FontSize', 10);
end

colormap(jet);
colorbar;
xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Molality (mol LiCl / kg H_2O)', 'FontSize', 12, 'FontWeight', 'bold');
title('LiCl Molality vs RH and Temperature', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
view(-30, 30);
set(gca, 'FontSize', 11);
print(fullfile(fig_out_dir, 'LiCl_Molality_3D'), '-dpng', '-r300');

if atacama_data_available
    fprintf('All plots saved to: %s\n', fig_out_dir);
    fprintf('Atacama desert daily cycle overlaid on all plots (red line)\n');
else
    fprintf('All plots saved to: %s\n', fig_out_dir);
end

%% Helper Functions

function param_val = calc_temp_param(a1, a2, a3, a4, T, Tr)
    % Calculate temperature-dependent parameter
    % Formula: P(T) = a1 + a2*(1/T - 1/Tr) + a3*ln(T/Tr) + a4*(T - Tr)
    % where T and Tr are in Kelvin
    
    % Handle NaN values
    if isnan(a1), a1 = 0; end
    if isnan(a2), a2 = 0; end
    if isnan(a3), a3 = 0; end
    if isnan(a4), a4 = 0; end
    
    param_val = a1 + a2*(1/T - 1/Tr) + a3*log(T/Tr) + a4*(T - Tr);
end

function residual = pitzer_residual(m, RH_target, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T)
    % Calculate residual: (calculated_water_activity - RH_target)
    % m = molality (mol/kg water)
    
    if m <= 0
        residual = 1e6;  % Penalty for negative molality
        return;
    end
    
    % Calculate water activity using Pitzer model
    aw_calc = pitzer_water_activity(m, nu, beta0, beta1, beta2, cphi, alpha1, alpha2, T);
    
    residual = aw_calc - RH_target;
end
