close all
clear
clc

% Script: calculate_hcl_licl_atacama_swing.m
% Calculate water uptake swing for HCl + LiCl mixture at Atacama conditions
% This is the only mixture with complete Pitzer parameters that can deliquesce
% at Atacama RH levels

[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

%% 1. Load Atacama Climate Data
fprintf('=== HCl + LiCl MIXTURE ATACAMA SWING CALCULATION ===\n\n');

atacama_rh_file = fullfile(filepath, '..', 'data', 'Atacama_RH.csv');
atacama_temp_file = fullfile(filepath, '..', 'data', 'Atacama_Temp.csv');

atacama_rh_data = readmatrix(atacama_rh_file);
atacama_temp_data = readmatrix(atacama_temp_file);

atacama_RH_full = atacama_rh_data(:, 2);
atacama_T_full = atacama_temp_data(:, 2);

% Clean data
atacama_RH_full = atacama_RH_full(~isnan(atacama_RH_full));
atacama_T_full = atacama_T_full(~isnan(atacama_T_full));

min_length = min(length(atacama_RH_full), length(atacama_T_full));
atacama_RH = atacama_RH_full(1:min_length);
atacama_T = atacama_T_full(1:min_length);

fprintf('Loaded %d Atacama data points\n', length(atacama_RH));
fprintf('RH range: %.1f%% - %.1f%%\n', min(atacama_RH)*100, max(atacama_RH)*100);
fprintf('Temp range: %.1f°C - %.1f°C\n\n', min(atacama_T), max(atacama_T));

%% 2. Load Pitzer Parameters
fprintf('Loading Pitzer parameters...\n');
pitzer_dir = fullfile(filepath, '..', 'data', 'parsed_thermodb');
data_dir = fullfile(filepath, '..', 'data');

binary_data = readtable(fullfile(pitzer_dir, 'pitzer_binary.csv'));
theta_data = readtable(fullfile(data_dir, 'pitzer_theta_updated.csv'));
psi_data = readtable(fullfile(data_dir, 'pitzer_psi_updated.csv'));

% Load additional parameters from Lassin 2018 to fill gaps
theta_lassin = readtable(fullfile(data_dir, 'pitzer_theta_lassin2018.csv'));
psi_lassin = readtable(fullfile(data_dir, 'pitzer_psi_lassin2018.csv'));

fprintf('  Binary: %d, Theta: %d + %d (Lassin), Psi: %d + %d (Lassin)\n', ...
    height(binary_data), height(theta_data), height(theta_lassin), ...
    height(psi_data), height(psi_lassin));

params = struct();
params.binary = create_binary_map(binary_data);
params.theta = merge_parameter_maps(create_theta_map(theta_data), create_theta_map(theta_lassin));
params.psi = merge_parameter_maps(create_psi_map(psi_data), create_psi_map(psi_lassin));
fprintf('Parameters loaded and merged successfully.\n\n');

%% 3. Define Mixture
% HCl + LiCl (50:50 mass ratio)
MW_HCl = 36.461;
MW_LiCl = 42.394;
MW_water = 18.015;

system = struct(...
    'name', 'HCl_LiCl', ...
    'mass_fractions', [0.5, 0.5], ...
    'salt_MWs', [MW_HCl, MW_LiCl], ...
    'salt_ions', {{containers.Map({'H+', 'Cl-'}, {1, 1}), ...
                   containers.Map({'Li+', 'Cl-'}, {1, 1})}});

fprintf('Mixture: HCl + LiCl (50:50 mass ratio)\n');
fprintf('DRH: HCl=17%%, LiCl=12%%\n');
fprintf('Theta(H+-Li+) = 0.0150\n');
fprintf('Psi(H+-Li+-Cl-) = 0.0000\n\n');

%% 4. Calculate Water Uptake for HCl + LiCl Mixture
fprintf('Calculating water uptake for HCl + LiCl mixture...\n');

water_uptake_mixture = nan(size(atacama_RH));
molality_mixture = nan(size(atacama_RH));

options = optimset('Display', 'off', 'TolX', 1e-6, 'MaxIter', 100);

n_success_mix = 0;
n_failed_mix = 0;

for i = 1:length(atacama_RH)
    RH = atacama_RH(i);
    T = atacama_T(i);
    
    % Skip if below DRH of either salt
    if RH < 0.12  % Below LiCl DRH
        n_failed_mix = n_failed_mix + 1;
        continue;
    end
    
    try
        % Solve for molality that gives this RH
        m_result = fzero(@(m) calculate_aw_error(m, RH, system, params, T, MW_water), ...
                         2.0, options);
        
        if m_result > 0 && m_result < 100
            molality_mixture(i) = m_result;
            
            % Calculate water mass per 100g salt
            total_salt_mass = 100;  % grams
            mass_water = 1000 * (total_salt_mass / sum(system.mass_fractions .* system.salt_MWs)) / m_result;
            
            % Water uptake = kg_water / (kg_water + kg_salt)
            water_uptake_mixture(i) = mass_water / (mass_water + total_salt_mass);
            
            n_success_mix = n_success_mix + 1;
        else
            n_failed_mix = n_failed_mix + 1;
        end
    catch
        n_failed_mix = n_failed_mix + 1;
    end
    
    if mod(i, 20) == 0
        fprintf('  Processed %d/%d points...\n', i, length(atacama_RH));
    end
end

fprintf('Successfully calculated %d/%d points (%.1f%%)\n', n_success_mix, length(atacama_RH), ...
    n_success_mix/length(atacama_RH)*100);
fprintf('Failed/below DRH: %d points\n\n', n_failed_mix);

%% 5. Calculate Water Uptake for Pure LiCl
fprintf('Calculating water uptake for pure LiCl (for comparison)...\n');

water_uptake_licl = nan(size(atacama_RH));
molality_licl = nan(size(atacama_RH));

% Add path to calculate_mf functions
addpath(fullfile(filepath, '..', 'calculate_mf'));

n_success_licl = 0;
n_failed_licl = 0;

for i = 1:length(atacama_RH)
    RH = atacama_RH(i);
    T = atacama_T(i);
    
    if RH < 0.12  % Below LiCl DRH
        n_failed_licl = n_failed_licl + 1;
        continue;
    end
    
    try
        % Use the analytical function for pure LiCl
        mf_salt = calculate_mf_LiCl(RH, T);
        
        if ~isnan(mf_salt) && mf_salt > 0 && mf_salt < 1
            % Water uptake = mf_water = 1 - mf_salt
            water_uptake_licl(i) = 1 - mf_salt;
            
            % Calculate molality from mass fraction
            % mf_salt = m_salt / (m_salt + m_water)
            % For 1 kg water: m_salt = mf_salt / (1 - mf_salt) kg
            % molality = moles_salt / kg_water = (m_salt / MW_salt) / 1
            molality_licl(i) = (mf_salt / (1 - mf_salt)) / (MW_LiCl / 1000);
            
            n_success_licl = n_success_licl + 1;
        else
            n_failed_licl = n_failed_licl + 1;
        end
    catch
        n_failed_licl = n_failed_licl + 1;
    end
end

fprintf('Successfully calculated %d/%d points (%.1f%%)\n', n_success_licl, length(atacama_RH), ...
    n_success_licl/length(atacama_RH)*100);
fprintf('Failed/below DRH: %d points\n\n', n_failed_licl);

%% 6. Calculate Swings and Compare
fprintf('=== WATER UPTAKE SWING COMPARISON ===\n\n');

% HCl + LiCl Mixture
valid_uptake_mix = water_uptake_mixture(~isnan(water_uptake_mixture));
if length(valid_uptake_mix) < 5
    error('Insufficient valid data points for mixture swing calculation');
end
uptake_min_mix = min(valid_uptake_mix);
uptake_max_mix = max(valid_uptake_mix);
swing_mix = uptake_max_mix - uptake_min_mix;
[~, idx_min_mix] = min(water_uptake_mixture);
[~, idx_max_mix] = max(water_uptake_mixture);

fprintf('HCl + LiCl Mixture (50:50 mass):\n');
fprintf('  Minimum: %.4f (%.2f%% water) at RH=%.1f%%, T=%.1f°C\n', ...
    uptake_min_mix, uptake_min_mix*100, atacama_RH(idx_min_mix)*100, atacama_T(idx_min_mix));
fprintf('  Maximum: %.4f (%.2f%% water) at RH=%.1f%%, T=%.1f°C\n', ...
    uptake_max_mix, uptake_max_mix*100, atacama_RH(idx_max_mix)*100, atacama_T(idx_max_mix));
fprintf('  *** SWING: %.4f (%.2f%% range) ***\n\n', swing_mix, swing_mix*100);

% Pure LiCl
valid_uptake_licl = water_uptake_licl(~isnan(water_uptake_licl));
if length(valid_uptake_licl) < 5
    error('Insufficient valid data points for LiCl swing calculation');
end
uptake_min_licl = min(valid_uptake_licl);
uptake_max_licl = max(valid_uptake_licl);
swing_licl = uptake_max_licl - uptake_min_licl;
[~, idx_min_licl] = min(water_uptake_licl);
[~, idx_max_licl] = max(water_uptake_licl);

fprintf('Pure LiCl:\n');
fprintf('  Minimum: %.4f (%.2f%% water) at RH=%.1f%%, T=%.1f°C\n', ...
    uptake_min_licl, uptake_min_licl*100, atacama_RH(idx_min_licl)*100, atacama_T(idx_min_licl));
fprintf('  Maximum: %.4f (%.2f%% water) at RH=%.1f%%, T=%.1f°C\n', ...
    uptake_max_licl, uptake_max_licl*100, atacama_RH(idx_max_licl)*100, atacama_T(idx_max_licl));
fprintf('  *** SWING: %.4f (%.2f%% range) ***\n\n', swing_licl, swing_licl*100);

% Comparison
improvement = ((swing_mix - swing_licl) / swing_licl) * 100;
if improvement > 0
    fprintf('==> Mixture has %.1f%% BETTER swing than pure LiCl\n\n', improvement);
else
    fprintf('==> Mixture has %.1f%% WORSE swing than pure LiCl\n\n', abs(improvement));
end

%% 7. Create Visualization
fprintf('Creating comparison visualization...\n');

valid_idx_mix = ~isnan(water_uptake_mixture);
valid_idx_licl = ~isnan(water_uptake_licl);

fig = figure('Position', [100, 100, 1600, 900]);

% Plot 1: Water uptake vs time - comparison
subplot(2, 3, 1);
hold on;
plot(find(valid_idx_mix), water_uptake_mixture(valid_idx_mix) * 100, 'r-', 'LineWidth', 2, 'DisplayName', 'HCl+LiCl Mix');
plot(find(valid_idx_licl), water_uptake_licl(valid_idx_licl) * 100, 'b-', 'LineWidth', 2, 'DisplayName', 'Pure LiCl');
xlabel('Atacama Data Point Index', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Water Uptake (% water by mass)', 'FontSize', 11, 'FontWeight', 'bold');
title('Water Uptake Over Time', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% Plot 2: Swing comparison bar chart
subplot(2, 3, 2);
bar([swing_mix * 100, swing_licl * 100]);
set(gca, 'XTickLabel', {'HCl+LiCl', 'Pure LiCl'});
ylabel('Water Uptake Swing (%)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Swing Comparison (%.1f%% improvement)', improvement), 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 3: Water uptake vs RH - both systems
subplot(2, 3, 3);
hold on;
scatter(atacama_RH(valid_idx_mix) * 100, water_uptake_mixture(valid_idx_mix) * 100, ...
    40, 'r', 'filled', 'DisplayName', 'HCl+LiCl Mix');
scatter(atacama_RH(valid_idx_licl) * 100, water_uptake_licl(valid_idx_licl) * 100, ...
    40, 'b', 'DisplayName', 'Pure LiCl');
xlabel('Relative Humidity (%)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Water Uptake (% water by mass)', 'FontSize', 11, 'FontWeight', 'bold');
title('Water Uptake vs RH', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% Plot 4: Water uptake vs RH - mixture colored by T
subplot(2, 3, 4);
scatter(atacama_RH(valid_idx_mix) * 100, water_uptake_mixture(valid_idx_mix) * 100, ...
    50, atacama_T(valid_idx_mix), 'filled');
hold on;
plot(atacama_RH(idx_min_mix)*100, uptake_min_mix*100, 'ro', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
plot(atacama_RH(idx_max_mix)*100, uptake_max_mix*100, 'go', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel('Relative Humidity (%)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Water Uptake (% water by mass)', 'FontSize', 11, 'FontWeight', 'bold');
title('HCl+LiCl Mixture (colored by T)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;
ylabel(colorbar, 'Temperature (°C)', 'FontSize', 10);
grid on;
hold off;

% Plot 5: Water uptake vs Temperature - mixture
subplot(2, 3, 5);
scatter(atacama_T(valid_idx_mix), water_uptake_mixture(valid_idx_mix) * 100, ...
    50, atacama_RH(valid_idx_mix) * 100, 'filled');
xlabel('Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Water Uptake (% water by mass)', 'FontSize', 11, 'FontWeight', 'bold');
title('HCl+LiCl Mixture (colored by RH)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;
ylabel(colorbar, 'Relative Humidity (%)', 'FontSize', 10);
grid on;

% Plot 6: Water uptake vs RH - pure LiCl colored by T
subplot(2, 3, 6);
scatter(atacama_RH(valid_idx_licl) * 100, water_uptake_licl(valid_idx_licl) * 100, ...
    50, atacama_T(valid_idx_licl), 'filled');
hold on;
plot(atacama_RH(idx_min_licl)*100, uptake_min_licl*100, 'ro', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
plot(atacama_RH(idx_max_licl)*100, uptake_max_licl*100, 'go', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel('Relative Humidity (%)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Water Uptake (% water by mass)', 'FontSize', 11, 'FontWeight', 'bold');
title('Pure LiCl (colored by T)', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;
ylabel(colorbar, 'Temperature (°C)', 'FontSize', 10);
grid on;
hold off;

% Save figure
fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end
saveas(fig, fullfile(fig_out_dir, 'hcl_licl_mixture_vs_pure_atacama.png'));
fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'hcl_licl_mixture_vs_pure_atacama.png'));

%% Helper Functions

function err = calculate_aw_error(m_total, RH_target, system, params, T, MWw)
    if m_total <= 0 || m_total > 100
        err = 1e10;
        return;
    end
    
    try
        composition = get_ion_composition_from_mass_basis(system, m_total);
        aw_calc = calculate_water_activity_pitzer(composition, T, params, MWw);
        err = aw_calc - RH_target;
    catch
        err = 1e10;
    end
end

function composition = get_ion_composition_from_mass_basis(system, total_molality)
    composition = containers.Map();
    
    for s = 1:length(system.salt_ions)
        mf_salt = system.mass_fractions(s);
        MW_salt = system.salt_MWs(s);
        ion_map = system.salt_ions{s};
        
        m_salt = total_molality * (mf_salt / MW_salt) / sum(system.mass_fractions ./ system.salt_MWs);
        
        ion_names = keys(ion_map);
        for i = 1:length(ion_names)
            ion = ion_names{i};
            stoich = ion_map(ion);
            
            if isKey(composition, ion)
                composition(ion) = composition(ion) + m_salt * stoich;
            else
                composition(ion) = m_salt * stoich;
            end
        end
    end
end

function aw = calculate_water_activity_pitzer(composition, T, params, MWw)
    ion_names = keys(composition);
    charges = get_ion_charges();
    
    cations = {};
    anions = {};
    for i = 1:length(ion_names)
        ion = ion_names{i};
        if charges(ion) > 0
            cations{end+1} = ion;
        else
            anions{end+1} = ion;
        end
    end
    
    I = 0;
    for i = 1:length(ion_names)
        m_i = composition(ion_names{i});
        z_i = charges(ion_names{i});
        I = I + 0.5 * m_i * z_i^2;
    end
    
    A_phi = 0.3915;
    
    phi = calculate_osmotic_coefficient(composition, I, A_phi, params, charges, cations, anions);
    
    sum_m = 0;
    for i = 1:length(ion_names)
        sum_m = sum_m + composition(ion_names{i});
    end
    
    ln_aw = -phi * sum_m * MWw / 1000;
    aw = exp(ln_aw);
    aw = max(0, min(1, aw));
end

function phi = calculate_osmotic_coefficient(comp, I, A_phi, params, charges, cations, anions)
    sqrt_I = sqrt(I);
    b = 1.2;
    f_phi = -A_phi * sqrt_I / (1 + b * sqrt_I);
    
    % Binary term
    term_ca = 0;
    for i = 1:length(cations)
        cat = cations{i};
        m_c = comp(cat);
        z_c = charges(cat);
        
        for j = 1:length(anions)
            an = anions{j};
            m_a = comp(an);
            z_a = charges(an);
            
            key = [cat '_' an];
            if isKey(params.binary, key)
                p = params.binary(key);
                
                if abs(z_c) == 2 && abs(z_a) == 2
                    alpha1 = 1.4; alpha2 = 12.0;
                elseif abs(z_c) >= 3 || abs(z_a) >= 3
                    alpha1 = 2.0; alpha2 = 12.0;
                else
                    alpha1 = 2.0; alpha2 = 12.0;
                end
                
                B_phi = p.beta0 + p.beta1 * exp(-alpha1*sqrt_I);
                if abs(p.beta2) > 1e-10
                    B_phi = B_phi + p.beta2 * exp(-alpha2*sqrt_I);
                end
                
                term_ca = term_ca + m_c * m_a * (B_phi + sqrt_I * p.cphi);
            end
        end
    end
    
    % Theta term
    term_theta = 0;
    for i = 1:(length(cations)-1)
        for j = (i+1):length(cations)
            cat1 = cations{i};
            cat2 = cations{j};
            m_c1 = comp(cat1);
            m_c2 = comp(cat2);
            
            key = [cat1 '_' cat2];
            if isKey(params.theta, key)
                theta = params.theta(key);
                term_theta = term_theta + m_c1 * m_c2 * theta;
            end
        end
    end
    
    for i = 1:(length(anions)-1)
        for j = (i+1):length(anions)
            an1 = anions{i};
            an2 = anions{j};
            m_a1 = comp(an1);
            m_a2 = comp(an2);
            
            key = [an1 '_' an2];
            if isKey(params.theta, key)
                theta = params.theta(key);
                term_theta = term_theta + m_a1 * m_a2 * theta;
            end
        end
    end
    
    % Psi term
    term_psi = 0;
    for i = 1:(length(cations)-1)
        for j = (i+1):length(cations)
            for k = 1:length(anions)
                cat1 = cations{i};
                cat2 = cations{j};
                an = anions{k};
                m_c1 = comp(cat1);
                m_c2 = comp(cat2);
                m_a = comp(an);
                
                key = [cat1 '_' cat2 '_' an];
                if isKey(params.psi, key)
                    psi = params.psi(key);
                    term_psi = term_psi + m_c1 * m_c2 * m_a * psi;
                end
            end
        end
    end
    
    for i = 1:length(cations)
        for j = 1:(length(anions)-1)
            for k = (j+1):length(anions)
                cat = cations{i};
                an1 = anions{j};
                an2 = anions{k};
                m_c = comp(cat);
                m_a1 = comp(an1);
                m_a2 = comp(an2);
                
                key = [cat '_' an1 '_' an2];
                if isKey(params.psi, key)
                    psi = params.psi(key);
                    term_psi = term_psi + m_c * m_a1 * m_a2 * psi;
                end
            end
        end
    end
    
    sum_m = 0;
    ion_names = keys(comp);
    for i = 1:length(ion_names)
        sum_m = sum_m + comp(ion_names{i});
    end
    
    phi = 1 + (2/sum_m) * (I * f_phi + term_ca + term_theta + term_psi);
end

function params_map = create_binary_map(table_data)
    params_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = table_data.species1{i};
        sp2 = table_data.species2{i};
        
        key1 = [sp1 '_' sp2];
        key2 = [sp2 '_' sp1];
        
        beta0 = table_data.beta0_a1(i);
        beta1 = table_data.beta1_a1(i);
        beta2 = table_data.beta2_a1(i);
        cphi = table_data.cphi_a1(i);
        
        if isnan(beta0), beta0 = 0; end
        if isnan(beta1), beta1 = 0; end
        if isnan(beta2), beta2 = 0; end
        if isnan(cphi), cphi = 0; end
        
        param_struct = struct('beta0', beta0, 'beta1', beta1, 'beta2', beta2, 'cphi', cphi);
        params_map(key1) = param_struct;
        params_map(key2) = param_struct;
    end
end

function merged_map = merge_parameter_maps(base_map, supplement_map)
    % Merge two parameter maps, with supplement_map adding/overwriting entries
    merged_map = base_map;
    
    if ~isempty(supplement_map)
        supp_keys = keys(supplement_map);
        for i = 1:length(supp_keys)
            key = supp_keys{i};
            merged_map(key) = supplement_map(key);
        end
    end
end

function theta_map = create_theta_map(table_data)
    theta_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = table_data.species1{i};
        sp2 = table_data.species2{i};
        
        key1 = [sp1 '_' sp2];
        key2 = [sp2 '_' sp1];
        
        theta = table_data.theta_a1(i);
        if ~isnan(theta)
            theta_map(key1) = theta;
            theta_map(key2) = theta;
        end
    end
end

function psi_map = create_psi_map(table_data)
    psi_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = table_data.species1{i};
        sp2 = table_data.species2{i};
        sp3 = table_data.species3{i};
        
        keys = {[sp1 '_' sp2 '_' sp3], [sp1 '_' sp3 '_' sp2], ...
                [sp2 '_' sp1 '_' sp3], [sp2 '_' sp3 '_' sp1], ...
                [sp3 '_' sp1 '_' sp2], [sp3 '_' sp2 '_' sp1]};
        
        psi = table_data.psi_a1(i);
        if ~isnan(psi)
            for k = 1:length(keys)
                psi_map(keys{k}) = psi;
            end
        end
    end
end

function charges = get_ion_charges()
    charges = containers.Map();
    charges('H+') = 1;
    charges('Li+') = 1;
    charges('Na+') = 1;
    charges('K+') = 1;
    charges('Mg++') = 2;
    charges('Ca++') = 2;
    charges('Cl-') = -1;
    charges('Br-') = -1;
    charges('I-') = -1;
    charges('SO4--') = -2;
end
