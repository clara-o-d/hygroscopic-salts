close all
clear
clc

% Script to calculate and plot molar water uptake vs RH for LiCl 
% and LiCl combined with other salts using Pitzer equations

%% Setup paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

% Define output directory
fig_out_dir = fullfile(filepath, '..', 'figures', 'multi_salt');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% Load Pitzer parameters
pitzer_dir = fullfile(filepath, '..', 'data', 'parsed_thermodb');
binary_data = readtable(fullfile(pitzer_dir, 'pitzer_binary.csv'), 'TextType', 'string');
theta_data = readtable(fullfile(pitzer_dir, 'pitzer_theta.csv'), 'TextType', 'string');
psi_data = readtable(fullfile(pitzer_dir, 'pitzer_psi.csv'), 'TextType', 'string');

% Convert to maps for easier lookup
params = struct();
params.binary = create_binary_map(binary_data);
params.theta = create_theta_map(theta_data);
params.psi = create_psi_map(psi_data);

%% Constants
T = 25; % Temperature in Celsius
MWw = 18.015; % Molecular weight of water (g/mol)

%% Define salt systems to simulate
% Each system has: name, composition (molality ratios), MW, ions
salt_systems = {};

% === PURE SALTS (Reference systems) ===
% 1. Pure LiCl
ions_1 = containers.Map({'Li+', 'Cl-'}, {1, 1});
salt_systems{end+1} = struct(...
    'name', 'LiCl', ...
    'display_name', 'LiCl (pure)', ...
    'ions', ions_1, ...
    'MW_salt', 42.394, ...
    'color', [0.0, 0.0, 0.8], ...  % Dark blue
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_LiCl', ...
    'func_needs_T', true);

% Pure NaCl
ions_nacl = containers.Map({'Na+', 'Cl-'}, {1, 1});
salt_systems{end+1} = struct(...
    'name', 'NaCl', ...
    'display_name', 'NaCl (pure)', ...
    'ions', ions_nacl, ...
    'MW_salt', 58.443, ...
    'color', [0.8, 0.0, 0.0], ...  % Red
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_NaCl', ...
    'func_needs_T', false);

% Pure KCl
ions_kcl = containers.Map({'K+', 'Cl-'}, {1, 1});
salt_systems{end+1} = struct(...
    'name', 'KCl', ...
    'display_name', 'KCl (pure)', ...
    'ions', ions_kcl, ...
    'MW_salt', 74.551, ...
    'color', [0.0, 0.7, 0.0], ...  % Green
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_KCl', ...
    'func_needs_T', false);

% Pure MgCl2
ions_mgcl2 = containers.Map({'Mg++', 'Cl-'}, {1, 2});
salt_systems{end+1} = struct(...
    'name', 'MgCl2', ...
    'display_name', 'MgCl_2 (pure)', ...
    'ions', ions_mgcl2, ...
    'MW_salt', 95.211, ...
    'color', [0.8, 0.5, 0.0], ...  % Orange
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_MgCl', ...
    'func_needs_T', false);

% Pure CaCl2
ions_cacl2 = containers.Map({'Ca++', 'Cl-'}, {1, 2});
salt_systems{end+1} = struct(...
    'name', 'CaCl2', ...
    'display_name', 'CaCl_2 (pure)', ...
    'ions', ions_cacl2, ...
    'MW_salt', 110.984, ...
    'color', [0.6, 0.0, 0.8], ...  % Purple
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_CaCl', ...
    'func_needs_T', true);

% Pure LiBr
ions_libr = containers.Map({'Li+', 'Br-'}, {1, 1});
salt_systems{end+1} = struct(...
    'name', 'LiBr', ...
    'display_name', 'LiBr (pure)', ...
    'ions', ions_libr, ...
    'MW_salt', 86.845, ...
    'color', [0.0, 0.6, 0.6], ...  % Cyan
    'line_style', '--', ...
    'line_width', 3.5, ...
    'is_pure', true, ...
    'calc_func', 'calculate_mf_LiBr', ...
    'func_needs_T', false);

% === BINARY MIXTURES ===

% 2. LiCl + NaCl (equal molality)
ions_2 = containers.Map({'Li+', 'Na+', 'Cl-'}, {1, 1, 2});
salt_systems{end+1} = struct(...
    'name', 'LiCl_NaCl', ...
    'display_name', 'LiCl + NaCl', ...
    'ions', ions_2, ...
    'MW_salt', (42.394 + 58.443) / 2, ...
    'color', [0.8, 0.0, 0.4], ...
    'line_style', '-', ...
    'line_width', 2.5, ...
    'is_pure', false);

% 3. LiCl + KCl (equal molality)
ions_3 = containers.Map({'Li+', 'K+', 'Cl-'}, {1, 1, 2});
salt_systems{end+1} = struct(...
    'name', 'LiCl_KCl', ...
    'display_name', 'LiCl + KCl', ...
    'ions', ions_3, ...
    'MW_salt', (42.394 + 74.551) / 2, ...
    'color', [0.0, 0.8, 0.4], ...
    'line_style', '-', ...
    'line_width', 2.5, ...
    'is_pure', false);

% 4. LiCl + MgCl2
ions_4 = containers.Map({'Li+', 'Mg++', 'Cl-'}, {1, 1, 3});
salt_systems{end+1} = struct(...
    'name', 'LiCl_MgCl2', ...
    'display_name', 'LiCl + MgCl_2', ...
    'ions', ions_4, ...
    'MW_salt', (42.394 + 95.211) / 2, ...
    'color', [0.8, 0.4, 0.0], ...
    'line_style', '-', ...
    'line_width', 2.5, ...
    'is_pure', false);

% 5. LiCl + CaCl2
ions_5 = containers.Map({'Li+', 'Ca++', 'Cl-'}, {1, 1, 3});
salt_systems{end+1} = struct(...
    'name', 'LiCl_CaCl2', ...
    'display_name', 'LiCl + CaCl_2', ...
    'ions', ions_5, ...
    'MW_salt', (42.394 + 110.984) / 2, ...
    'color', [0.6, 0.0, 0.8], ...
    'line_style', '-', ...
    'line_width', 2.5, ...
    'is_pure', false);

% 6. LiCl + LiBr
ions_6 = containers.Map({'Li+', 'Cl-', 'Br-'}, {2, 1, 1});
salt_systems{end+1} = struct(...
    'name', 'LiCl_LiBr', ...
    'display_name', 'LiCl + LiBr', ...
    'ions', ions_6, ...
    'MW_salt', (42.394 + 86.845) / 2, ...
    'color', [0.4, 0.8, 0.0], ...
    'line_style', '-', ...
    'line_width', 2.5, ...
    'is_pure', false);

%% Calculate water uptake for each system
RH_vec = linspace(0.15, 0.95, 100);

results = struct();

fprintf('Calculating water uptake for salt systems...\n');

for sys_idx = 1:length(salt_systems)
    system = salt_systems{sys_idx};
    fprintf('  Processing: %s\n', system.display_name);
    
    % Initialize arrays
    x_water = zeros(size(RH_vec));
    molality_total = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    
    % Calculate for each RH
    for i = 1:length(RH_vec)
        RH = RH_vec(i);
        
        % For pure salts, use the analytical functions
        if system.is_pure
            try
                % Call the appropriate calculate_mf function
                if system.func_needs_T
                    mf_salt = feval(system.calc_func, RH, T);
                else
                    mf_salt = feval(system.calc_func, RH);
                end
                mf_water(i) = 1 - mf_salt;
                
                % Convert to molality
                m_salt = 1000 * mf_salt / (system.MW_salt * (1 - mf_salt));
                molality_total(i) = m_salt;
                
                % Calculate number of ions (based on composition)
                ion_names_sys = keys(system.ions);
                n_ions_per_formula = 0;
                for j = 1:length(ion_names_sys)
                    n_ions_per_formula = n_ions_per_formula + system.ions(ion_names_sys{j});
                end
                
                % Calculate mole fraction
                n_water = (1 - mf_salt) / MWw;
                n_ions = m_salt / 1000 * n_ions_per_formula;
                x_water(i) = n_water / (n_water + n_ions);
            catch ME
                % Some functions have restricted RH ranges
                x_water(i) = NaN;
                molality_total(i) = NaN;
                mf_water(i) = NaN;
            end
        else
            % For mixed salt systems, use Pitzer model iteratively
            try
                % Initial guess for total molality
                m_guess = 1.0;
                
                % Solve for molality that gives the target RH
                % RH = aw = exp(-phi * nu * m_total * Mw)
                % Use optimization
                options = optimset('Display', 'off', 'TolX', 1e-6);
                
                m_result = fzero(@(m) calculate_aw_error(m, RH, system, params, T, MWw), ...
                                 m_guess, options);
                
                if m_result > 0 && m_result < 100
                    molality_total(i) = m_result;
                    
                    % Calculate water activity and mole fraction
                    composition = scale_composition(system.ions, m_result);
                    aw_calc = calculate_water_activity_pitzer(composition, T, params, MWw);
                    
                    % Calculate total moles of ions
                    n_ions_total = 0;
                    ion_names = keys(composition);
                    for j = 1:length(ion_names)
                        n_ions_total = n_ions_total + composition(ion_names{j});
                    end
                    
                    % Mole fraction of water (ionic basis)
                    n_water = 55.51; % moles of water per kg
                    x_water(i) = n_water / (n_water + n_ions_total);
                    
                    % Mass fraction of water
                    % Total mass = 1000g water + m * MW_salt
                    mass_salt = m_result * system.MW_salt;
                    mf_water(i) = 1000 / (1000 + mass_salt);
                else
                    x_water(i) = NaN;
                    molality_total(i) = NaN;
                    mf_water(i) = NaN;
                end
            catch ME
                warning('Error for %s at RH=%.3f: %s', system.name, RH, ME.message);
                x_water(i) = NaN;
                molality_total(i) = NaN;
                mf_water(i) = NaN;
            end
        end
    end
    
    % Store results
    results.(system.name) = struct(...
        'RH', RH_vec, ...
        'x_water', x_water, ...
        'molality_total', molality_total, ...
        'mf_water', mf_water, ...
        'display_name', system.display_name, ...
        'color', system.color, ...
        'line_style', system.line_style, ...
        'line_width', system.line_width, ...
        'is_pure', system.is_pure);
end

fprintf('Calculations complete!\n');

%% Plotting

% Plot 1: Molar Water Uptake (x_water) vs RH
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

system_names = fieldnames(results);
for i = 1:length(system_names)
    data = results.(system_names{i});
    plot(data.RH * 100, data.x_water, ...
         'LineWidth', data.line_width, ...
         'LineStyle', data.line_style, ...
         'Color', data.color, ...
         'DisplayName', data.display_name);
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Molar Water Uptake (x_w, ionic basis)', 'FontSize', 20, 'FontWeight', 'bold');
title('Molar Water Uptake vs RH: LiCl and Multi-Salt Systems', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 14);
xlim([15 95]);
ylim([0 1]);
set(gca, 'FontSize', 18);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'binary_molar_water_uptake_vs_RH.png'));
savefig(fullfile(fig_out_dir, 'binary_molar_water_uptake_vs_RH.fig'));

% Plot 2: Mass-based Water Uptake vs RH
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for i = 1:length(system_names)
    data = results.(system_names{i});
    plot(data.RH * 100, data.mf_water, ...
         'LineWidth', data.line_width, ...
         'LineStyle', data.line_style, ...
         'Color', data.color, ...
         'DisplayName', data.display_name);
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Mass Fraction of Water', 'FontSize', 20, 'FontWeight', 'bold');
title('Mass-Based Water Uptake vs RH: LiCl and Multi-Salt Systems', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 16);
xlim([15 95]);
ylim([0 1]);
set(gca, 'FontSize', 18);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'binary_mass_water_uptake_vs_RH.png'));
savefig(fullfile(fig_out_dir, 'binary_mass_water_uptake_vs_RH.fig'));

% Plot 3: Total Molality vs RH
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for i = 1:length(system_names)
    data = results.(system_names{i});
    plot(data.RH * 100, data.molality_total, ...
         'LineWidth', data.line_width, ...
         'LineStyle', data.line_style, ...
         'Color', data.color, ...
         'DisplayName', data.display_name);
end

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Total Molality (mol/kg H_2O)', 'FontSize', 20, 'FontWeight', 'bold');
title('Total Salt Molality vs RH: LiCl and Multi-Salt Systems', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 16);
xlim([15 95]);
set(gca, 'FontSize', 18);
set(gcf, 'color', 'w');
set(gca, 'YScale', 'log');

% Save
saveas(gcf, fullfile(fig_out_dir, 'molality_vs_RH.png'));
savefig(fullfile(fig_out_dir, 'molality_vs_RH.fig'));

% Plot 4: Dual Y-axis plot (Molar and Mass-based)
figure('Position', [100, 100, 1400, 800]);

yyaxis left
for i = 1:length(system_names)
    data = results.(system_names{i});
    plot(data.RH * 100, data.x_water, ...
         'LineWidth', data.line_width, ...
         'LineStyle', data.line_style, ...
         'Color', data.color, ...
         'DisplayName', [data.display_name ' (Molar)']);
    hold on;
end
ylabel('Molar Water Uptake (x_w)', 'FontSize', 20, 'FontWeight', 'bold');
ylim([0 1]);
ax = gca;
ax.YColor = 'k';

yyaxis right
for i = 1:length(system_names)
    data = results.(system_names{i});
    plot(data.RH * 100, data.mf_water, ...
         'LineWidth', data.line_width, ...
         'LineStyle', ':', ...
         'Color', data.color, ...
         'DisplayName', [data.display_name ' (Mass)']);
end
ylabel('Mass Fraction of Water', 'FontSize', 20, 'FontWeight', 'bold');
ylim([0 1]);
ax.YColor = 'k';

xlabel('Relative Humidity (%)', 'FontSize', 20, 'FontWeight', 'bold');
title('Molar and Mass-Based Water Uptake vs RH', 'FontSize', 24, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12, 'NumColumns', 2);
xlim([15 95]);
grid on; box on;
set(gca, 'FontSize', 18);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'dual_water_uptake_vs_RH.png'));
savefig(fullfile(fig_out_dir, 'dual_water_uptake_vs_RH.fig'));

fprintf('\nAll plots saved to: %s\n', fig_out_dir);

%% Helper Functions

function params_map = create_binary_map(table_data)
    % Create a map for binary interaction parameters
    params_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = char(table_data.species1(i));
        sp2 = char(table_data.species2(i));
        
        % Skip empty or invalid rows
        if isempty(sp1) || isempty(sp2) || strcmp(sp1, 'x')
            continue;
        end
        
        % Create key (both orderings)
        key1 = [sp1 '_' sp2];
        key2 = [sp2 '_' sp1];
        
        % Extract parameters (handle NaN)
        beta0 = table_data.beta0_a1(i);
        beta1 = table_data.beta1_a1(i);
        beta2 = table_data.beta2_a1(i);
        cphi = table_data.cphi_a1(i);
        
        if isnan(beta0), beta0 = 0; end
        if isnan(beta1), beta1 = 0; end
        if isnan(beta2), beta2 = 0; end
        if isnan(cphi), cphi = 0; end
        
        param_struct = struct('beta0', beta0, 'beta1', beta1, ...
                             'beta2', beta2, 'cphi', cphi);
        
        params_map(key1) = param_struct;
        params_map(key2) = param_struct;
    end
end

function theta_map = create_theta_map(table_data)
    % Create a map for theta parameters (like-charge interactions)
    theta_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = char(table_data.species1(i));
        sp2 = char(table_data.species2(i));
        
        if isempty(sp1) || isempty(sp2)
            continue;
        end
        
        key1 = [sp1 '_' sp2];
        key2 = [sp2 '_' sp1];
        
        theta = table_data.theta_a1(i);
        if isnan(theta), theta = 0; end
        
        theta_map(key1) = theta;
        theta_map(key2) = theta;
    end
end

function psi_map = create_psi_map(table_data)
    % Create a map for psi parameters (ternary interactions)
    psi_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = char(table_data.species1(i));
        sp2 = char(table_data.species2(i));
        sp3 = char(table_data.species3(i));
        
        if isempty(sp1) || isempty(sp2) || isempty(sp3)
            continue;
        end
        
        % Create all possible orderings for lookup
        keys = {[sp1 '_' sp2 '_' sp3], [sp1 '_' sp3 '_' sp2], ...
                [sp2 '_' sp1 '_' sp3], [sp2 '_' sp3 '_' sp1], ...
                [sp3 '_' sp1 '_' sp2], [sp3 '_' sp2 '_' sp1]};
        
        psi = table_data.psi_a1(i);
        if isnan(psi), psi = 0; end
        
        for k = 1:length(keys)
            psi_map(keys{k}) = psi;
        end
    end
end

function comp = scale_composition(ions, total_molality)
    % Scale ion composition to match total molality
    % ions is a containers.Map
    ion_names = keys(ions);
    total_ratio = 0;
    for i = 1:length(ion_names)
        total_ratio = total_ratio + ions(ion_names{i});
    end
    
    comp = containers.Map();
    for i = 1:length(ion_names)
        comp(ion_names{i}) = ions(ion_names{i}) * total_molality / total_ratio;
    end
end

function aw = calculate_water_activity_pitzer(composition, T, params, MWw)
    % Calculate water activity using Pitzer model for multi-salt solutions
    % composition: containers.Map with ion names as keys and molalities as values
    
    % Get ion list and separate cations/anions
    ion_names = keys(composition);
    charges = get_ion_charges();
    
    cations = {};
    anions = {};
    for i = 1:length(ion_names)
        ion = ion_names{i};
        if ~isKey(charges, ion)
            error('Unknown ion: %s', ion);
        end
        if charges(ion) > 0
            cations{end+1} = ion;
        else
            anions{end+1} = ion;
        end
    end
    
    % Calculate ionic strength
    I = 0;
    for i = 1:length(ion_names)
        m_i = composition(ion_names{i});
        z_i = charges(ion_names{i});
        I = I + 0.5 * m_i * z_i^2;
    end
    
    % Debye-Huckel parameter
    A_phi = 0.3915; % at 25°C
    
    % Calculate osmotic coefficient
    phi = calculate_osmotic_coefficient(composition, I, A_phi, params, charges, cations, anions);
    
    % Calculate water activity
    % ln(aw) = -phi * sum(m_i) * MWw / 1000
    sum_m = 0;
    for i = 1:length(ion_names)
        sum_m = sum_m + composition(ion_names{i});
    end
    
    ln_aw = -phi * sum_m * MWw / 1000;
    aw = exp(ln_aw);
    
    % Bound between 0 and 1
    aw = max(0, min(1, aw));
end

function phi = calculate_osmotic_coefficient(comp, I, A_phi, params, charges, cations, anions)
    % Calculate osmotic coefficient using Pitzer model
    % comp is a containers.Map
    
    sqrt_I = sqrt(I);
    b = 1.2;
    
    % Debye-Huckel term
    f_phi = -A_phi * sqrt_I / (1 + b * sqrt_I);
    
    % Binary interaction terms (cation-anion)
    term_ca = 0;
    for i = 1:length(cations)
        cat = cations{i};
        m_c = comp(cat);
        z_c = charges(cat);
        
        for j = 1:length(anions)
            an = anions{j};
            m_a = comp(an);
            z_a = charges(an);
            
            % Get binary parameters
            key = [cat '_' an];
            if isKey(params.binary, key)
                p = params.binary(key);
                
                % Determine alpha values based on charge
                if abs(z_c) == 2 && abs(z_a) == 2
                    alpha1 = 1.4;
                    alpha2 = 12.0;
                elseif abs(z_c) >= 3 || abs(z_a) >= 3
                    alpha1 = 2.0;
                    alpha2 = 12.0;
                else
                    alpha1 = 2.0;
                    alpha2 = 12.0;
                end
                
                % B_phi calculation
                g1 = 2 * (1 - (1 + alpha1*sqrt_I) * exp(-alpha1*sqrt_I)) / (alpha1^2 * I);
                B_phi = p.beta0 + p.beta1 * exp(-alpha1*sqrt_I);
                
                if abs(p.beta2) > 1e-10
                    B_phi = B_phi + p.beta2 * exp(-alpha2*sqrt_I);
                end
                
                term_ca = term_ca + m_c * m_a * (B_phi + sqrt_I * p.cphi);
            end
        end
    end
    
    % Theta terms (cation-cation and anion-anion)
    term_theta = 0;
    
    % Cation-cation
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
    
    % Anion-anion
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
    
    % Psi terms (ternary interactions)
    term_psi = 0;
    
    % Cation-cation-anion
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
    
    % Cation-anion-anion
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
    
    % Calculate total sum of molalities
    sum_m = 0;
    ion_names = keys(comp);
    for i = 1:length(ion_names)
        sum_m = sum_m + comp(ion_names{i});
    end
    
    % Osmotic coefficient
    phi = 1 + (2/sum_m) * (I * f_phi + term_ca + term_theta + term_psi);
end

function charges = get_ion_charges()
    % Define charges for common ions
    charges = containers.Map();
    
    % Cations
    charges('H+') = 1;
    charges('Li+') = 1;
    charges('Na+') = 1;
    charges('K+') = 1;
    charges('Rb+') = 1;
    charges('Cs+') = 1;
    charges('NH4+') = 1;
    
    charges('Mg++') = 2;
    charges('Ca++') = 2;
    charges('Sr++') = 2;
    charges('Ba++') = 2;
    charges('Mn++') = 2;
    charges('Fe++') = 2;
    charges('Co++') = 2;
    charges('Ni++') = 2;
    charges('Cu++') = 2;
    charges('Zn++') = 2;
    
    charges('Al+++') = 3;
    charges('Cr+++') = 3;
    charges('Fe+++') = 3;
    
    % Anions
    charges('F-') = -1;
    charges('Cl-') = -1;
    charges('Br-') = -1;
    charges('I-') = -1;
    charges('OH-') = -1;
    charges('NO3-') = -1;
    charges('ClO4-') = -1;
    charges('HSO4-') = -1;
    charges('H2PO4-') = -1;
    charges('HCO3-') = -1;
    
    charges('SO4--') = -2;
    charges('CO3--') = -2;
    charges('HPO4--') = -2;
    
    charges('PO4---') = -3;
end

function err = calculate_aw_error(m_total, RH_target, system, params, T, MWw)
    % Calculate error between target RH and calculated water activity
    
    if m_total <= 0 || m_total > 100
        err = 1e10;
        return;
    end
    
    try
        composition = scale_composition(system.ions, m_total);
        aw_calc = calculate_water_activity_pitzer(composition, T, params, MWw);
        err = aw_calc - RH_target;
    catch
        err = 1e10;
    end
end
