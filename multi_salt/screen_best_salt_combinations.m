close all
clear
clc

% Script to systematically screen all possible binary salt combinations
% at 70% RH with constant total salt mass to find the best water uptake

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
RH_target = 0.70; % 70% RH
MWw = 18.015; % Molecular weight of water (g/mol)
total_salt_mass = 100; % grams (constant for all comparisons)

%% Define available salts
% Only include salts with Pitzer data and calculate_mf functions
available_salts = {
    struct('name', 'LiCl', 'MW', 42.394, 'ions', containers.Map({'Li+', 'Cl-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_LiCl', 'needs_T', true);
    struct('name', 'LiBr', 'MW', 86.845, 'ions', containers.Map({'Li+', 'Br-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_LiBr', 'needs_T', false);
    struct('name', 'NaCl', 'MW', 58.443, 'ions', containers.Map({'Na+', 'Cl-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_NaCl', 'needs_T', false);
    struct('name', 'NaBr', 'MW', 102.894, 'ions', containers.Map({'Na+', 'Br-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_NaBr', 'needs_T', false);
    struct('name', 'KCl', 'MW', 74.551, 'ions', containers.Map({'K+', 'Cl-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_KCl', 'needs_T', false);
    struct('name', 'KBr', 'MW', 119.002, 'ions', containers.Map({'K+', 'Br-'}, {1, 1}), 'has_calc_func', true, 'func', 'calculate_mf_KBr', 'needs_T', false);
    struct('name', 'MgCl2', 'MW', 95.211, 'ions', containers.Map({'Mg++', 'Cl-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_MgCl', 'needs_T', false);
    struct('name', 'CaCl2', 'MW', 110.984, 'ions', containers.Map({'Ca++', 'Cl-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_CaCl', 'needs_T', true);
    struct('name', 'CaBr2', 'MW', 199.886, 'ions', containers.Map({'Ca++', 'Br-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_CaBr2', 'needs_T', false);
    struct('name', 'SrCl2', 'MW', 158.526, 'ions', containers.Map({'Sr++', 'Cl-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_SrCl2', 'needs_T', false);
    struct('name', 'BaCl2', 'MW', 208.23, 'ions', containers.Map({'Ba++', 'Cl-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_BaCl2', 'needs_T', false);
    struct('name', 'ZnCl2', 'MW', 136.286, 'ions', containers.Map({'Zn++', 'Cl-'}, {1, 2}), 'has_calc_func', true, 'func', 'calculate_mf_ZnCl', 'needs_T', false);
};

fprintf('=== SCREENING SALT COMBINATIONS FOR MAXIMUM WATER UPTAKE ===\n');
fprintf('Target RH: %.0f%%\n', RH_target*100);
fprintf('Total salt mass: %.0f g\n', total_salt_mass);
fprintf('Temperature: %.0f°C\n\n', T);

%% 1. Screen Pure Salts
fprintf('Screening pure salts...\n');
pure_results = [];

for i = 1:length(available_salts)
    salt = available_salts{i};
    
    try
        % Calculate using analytical function
        if salt.needs_T
            mf_salt = feval(salt.func, RH_target, T);
        else
            mf_salt = feval(salt.func, RH_target);
        end
        
        % Calculate water uptake
        mass_water = (1 - mf_salt) / mf_salt * total_salt_mass;
        mf_water = mass_water / (mass_water + total_salt_mass);
        
        % Calculate molality
        molality = (total_salt_mass / salt.MW) / (mass_water / 1000);
        
        % Calculate mole fraction
        ion_names = keys(salt.ions);
        n_ions_total = 0;
        for j = 1:length(ion_names)
            n_ions_total = n_ions_total + salt.ions(ion_names{j});
        end
        n_salt = total_salt_mass / salt.MW;
        n_water = mass_water / MWw;
        x_water = n_water / (n_water + n_ions_total * n_salt);
        
        pure_results = [pure_results; struct(...
            'name', salt.name, ...
            'type', 'pure', ...
            'mass_water', mass_water, ...
            'mf_water', mf_water, ...
            'x_water', x_water, ...
            'molality', molality, ...
            'composition', salt.name)];
        
        fprintf('  %10s: %.2f g water (%.1f%% water, x_w=%.4f)\n', ...
            salt.name, mass_water, mf_water*100, x_water);
    catch ME
        fprintf('  %10s: Failed (%s)\n', salt.name, ME.message);
    end
end

%% 2. Screen Binary Combinations (50:50 mass ratio)
fprintf('\nScreening binary combinations (50:50 mass)...\n');
binary_results = [];
n_attempted = 0;
n_success = 0;

for i = 1:length(available_salts)
    for j = (i+1):length(available_salts)
        salt1 = available_salts{i};
        salt2 = available_salts{j};
        n_attempted = n_attempted + 1;
        
        % Check if we have Pitzer parameters for this combination
        if ~check_pitzer_available(salt1, salt2, params)
            continue;
        end
        
        try
            % Define the binary mixture (50:50 mass)
            system = struct(...
                'name', [salt1.name '_' salt2.name], ...
                'mass_fractions', [0.5, 0.5], ...
                'salt_MWs', [salt1.MW, salt2.MW], ...
                'salt_ions', {{salt1.ions, salt2.ions}});
            
            % Find the molality that gives RH_target
            options = optimset('Display', 'off', 'TolX', 1e-6, 'MaxIter', 100);
            
            m_result = fzero(@(m) calculate_aw_error_mass_basis(m, RH_target, system, params, T, MWw), ...
                             1.0, options);
            
            if m_result > 0 && m_result < 100
                % Calculate water mass
                mass_water = 1000 * (total_salt_mass / (sum(system.mass_fractions .* system.salt_MWs))) / m_result;
                mf_water = mass_water / (mass_water + total_salt_mass);
                
                % Calculate mole fraction
                composition = get_ion_composition_from_mass_basis(system, m_result);
                ion_names = keys(composition);
                n_ions_total = 0;
                for k = 1:length(ion_names)
                    n_ions_total = n_ions_total + composition(ion_names{k});
                end
                n_water = (mass_water / 1000) * 55.51;
                x_water = n_water / (n_water + n_ions_total);
                
                binary_results = [binary_results; struct(...
                    'name', [salt1.name '+' salt2.name], ...
                    'type', 'binary', ...
                    'mass_water', mass_water, ...
                    'mf_water', mf_water, ...
                    'x_water', x_water, ...
                    'molality', m_result, ...
                    'composition', sprintf('%s (50%%) + %s (50%%)', salt1.name, salt2.name))];
                
                n_success = n_success + 1;
            end
        catch ME
            % Silent failure for combinations that don't converge
        end
    end
end

fprintf('  Successfully evaluated %d out of %d binary combinations\n', n_success, n_attempted);

%% 3. Combine and Rank Results
all_results = [pure_results; binary_results];

% Sort by mass of water (descending)
[~, sort_idx] = sort([all_results.mass_water], 'descend');
all_results = all_results(sort_idx);

%% 4. Display Top Results
fprintf('\n=== TOP 20 SALT COMBINATIONS ===\n');
fprintf('%-5s %-30s %-12s %-12s %-12s %-10s\n', 'Rank', 'System', 'Water (g)', 'MF Water', 'x_water', 'Molality');
fprintf('%s\n', repmat('-', 1, 90));

n_display = min(20, length(all_results));
for i = 1:n_display
    r = all_results(i);
    fprintf('%-5d %-30s %10.2f g %10.1f%% %11.4f %10.2f\n', ...
        i, r.name, r.mass_water, r.mf_water*100, r.x_water, r.molality);
end

%% 5. Statistics and Analysis
fprintf('\n=== STATISTICS ===\n');

% Best pure salt
pure_idx = strcmp({all_results.type}, 'pure');
best_pure = all_results(find(pure_idx, 1, 'first'));
fprintf('Best Pure Salt: %s (%.2f g water, %.1f%%)\n', ...
    best_pure.name, best_pure.mass_water, best_pure.mf_water*100);

% Best binary
binary_idx = strcmp({all_results.type}, 'binary');
if any(binary_idx)
    best_binary = all_results(find(binary_idx, 1, 'first'));
    fprintf('Best Binary Mix: %s (%.2f g water, %.1f%%)\n', ...
        best_binary.name, best_binary.mass_water, best_binary.mf_water*100);
    
    improvement = (best_binary.mass_water - best_pure.mass_water) / best_pure.mass_water * 100;
    fprintf('Improvement: %.1f%% more water than best pure salt\n', improvement);
end

%% 6. Create Visualizations

% Figure 1: Top 15 systems comparison
n_plot = min(15, length(all_results));
top_systems = all_results(1:n_plot);

figure('Position', [100, 100, 1400, 800]);
subplot(2,1,1);
hold on; grid on; box on;

x_pos = 1:n_plot;
colors_plot = zeros(n_plot, 3);
for i = 1:n_plot
    if strcmp(top_systems(i).type, 'pure')
        colors_plot(i,:) = [0.3, 0.5, 0.8]; % Blue for pure
    else
        colors_plot(i,:) = [0.8, 0.4, 0.2]; % Orange for binary
    end
end

b = bar(x_pos, [top_systems.mass_water], 'FaceColor', 'flat');
b.CData = colors_plot;

ylabel('Water Uptake (g)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Top %d Salt Systems at %.0f%% RH (%.0f g salt)', n_plot, RH_target*100, total_salt_mass), ...
      'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'XTick', x_pos, 'XTickLabel', {top_systems.name}, 'XTickLabelRotation', 45);
set(gca, 'FontSize', 12);
legend({'Pure Salt', 'Binary Mix'}, 'Location', 'northeast');

subplot(2,1,2);
hold on; grid on; box on;

b2 = bar(x_pos, [top_systems.mf_water]*100, 'FaceColor', 'flat');
b2.CData = colors_plot;

ylabel('Water Mass Fraction (%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('System', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XTick', x_pos, 'XTickLabel', {top_systems.name}, 'XTickLabelRotation', 45);
set(gca, 'FontSize', 12);

set(gcf, 'color', 'w');
saveas(gcf, fullfile(fig_out_dir, 'top_salt_combinations_screening.png'));
savefig(fullfile(fig_out_dir, 'top_salt_combinations_screening.fig'));

% Figure 2: Pure vs Binary comparison
figure('Position', [100, 100, 1200, 600]);

subplot(1,2,1);
pure_waters = [all_results(pure_idx).mass_water];
binary_waters = [all_results(binary_idx).mass_water];

boxplot([pure_waters'; binary_waters'], [ones(length(pure_waters),1); 2*ones(length(binary_waters),1)], ...
        'Labels', {'Pure Salts', 'Binary Mixtures'});
ylabel('Water Uptake (g)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Water Uptake Distribution at %.0f%% RH', RH_target*100), 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);

subplot(1,2,2);
histogram(pure_waters, 10, 'FaceColor', [0.3, 0.5, 0.8], 'FaceAlpha', 0.7, 'DisplayName', 'Pure');
hold on;
histogram(binary_waters, 10, 'FaceColor', [0.8, 0.4, 0.2], 'FaceAlpha', 0.7, 'DisplayName', 'Binary');
xlabel('Water Uptake (g)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Count', 'FontSize', 14, 'FontWeight', 'bold');
title('Distribution Comparison', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);

set(gcf, 'color', 'w');
saveas(gcf, fullfile(fig_out_dir, 'pure_vs_binary_screening.png'));
savefig(fullfile(fig_out_dir, 'pure_vs_binary_screening.fig'));

%% 7. Save Results to CSV
results_table = struct2table(all_results);
writetable(results_table, fullfile(fig_out_dir, 'salt_screening_results.csv'));
fprintf('\nResults saved to: %s\n', fullfile(fig_out_dir, 'salt_screening_results.csv'));

fprintf('\nScreening complete!\n');

%% Helper Functions

function available = check_pitzer_available(salt1, salt2, params)
    % Check if Pitzer parameters are available for this salt combination
    % Need to check all ion-ion interactions
    
    ion_names1 = keys(salt1.ions);
    ion_names2 = keys(salt2.ions);
    
    charges = get_ion_charges();
    
    % Check all cation-anion pairs
    for i = 1:length(ion_names1)
        ion1 = ion_names1{i};
        z1 = charges(ion1);
        
        for j = 1:length(ion_names2)
            ion2 = ion_names2{j};
            z2 = charges(ion2);
            
            % Only check cation-anion pairs
            if sign(z1) ~= sign(z2)
                key = [ion1 '_' ion2];
                if ~isKey(params.binary, key)
                    available = false;
                    return;
                end
            end
        end
    end
    
    available = true;
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

function err = calculate_aw_error_mass_basis(m_total, RH_target, system, params, T, MWw)
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

function params_map = create_binary_map(table_data)
    params_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = char(table_data.species1(i));
        sp2 = char(table_data.species2(i));
        
        if isempty(sp1) || isempty(sp2) || strcmp(sp1, 'x')
            continue;
        end
        
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
        
        param_struct = struct('beta0', beta0, 'beta1', beta1, ...
                             'beta2', beta2, 'cphi', cphi);
        
        params_map(key1) = param_struct;
        params_map(key2) = param_struct;
    end
end

function theta_map = create_theta_map(table_data)
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
    psi_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    for i = 1:height(table_data)
        sp1 = char(table_data.species1(i));
        sp2 = char(table_data.species2(i));
        sp3 = char(table_data.species3(i));
        
        if isempty(sp1) || isempty(sp2) || isempty(sp3)
            continue;
        end
        
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

function aw = calculate_water_activity_pitzer(composition, T, params, MWw)
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
                    alpha1 = 1.4;
                    alpha2 = 12.0;
                elseif abs(z_c) >= 3 || abs(z_a) >= 3
                    alpha1 = 2.0;
                    alpha2 = 12.0;
                else
                    alpha1 = 2.0;
                    alpha2 = 12.0;
                end
                
                B_phi = p.beta0 + p.beta1 * exp(-alpha1*sqrt_I);
                
                if abs(p.beta2) > 1e-10
                    B_phi = B_phi + p.beta2 * exp(-alpha2*sqrt_I);
                end
                
                term_ca = term_ca + m_c * m_a * (B_phi + sqrt_I * p.cphi);
            end
        end
    end
    
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

function charges = get_ion_charges()
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
