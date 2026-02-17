close all
clear
clc

% Script: calculate_sb_mixture_swings.m
% Calculate water uptake swings for all viable Santa Barbara salt mixtures
% Santa Barbara RH range: 42-88%

[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

fprintf('=== SANTA BARBARA SALT MIXTURE SWING CALCULATION ===\n\n');

%% 1. Define Santa Barbara Conditions
sb_min_RH = 0.42;
sb_max_RH = 0.88;
sb_avg_T = 16;  % °C

% Create RH sweep across Santa Barbara range
n_points = 50;
RH_vec = linspace(sb_min_RH, sb_max_RH, n_points);

fprintf('Santa Barbara conditions:\n');
fprintf('  RH range: %.1f%% - %.1f%%\n', sb_min_RH*100, sb_max_RH*100);
fprintf('  Temperature: %.1f°C (assumed constant)\n', sb_avg_T);
fprintf('  Evaluating %d RH points\n\n', n_points);

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

%% 3. Define Viable Mixtures
MW_water = 18.015;

% Molecular weights
MW_map = containers.Map();
MW_map('HCl') = 36.461;
MW_map('LiCl') = 42.394;
MW_map('NaCl') = 58.443;
MW_map('KCl') = 74.551;
MW_map('CsCl') = 168.36;
MW_map('CaCl2') = 110.98;
MW_map('MgCl2') = 95.2;
MW_map('NaBr') = 102.89;
MW_map('KBr') = 119.00;
MW_map('MgSO4') = 120.37;
MW_map('NaOH') = 40.00;

% Define all 12 viable mixtures from screening results
mixtures = {};

% 1. CaCl2 + KCl
mixtures{end+1} = struct('name', 'CaCl2+KCl', 'salts', {{'CaCl2', 'KCl'}}, ...
    'ions', {{containers.Map({'Ca++', 'Cl-'}, {1, 2}), containers.Map({'K+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.31, 'drh2', 0.855);

% 2. CaCl2 + NaCl
mixtures{end+1} = struct('name', 'CaCl2+NaCl', 'salts', {{'CaCl2', 'NaCl'}}, ...
    'ions', {{containers.Map({'Ca++', 'Cl-'}, {1, 2}), containers.Map({'Na+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.31, 'drh2', 0.765);

% 3. MgCl2 + MgSO4 (anion-anion mixture)
mixtures{end+1} = struct('name', 'MgCl2+MgSO4', 'salts', {{'MgCl2', 'MgSO4'}}, ...
    'ions', {{containers.Map({'Mg++', 'Cl-'}, {1, 2}), containers.Map({'Mg++', 'SO4--'}, {1, 1})}}, ...
    'drh1', 0.33, 'drh2', 0.86);

% 4. NaCl + NaOH (anion-anion mixture)
mixtures{end+1} = struct('name', 'NaCl+NaOH', 'salts', {{'NaCl', 'NaOH'}}, ...
    'ions', {{containers.Map({'Na+', 'Cl-'}, {1, 1}), containers.Map({'Na+', 'OH-'}, {1, 1})}}, ...
    'drh1', 0.765, 'drh2', 0.23);

% 5. HCl + CsCl
mixtures{end+1} = struct('name', 'HCl+CsCl', 'salts', {{'HCl', 'CsCl'}}, ...
    'ions', {{containers.Map({'H+', 'Cl-'}, {1, 1}), containers.Map({'Cs+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.17, 'drh2', 0.82);

% 6. HCl + KCl
mixtures{end+1} = struct('name', 'HCl+KCl', 'salts', {{'HCl', 'KCl'}}, ...
    'ions', {{containers.Map({'H+', 'Cl-'}, {1, 1}), containers.Map({'K+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.17, 'drh2', 0.855);

% 7. HCl + LiCl
mixtures{end+1} = struct('name', 'HCl+LiCl', 'salts', {{'HCl', 'LiCl'}}, ...
    'ions', {{containers.Map({'H+', 'Cl-'}, {1, 1}), containers.Map({'Li+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.17, 'drh2', 0.12);

% 8. HCl + NaCl
mixtures{end+1} = struct('name', 'HCl+NaCl', 'salts', {{'HCl', 'NaCl'}}, ...
    'ions', {{containers.Map({'H+', 'Cl-'}, {1, 1}), containers.Map({'Na+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.17, 'drh2', 0.765);

% 9. LiCl + CsCl
mixtures{end+1} = struct('name', 'LiCl+CsCl', 'salts', {{'LiCl', 'CsCl'}}, ...
    'ions', {{containers.Map({'Li+', 'Cl-'}, {1, 1}), containers.Map({'Cs+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.12, 'drh2', 0.82);

% 10. LiCl + KCl
mixtures{end+1} = struct('name', 'LiCl+KCl', 'salts', {{'LiCl', 'KCl'}}, ...
    'ions', {{containers.Map({'Li+', 'Cl-'}, {1, 1}), containers.Map({'K+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.12, 'drh2', 0.855);

% 11. LiCl + NaCl
mixtures{end+1} = struct('name', 'LiCl+NaCl', 'salts', {{'LiCl', 'NaCl'}}, ...
    'ions', {{containers.Map({'Li+', 'Cl-'}, {1, 1}), containers.Map({'Na+', 'Cl-'}, {1, 1})}}, ...
    'drh1', 0.12, 'drh2', 0.765);

% 12. NaBr + KBr
mixtures{end+1} = struct('name', 'NaBr+KBr', 'salts', {{'NaBr', 'KBr'}}, ...
    'ions', {{containers.Map({'Na+', 'Br-'}, {1, 1}), containers.Map({'K+', 'Br-'}, {1, 1})}}, ...
    'drh1', 0.614, 'drh2', 0.833);

%% 4. Calculate Swings for Each Mixture
fprintf('%s\n', repmat('=', 1, 80));
fprintf('CALCULATING WATER UPTAKE SWINGS (50:50 mass ratio)\n');
fprintf('%s\n', repmat('=', 1, 80));

results = struct('name', {}, 'swing', {}, 'min_uptake', {}, 'max_uptake', {}, ...
    'min_RH', {}, 'max_RH', {}, 'success_rate', {}, 'drh_eff', {}, ...
    'swing_molar', {}, 'min_uptake_molar', {}, 'max_uptake_molar', {});

for m = 1:length(mixtures)
    mix = mixtures{m};
    fprintf('\n[%d/%d] Processing %s...\n', m, length(mixtures), mix.name);
    
    % Set up system
    system = struct(...
        'name', mix.name, ...
        'mass_fractions', [0.5, 0.5], ...
        'salt_MWs', [MW_map(mix.salts{1}), MW_map(mix.salts{2})], ...
        'salt_ions', {mix.ions});
    
    % Determine effective DRH (max of the two salts)
    drh_eff = max(mix.drh1, mix.drh2);
    
    % Calculate water uptake across RH range
    water_uptake = nan(size(RH_vec));
    water_uptake_molar = nan(size(RH_vec));  % moles H2O per mole salt
    
    options = optimset('Display', 'off', 'TolX', 1e-6, 'MaxIter', 100);
    n_success = 0;
    
    for i = 1:length(RH_vec)
        RH = RH_vec(i);
        
        % Skip if below effective DRH
        if RH < drh_eff
            continue;
        end
        
        try
            % Solve for molality
            m_result = fzero(@(m) calculate_aw_error(m, RH, system, params, sb_avg_T, MW_water), ...
                             2.0, options);
            
            if m_result > 0 && m_result < 100
                % Calculate water uptake (mass basis)
                total_salt_mass = 100;
                mass_water = 1000 * (total_salt_mass / sum(system.mass_fractions .* system.salt_MWs)) / m_result;
                water_uptake(i) = mass_water / (mass_water + total_salt_mass);
                
                % Calculate water uptake (molar basis: moles H2O per mole salt)
                moles_water = mass_water / MW_water;
                moles_salt = sum(system.mass_fractions * total_salt_mass ./ system.salt_MWs);
                water_uptake_molar(i) = moles_water / moles_salt;
                
                n_success = n_success + 1;
            end
        catch
            % Failed to converge
        end
    end
    
    % Calculate swing if we have enough data
    valid_uptake = water_uptake(~isnan(water_uptake));
    valid_uptake_molar = water_uptake_molar(~isnan(water_uptake_molar));
    valid_RH = RH_vec(~isnan(water_uptake));
    
    if length(valid_uptake) >= 5
        min_uptake = min(valid_uptake);
        max_uptake = max(valid_uptake);
        swing = max_uptake - min_uptake;
        
        min_uptake_molar = min(valid_uptake_molar);
        max_uptake_molar = max(valid_uptake_molar);
        swing_molar = max_uptake_molar - min_uptake_molar;
        
        [~, idx_min] = min(water_uptake);
        [~, idx_max] = max(water_uptake);
        
        success_rate = n_success / length(RH_vec) * 100;
        
        % Store results
        results(end+1).name = mix.name;
        results(end).swing = swing;
        results(end).min_uptake = min_uptake;
        results(end).max_uptake = max_uptake;
        results(end).min_RH = RH_vec(idx_min);
        results(end).max_RH = RH_vec(idx_max);
        results(end).success_rate = success_rate;
        results(end).drh_eff = drh_eff;
        results(end).swing_molar = swing_molar;
        results(end).min_uptake_molar = min_uptake_molar;
        results(end).max_uptake_molar = max_uptake_molar;
        
        fprintf('  Swing (mass): %.4f (%.2f%% range)\n', swing, swing*100);
        fprintf('  Swing (molar): %.2f mol H2O/mol salt\n', swing_molar);
        fprintf('  Min: %.4f at RH=%.1f%%, Max: %.4f at RH=%.1f%%\n', ...
            min_uptake, RH_vec(idx_min)*100, max_uptake, RH_vec(idx_max)*100);
        fprintf('  Success rate: %.1f%%\n', success_rate);
    else
        fprintf('  FAILED: Insufficient data points (%d valid)\n', length(valid_uptake));
    end
end

%% 5. Calculate Pure Salt Swings for Comparison
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('CALCULATING PURE SALT SWINGS FOR COMPARISON\n');
fprintf('%s\n', repmat('=', 1, 80));

% Get unique salts from all mixtures
unique_salts = {};
salt_drh_map = containers.Map();
for m = 1:length(mixtures)
    for s = 1:length(mixtures{m}.salts)
        salt = mixtures{m}.salts{s};
        if ~any(strcmp(unique_salts, salt))
            unique_salts{end+1} = salt;
            if s == 1
                salt_drh_map(salt) = mixtures{m}.drh1;
            else
                salt_drh_map(salt) = mixtures{m}.drh2;
            end
        end
    end
end

fprintf('Evaluating %d unique salts...\n', length(unique_salts));

pure_results = struct('name', {}, 'swing', {}, 'min_uptake', {}, 'max_uptake', {}, ...
    'min_RH', {}, 'max_RH', {}, 'success_rate', {}, 'drh', {}, ...
    'swing_molar', {}, 'min_uptake_molar', {}, 'max_uptake_molar', {});

% Define ion stoichiometry for each salt
salt_ion_map = containers.Map();
salt_ion_map('HCl') = containers.Map({'H+', 'Cl-'}, {1, 1});
salt_ion_map('LiCl') = containers.Map({'Li+', 'Cl-'}, {1, 1});
salt_ion_map('NaCl') = containers.Map({'Na+', 'Cl-'}, {1, 1});
salt_ion_map('KCl') = containers.Map({'K+', 'Cl-'}, {1, 1});
salt_ion_map('CsCl') = containers.Map({'Cs+', 'Cl-'}, {1, 1});
salt_ion_map('CaCl2') = containers.Map({'Ca++', 'Cl-'}, {1, 2});
salt_ion_map('MgCl2') = containers.Map({'Mg++', 'Cl-'}, {1, 2});
salt_ion_map('NaBr') = containers.Map({'Na+', 'Br-'}, {1, 1});
salt_ion_map('KBr') = containers.Map({'K+', 'Br-'}, {1, 1});
salt_ion_map('MgSO4') = containers.Map({'Mg++', 'SO4--'}, {1, 1});
salt_ion_map('NaOH') = containers.Map({'Na+', 'OH-'}, {1, 1});

for s = 1:length(unique_salts)
    salt_name = unique_salts{s};
    fprintf('\n[%d/%d] Processing pure %s...\n', s, length(unique_salts), salt_name);
    
    % Set up pure salt system
    system = struct(...
        'name', salt_name, ...
        'mass_fractions', 1.0, ...
        'salt_MWs', MW_map(salt_name), ...
        'salt_ions', {{salt_ion_map(salt_name)}});
    
    drh = salt_drh_map(salt_name);
    
    % Calculate water uptake across RH range
    water_uptake = nan(size(RH_vec));
    water_uptake_molar = nan(size(RH_vec));  % moles H2O per mole salt
    
    options = optimset('Display', 'off', 'TolX', 1e-6, 'MaxIter', 100);
    n_success = 0;
    
    for i = 1:length(RH_vec)
        RH = RH_vec(i);
        
        % Skip if below DRH
        if RH < drh
            continue;
        end
        
        try
            % Solve for molality
            m_result = fzero(@(m) calculate_aw_error(m, RH, system, params, sb_avg_T, MW_water), ...
                             2.0, options);
            
            if m_result > 0 && m_result < 100
                % Calculate water uptake (mass basis)
                total_salt_mass = 100;
                mass_water = 1000 * (total_salt_mass / system.salt_MWs) / m_result;
                water_uptake(i) = mass_water / (mass_water + total_salt_mass);
                
                % Calculate water uptake (molar basis: moles H2O per mole salt)
                moles_water = mass_water / MW_water;
                moles_salt = total_salt_mass / system.salt_MWs;
                water_uptake_molar(i) = moles_water / moles_salt;
                
                n_success = n_success + 1;
            end
        catch
            % Failed to converge
        end
    end
    
    % Calculate swing if we have enough data
    valid_uptake = water_uptake(~isnan(water_uptake));
    valid_uptake_molar = water_uptake_molar(~isnan(water_uptake_molar));
    
    if length(valid_uptake) >= 5
        min_uptake = min(valid_uptake);
        max_uptake = max(valid_uptake);
        swing = max_uptake - min_uptake;
        
        min_uptake_molar = min(valid_uptake_molar);
        max_uptake_molar = max(valid_uptake_molar);
        swing_molar = max_uptake_molar - min_uptake_molar;
        
        [~, idx_min] = min(water_uptake);
        [~, idx_max] = max(water_uptake);
        
        success_rate = n_success / length(RH_vec) * 100;
        
        % Store results
        pure_results(end+1).name = salt_name;
        pure_results(end).swing = swing;
        pure_results(end).min_uptake = min_uptake;
        pure_results(end).max_uptake = max_uptake;
        pure_results(end).min_RH = RH_vec(idx_min);
        pure_results(end).max_RH = RH_vec(idx_max);
        pure_results(end).success_rate = success_rate;
        pure_results(end).drh = drh;
        pure_results(end).swing_molar = swing_molar;
        pure_results(end).min_uptake_molar = min_uptake_molar;
        pure_results(end).max_uptake_molar = max_uptake_molar;
        
        fprintf('  Swing (mass): %.4f (%.2f%% range)\n', swing, swing*100);
        fprintf('  Swing (molar): %.2f mol H2O/mol salt\n', swing_molar);
        fprintf('  Min: %.4f at RH=%.1f%%, Max: %.4f at RH=%.1f%%\n', ...
            min_uptake, RH_vec(idx_min)*100, max_uptake, RH_vec(idx_max)*100);
        fprintf('  Success rate: %.1f%%\n', success_rate);
    else
        fprintf('  FAILED: Insufficient data points (%d valid)\n', length(valid_uptake));
    end
end

%% 6. Sort and Display Results with Comparison
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('RESULTS SUMMARY - MIXTURES RANKED BY SWING\n');
fprintf('%s\n', repmat('=', 1, 80));

if isempty(results)
    fprintf('No mixtures successfully calculated.\n');
else
    % Sort by swing (descending)
    [~, sort_idx] = sort([results.swing], 'descend');
    results_sorted = results(sort_idx);
    
    fprintf('\n%-15s  Swing(%%w)  Min(%%w)  Max(%%w)  MinRH  MaxRH  DRH_eff  vs Pure\n', 'Mixture');
    fprintf('%s\n', repmat('-', 1, 88));
    
    for i = 1:length(results_sorted)
        r = results_sorted(i);
        
        % Find corresponding pure salts and calculate average pure swing
        parts = strsplit(r.name, '+');
        pure_swings = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_swings(end+1) = pure_results(idx).swing;
            end
        end
        
        if ~isempty(pure_swings)
            avg_pure_swing = mean(pure_swings);
            improvement = ((r.swing - avg_pure_swing) / avg_pure_swing) * 100;
            vs_pure_str = sprintf('%+.1f%%', improvement);
        else
            vs_pure_str = 'N/A';
        end
        
        fprintf('%-15s  %7.2f    %6.2f   %6.2f   %5.1f%% %5.1f%%  %5.1f%%  %s\n', ...
            r.name, r.swing*100, r.min_uptake*100, r.max_uptake*100, ...
            r.min_RH*100, r.max_RH*100, r.drh_eff*100, vs_pure_str);
    end
    
    fprintf('\n');
    fprintf('Top 3 mixtures by swing:\n');
    for i = 1:min(3, length(results_sorted))
        fprintf('  %d. %s: %.2f%% swing\n', i, results_sorted(i).name, results_sorted(i).swing*100);
    end
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('PURE SALT RESULTS - RANKED BY SWING\n');
fprintf('%s\n', repmat('=', 1, 80));

if ~isempty(pure_results)
    [~, sort_idx] = sort([pure_results.swing], 'descend');
    pure_sorted = pure_results(sort_idx);
    
    fprintf('\n%-10s  Swing(%%w)  Min(%%w)  Max(%%w)  MinRH  MaxRH  DRH\n', 'Salt');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for i = 1:length(pure_sorted)
        r = pure_sorted(i);
        fprintf('%-10s  %7.2f    %6.2f   %6.2f   %5.1f%% %5.1f%%  %5.1f%%\n', ...
            r.name, r.swing*100, r.min_uptake*100, r.max_uptake*100, ...
            r.min_RH*100, r.max_RH*100, r.drh*100);
    end
    
    fprintf('\n');
    fprintf('Top 3 pure salts by swing:\n');
    for i = 1:min(3, length(pure_sorted))
        fprintf('  %d. %s: %.2f%% swing\n', i, pure_sorted(i).name, pure_sorted(i).swing*100);
    end
end

%% 7. Create Visualization
if ~isempty(results)
    fprintf('\nCreating visualization...\n');
    
    fig = figure('Position', [100, 100, 1600, 1000]);
    
    % Plot 1: Mixture vs Pure Salt Swings
    subplot(2, 3, 1);
    hold on;
    
    % Plot mixture swings
    x_mix = 1:length(results_sorted);
    bar(x_mix, [results_sorted.swing] * 100, 'FaceColor', [0.2 0.6 0.8]);
    
    % Plot pure salt swings for comparison (as lines)
    for i = 1:length(results_sorted)
        parts = strsplit(results_sorted(i).name, '+');
        pure_vals = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_vals(end+1) = pure_results(idx).swing * 100;
            end
        end
        if ~isempty(pure_vals)
            plot([i-0.3, i+0.3], [pure_vals(1), pure_vals(1)], 'r-', 'LineWidth', 2);
            if length(pure_vals) > 1
                plot([i-0.3, i+0.3], [pure_vals(2), pure_vals(2)], 'g-', 'LineWidth', 2);
            end
        end
    end
    
    set(gca, 'XTick', 1:length(results_sorted), 'XTickLabel', {results_sorted.name}, 'XTickLabelRotation', 45);
    ylabel('Water Uptake Swing (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Mixture vs Pure Salt Swings', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Mixture', 'Pure Salt 1', 'Pure Salt 2', 'Location', 'best');
    grid on;
    hold off;
    
    % Plot 2: Pure Salt Swings
    subplot(2, 3, 2);
    if ~isempty(pure_results)
        bar([pure_sorted.swing] * 100, 'FaceColor', [0.8 0.4 0.2]);
        set(gca, 'XTickLabel', {pure_sorted.name}, 'XTickLabelRotation', 45);
        ylabel('Water Uptake Swing (%)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Pure Salt Performance', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
    end
    
    % Plot 3: Improvement over Pure Salts
    subplot(2, 3, 3);
    improvements = [];
    labels = {};
    for i = 1:length(results_sorted)
        parts = strsplit(results_sorted(i).name, '+');
        pure_swings = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_swings(end+1) = pure_results(idx).swing;
            end
        end
        if ~isempty(pure_swings)
            avg_pure = mean(pure_swings);
            improvement = ((results_sorted(i).swing - avg_pure) / avg_pure) * 100;
            improvements(end+1) = improvement;
            labels{end+1} = results_sorted(i).name;
        end
    end
    
    bar_colors = zeros(length(improvements), 3);
    for i = 1:length(improvements)
        if improvements(i) > 0
            bar_colors(i, :) = [0.2 0.8 0.2];  % Green for positive
        else
            bar_colors(i, :) = [0.8 0.2 0.2];  % Red for negative
        end
    end
    
    b = bar(improvements);
    b.FaceColor = 'flat';
    b.CData = bar_colors;
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Improvement over Avg Pure (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Mixture Improvement vs Pure Salts', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    hold on;
    plot([0 length(improvements)+1], [0 0], 'k--', 'LineWidth', 1);
    hold off;
    
    % Plot 4: Min vs Max uptake
    subplot(2, 3, 4);
    hold on;
    plot([results_sorted.min_uptake]*100, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Min');
    plot([results_sorted.max_uptake]*100, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Max');
    set(gca, 'XTickLabel', {results_sorted.name}, 'XTickLabelRotation', 45);
    ylabel('Water Uptake (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Min and Max Water Uptake', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Plot 5: Swing vs Effective DRH (Mixtures and Pure)
    subplot(2, 3, 5);
    hold on;
    scatter([results_sorted.drh_eff]*100, [results_sorted.swing]*100, 100, 'b', 'filled', 'DisplayName', 'Mixtures');
    if ~isempty(pure_results)
        scatter([pure_sorted.drh]*100, [pure_sorted.swing]*100, 100, 'r', 'filled', 'DisplayName', 'Pure Salts');
    end
    xlabel('DRH (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Swing (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Swing vs DRH', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Plot 6: RH range at min/max
    subplot(2, 3, 6);
    hold on;
    for i = 1:length(results_sorted)
        plot([results_sorted(i).min_RH results_sorted(i).max_RH]*100, [i i], 'b-', 'LineWidth', 3);
        plot(results_sorted(i).min_RH*100, i, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(results_sorted(i).max_RH*100, i, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    xlim([40 90]);
    ylim([0 length(results_sorted)+1]);
    set(gca, 'YTick', 1:length(results_sorted), 'YTickLabel', {results_sorted.name});
    xlabel('Relative Humidity (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('RH at Min (blue) and Max (red) Uptake', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    hold off;
    
    % Save figure
    fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    saveas(fig, fullfile(fig_out_dir, 'santa_barbara_mixture_swings.png'));
    fprintf('Figure saved to: %s\n', fullfile(fig_out_dir, 'santa_barbara_mixture_swings.png'));
    
    %% Create Molar Basis Figure
    fprintf('\nCreating molar basis visualization...\n');
    
    fig_molar = figure('Position', [150, 150, 1600, 1000]);
    
    % Plot 1: Mixture vs Pure Salt Swings (Molar Basis)
    subplot(2, 3, 1);
    hold on;
    
    % Plot mixture swings
    x_mix = 1:length(results_sorted);
    bar(x_mix, [results_sorted.swing_molar], 'FaceColor', [0.2 0.6 0.8]);
    
    % Plot pure salt swings for comparison (as lines)
    for i = 1:length(results_sorted)
        parts = strsplit(results_sorted(i).name, '+');
        pure_vals = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_vals(end+1) = pure_results(idx).swing_molar;
            end
        end
        if ~isempty(pure_vals)
            plot([i-0.3, i+0.3], [pure_vals(1), pure_vals(1)], 'r-', 'LineWidth', 2);
            if length(pure_vals) > 1
                plot([i-0.3, i+0.3], [pure_vals(2), pure_vals(2)], 'g-', 'LineWidth', 2);
            end
        end
    end
    
    set(gca, 'XTick', 1:length(results_sorted), 'XTickLabel', {results_sorted.name}, 'XTickLabelRotation', 45);
    ylabel('Water Uptake Swing (mol H_2O/mol salt)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Mixture vs Pure Salt Swings (Molar Basis)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Mixture', 'Pure Salt 1', 'Pure Salt 2', 'Location', 'best');
    grid on;
    hold off;
    
    % Plot 2: Pure Salt Swings (Molar Basis)
    subplot(2, 3, 2);
    if ~isempty(pure_results)
        bar([pure_sorted.swing_molar], 'FaceColor', [0.8 0.4 0.2]);
        set(gca, 'XTick', 1:length(pure_sorted), 'XTickLabel', {pure_sorted.name}, 'XTickLabelRotation', 45);
        ylabel('Water Uptake Swing (mol H_2O/mol salt)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Pure Salt Performance (Molar Basis)', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
    end
    
    % Plot 3: Improvement over Pure Salts (Molar Basis)
    subplot(2, 3, 3);
    improvements_molar = [];
    labels = {};
    for i = 1:length(results_sorted)
        parts = strsplit(results_sorted(i).name, '+');
        pure_swings = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_swings(end+1) = pure_results(idx).swing_molar;
            end
        end
        if ~isempty(pure_swings)
            avg_pure = mean(pure_swings);
            improvement = ((results_sorted(i).swing_molar - avg_pure) / avg_pure) * 100;
            improvements_molar(end+1) = improvement;
            labels{end+1} = results_sorted(i).name;
        end
    end
    
    bar_colors = zeros(length(improvements_molar), 3);
    for i = 1:length(improvements_molar)
        if improvements_molar(i) > 0
            bar_colors(i, :) = [0.2 0.8 0.2];  % Green for positive
        else
            bar_colors(i, :) = [0.8 0.2 0.2];  % Red for negative
        end
    end
    
    b = bar(improvements_molar);
    b.FaceColor = 'flat';
    b.CData = bar_colors;
    set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45);
    ylabel('Improvement over Avg Pure (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Mixture Improvement vs Pure Salts (Molar)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    hold on;
    plot([0 length(improvements_molar)+1], [0 0], 'k--', 'LineWidth', 1);
    hold off;
    
    % Plot 4: Min vs Max uptake (Molar Basis)
    subplot(2, 3, 4);
    hold on;
    plot([results_sorted.min_uptake_molar], 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Min');
    plot([results_sorted.max_uptake_molar], 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Max');
    set(gca, 'XTick', 1:length(results_sorted), 'XTickLabel', {results_sorted.name}, 'XTickLabelRotation', 45);
    ylabel('Water Uptake (mol H_2O/mol salt)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Min and Max Water Uptake (Molar Basis)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Plot 5: Molar Swing vs Effective DRH (Mixtures and Pure)
    subplot(2, 3, 5);
    hold on;
    scatter([results_sorted.drh_eff]*100, [results_sorted.swing_molar], 100, 'b', 'filled', 'DisplayName', 'Mixtures');
    if ~isempty(pure_results)
        scatter([pure_sorted.drh]*100, [pure_sorted.swing_molar], 100, 'r', 'filled', 'DisplayName', 'Pure Salts');
    end
    xlabel('DRH (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Swing (mol H_2O/mol salt)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Molar Swing vs DRH', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Plot 6: Comparison of Mass vs Molar Rankings
    subplot(2, 3, 6);
    hold on;
    x_pos = 1:length(results_sorted);
    
    % Sort by molar swing for comparison
    [~, molar_sort_idx] = sort([results.swing_molar], 'descend');
    results_molar_sorted = results(molar_sort_idx);
    
    % Find ranking positions
    mass_ranks = zeros(1, length(results_sorted));
    molar_ranks = zeros(1, length(results_sorted));
    for i = 1:length(results_sorted)
        % Find where this mixture appears in molar ranking
        molar_idx = find(strcmp({results_molar_sorted.name}, results_sorted(i).name));
        molar_ranks(i) = molar_idx;
        mass_ranks(i) = i;
    end
    
    plot(x_pos, mass_ranks, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Mass Basis Rank');
    plot(x_pos, molar_ranks, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Molar Basis Rank');
    set(gca, 'XTick', 1:length(results_sorted), 'XTickLabel', {results_sorted.name}, 'XTickLabelRotation', 45);
    ylabel('Rank (1 = best)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Mass vs Molar Basis Rankings', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;
    set(gca, 'YDir', 'reverse');  % Lower rank number at top
    hold off;
    
    % Save molar figure
    saveas(fig_molar, fullfile(fig_out_dir, 'santa_barbara_mixture_swings_molar.png'));
    fprintf('Molar basis figure saved to: %s\n', fullfile(fig_out_dir, 'santa_barbara_mixture_swings_molar.png'));
end

%% 8. Save Results to CSV
if ~isempty(results)
    % Save mixture results
    output_file = fullfile(filepath, 'sb_mixture_swing_results.csv');
    
    name_list = {results_sorted.name}';
    swing_list = [results_sorted.swing]' * 100;
    min_uptake_list = [results_sorted.min_uptake]' * 100;
    max_uptake_list = [results_sorted.max_uptake]' * 100;
    min_RH_list = [results_sorted.min_RH]' * 100;
    max_RH_list = [results_sorted.max_RH]' * 100;
    drh_eff_list = [results_sorted.drh_eff]' * 100;
    
    % Calculate improvement over pure salts
    improvement_list = zeros(length(results_sorted), 1);
    for i = 1:length(results_sorted)
        parts = strsplit(results_sorted(i).name, '+');
        pure_swings = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_swings(end+1) = pure_results(idx).swing;
            end
        end
        if ~isempty(pure_swings)
            avg_pure = mean(pure_swings);
            improvement_list(i) = ((results_sorted(i).swing - avg_pure) / avg_pure) * 100;
        else
            improvement_list(i) = NaN;
        end
    end
    
    result_table = table(name_list, swing_list, min_uptake_list, max_uptake_list, ...
        min_RH_list, max_RH_list, drh_eff_list, improvement_list, ...
        'VariableNames', {'Mixture', 'Swing_pct', 'MinUptake_pct', 'MaxUptake_pct', ...
        'MinRH_pct', 'MaxRH_pct', 'EffectiveDRH_pct', 'ImprovementOverPure_pct'});
    
    writetable(result_table, output_file);
    fprintf('Mixture results saved to: %s\n', output_file);
end

if ~isempty(pure_results)
    % Save pure salt results
    output_file_pure = fullfile(filepath, 'sb_pure_salt_swing_results.csv');
    
    name_list = {pure_sorted.name}';
    swing_list = [pure_sorted.swing]' * 100;
    min_uptake_list = [pure_sorted.min_uptake]' * 100;
    max_uptake_list = [pure_sorted.max_uptake]' * 100;
    min_RH_list = [pure_sorted.min_RH]' * 100;
    max_RH_list = [pure_sorted.max_RH]' * 100;
    drh_list = [pure_sorted.drh]' * 100;
    
    pure_table = table(name_list, swing_list, min_uptake_list, max_uptake_list, ...
        min_RH_list, max_RH_list, drh_list, ...
        'VariableNames', {'Salt', 'Swing_pct', 'MinUptake_pct', 'MaxUptake_pct', ...
        'MinRH_pct', 'MaxRH_pct', 'DRH_pct'});
    
    writetable(pure_table, output_file_pure);
    fprintf('Pure salt results saved to: %s\n', output_file_pure);
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('KEY FINDINGS\n');
fprintf('%s\n', repmat('=', 1, 80));

if ~isempty(results) && ~isempty(pure_results)
    % Find best mixture
    [max_mix_swing, idx_max_mix] = max([results.swing]);
    best_mixture = results(idx_max_mix);
    
    % Find best pure salt
    [max_pure_swing, idx_max_pure] = max([pure_results.swing]);
    best_pure = pure_results(idx_max_pure);
    
    fprintf('\nBest Mixture: %s with %.2f%% swing\n', best_mixture.name, max_mix_swing*100);
    fprintf('Best Pure Salt: %s with %.2f%% swing\n', best_pure.name, max_pure_swing*100);
    
    if max_mix_swing > max_pure_swing
        improvement = ((max_mix_swing - max_pure_swing) / max_pure_swing) * 100;
        fprintf('=> Best mixture is %.1f%% better than best pure salt!\n', improvement);
    else
        fprintf('=> Best pure salt outperforms mixtures\n');
    end
    
    % Count how many mixtures beat their component averages
    n_better = 0;
    n_worse = 0;
    for i = 1:length(results)
        parts = strsplit(results(i).name, '+');
        pure_swings = [];
        for p = 1:length(parts)
            idx = find(strcmp({pure_results.name}, parts{p}));
            if ~isempty(idx)
                pure_swings(end+1) = pure_results(idx).swing;
            end
        end
        if ~isempty(pure_swings)
            if results(i).swing > mean(pure_swings)
                n_better = n_better + 1;
            else
                n_worse = n_worse + 1;
            end
        end
    end
    
    fprintf('\nMixture Performance Summary:\n');
    fprintf('  %d mixtures outperform their component average\n', n_better);
    fprintf('  %d mixtures underperform their component average\n', n_worse);
    fprintf('  Success rate: %.1f%%\n', n_better/(n_better+n_worse)*100);
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Calculation complete!\n');
fprintf('%s\n', repmat('=', 1, 80));

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
    charges('Cs+') = 1;
    charges('Mg++') = 2;
    charges('Ca++') = 2;
    charges('Cl-') = -1;
    charges('Br-') = -1;
    charges('I-') = -1;
    charges('OH-') = -1;
    charges('SO4--') = -2;
end
