close all
clear
clc

% Script: evaluate_all_salts_atacama.m
% Evaluate all salts with temperature-dependent Pitzer parameters
% to find the best candidates for AWH in Atacama desert conditions
%
% Metric: Maximum water uptake swing (max - min) over daily cycle
% Water uptake = kg H2O / (kg H2O + kg salt)

% Add necessary paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

% Define output directory
fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% 1. Load Pitzer Parameters Database
fprintf('Loading Pitzer parameters database...\n');
csv_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_binary.csv');
pitzer_data = readtable(csv_file);

%% 2. Load Atacama Climate Data
fprintf('Loading Atacama desert climate data...\n');
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

fprintf('  Loaded %d Atacama data points\n', length(atacama_RH));
fprintf('  RH range: %.1f%% - %.1f%%\n', min(atacama_RH)*100, max(atacama_RH)*100);
fprintf('  Temp range: %.1f°C - %.1f°C\n', min(atacama_T), max(atacama_T));

% Get maximum RH in Atacama for DRH screening
atacama_max_RH = max(atacama_RH);
fprintf('  Maximum RH available: %.1f%%\n\n', atacama_max_RH*100);

%% 3. Define Salt Properties (MW, stoichiometry, DRH)
% Molecular weights (g/mol), stoichiometry, and deliquescence RH
salt_props = containers.Map();

% Format: {MW, nu, n_cations, n_anions, charge_cation, charge_anion, DRH}
salt_props('Li+_Cl-') = {42.4, 2, 1, 1, 1, 1, 0.12};  % LiCl
salt_props('Na+_Cl-') = {58.443, 2, 1, 1, 1, 1, 0.765};  % NaCl
salt_props('K+_Cl-') = {74.551, 2, 1, 1, 1, 1, 0.855};  % KCl
salt_props('Mg++_Cl-') = {95.2, 3, 1, 2, 2, 1, 0.33};  % MgCl2
salt_props('Ca++_Cl-') = {110.98, 3, 1, 2, 2, 1, 0.31};  % CaCl2
salt_props('Li+_Br-') = {86.85, 2, 1, 1, 1, 1, 0.07};  % LiBr
salt_props('Na+_Br-') = {102.89, 2, 1, 1, 1, 1, 0.614};  % NaBr
salt_props('K+_Br-') = {119.00, 2, 1, 1, 1, 1, 0.833};  % KBr
salt_props('Na+_I-') = {149.89, 2, 1, 1, 1, 1, 0.581};  % NaI
salt_props('K+_I-') = {166.00, 2, 1, 1, 1, 1, 0.97};  % KI
salt_props('H+_Cl-') = {36.46, 2, 1, 1, 1, 1, 0.17};  % HCl
salt_props('NH4+_Cl-') = {53.49, 2, 1, 1, 1, 1, 0.815};  % NH4Cl
salt_props('Cs+_Cl-') = {168.36, 2, 1, 1, 1, 1, 0.82};  % CsCl
salt_props('Cs+_Br-') = {212.81, 2, 1, 1, 1, 1, 0.848};  % CsBr
salt_props('Sr++_Cl-') = {158.53, 3, 1, 2, 2, 1, 0.50};  % SrCl2 (estimated)
salt_props('Ba++_Cl-') = {208.23, 3, 1, 2, 2, 1, 0.90};  % BaCl2 (estimated)
salt_props('Zn++_Cl-') = {136.3, 3, 1, 2, 2, 1, 0.07};  % ZnCl2
salt_props('Na+_SO4--') = {142.04, 3, 2, 1, 1, 2, 0.93};  % Na2SO4
salt_props('K+_SO4--') = {174.26, 3, 2, 1, 1, 2, 0.97};  % K2SO4
salt_props('Mg++_SO4--') = {120.37, 2, 1, 1, 2, 2, 0.86};  % MgSO4 (estimated)
salt_props('Ca++_SO4--') = {136.14, 2, 1, 1, 2, 2, 0.98};  % CaSO4 (estimated, low solubility)
salt_props('Zn++_SO4--') = {161.47, 2, 1, 1, 2, 2, 0.90};  % ZnSO4 (estimated)
salt_props('Cu++_SO4--') = {159.61, 2, 1, 1, 2, 2, 0.97};  % CuSO4 (estimated)

% Constants
MW_water = 18.015;
Tr = 298.15;

%% 4. Screen All Salts
fprintf('Screening salts with temperature-dependent parameters...\n');
fprintf('%s\n', repmat('=', 1, 80));

results = struct('salt_name', {}, 'cation', {}, 'anion', {}, ...
    'water_uptake_min', {}, 'water_uptake_max', {}, 'swing', {}, ...
    'n_valid_points', {}, 'MW', {}, 'nu', {});

salt_count = 0;
n_screened_drh = 0;  % Track salts screened due to high DRH
n_screened_nodata = 0;  % Track salts with insufficient data
n_screened_notemp = 0;  % Track salts without temperature dependence

for row = 1:height(pitzer_data)
    species1 = pitzer_data.species1{row};
    species2 = pitzer_data.species2{row};
    
    % Check if this is a valid cation-anion pair
    % Cations: +, ++, +++, etc. or H+, NH4+
    % Anions: -, --, ---, etc.
    is_cation_1 = contains(species1, '+') && ~contains(species1, '-');
    is_anion_1 = contains(species1, '-') && ~contains(species1, '+');
    is_cation_2 = contains(species2, '+') && ~contains(species2, '-');
    is_anion_2 = contains(species2, '-') && ~contains(species2, '+');
    
    if ~((is_cation_1 && is_anion_2) || (is_cation_2 && is_anion_1))
        continue;  % Not a cation-anion pair
    end
    
    % Determine cation and anion
    if is_cation_1
        cation = species1;
        anion = species2;
    else
        cation = species2;
        anion = species1;
    end
    
    % Check if has temperature dependence (non-zero a2, a3, or a4)
    has_temp_dep = false;
    if ~isnan(pitzer_data.beta0_a2(row)) && pitzer_data.beta0_a2(row) ~= 0
        has_temp_dep = true;
    end
    if ~isnan(pitzer_data.beta0_a3(row)) && pitzer_data.beta0_a3(row) ~= 0
        has_temp_dep = true;
    end
    if ~isnan(pitzer_data.beta0_a4(row)) && pitzer_data.beta0_a4(row) ~= 0
        has_temp_dep = true;
    end
    
    if ~has_temp_dep
        n_screened_notemp = n_screened_notemp + 1;
        continue;  % No temperature dependence
    end
    
    % Create salt key
    salt_key = sprintf('%s_%s', cation, anion);
    
    % Check if we have properties for this salt
    if ~isKey(salt_props, salt_key)
        continue;  % Unknown salt
    end
    
    props = salt_props(salt_key);
    MW = props{1};
    nu = props{2};
    DRH = props{7};  % Deliquescence RH
    
    % Check if salt can deliquesce at Atacama conditions
    if DRH > atacama_max_RH
        fprintf('  Skipping %s: DRH=%.1f%% exceeds max Atacama RH=%.1f%%\n', ...
            salt_key, DRH*100, atacama_max_RH*100);
        n_screened_drh = n_screened_drh + 1;
        continue;  % Cannot deliquesce at Atacama conditions
    end
    
    salt_count = salt_count + 1;
    fprintf('\n[%d] Processing %s (MW=%.1f, nu=%d, DRH=%.1f%%)...\n', ...
        salt_count, salt_key, MW, nu, DRH*100);
    
    % Extract Pitzer parameters
    params.beta0_a1 = pitzer_data.beta0_a1(row);
    params.beta0_a2 = pitzer_data.beta0_a2(row);
    params.beta0_a3 = pitzer_data.beta0_a3(row);
    params.beta0_a4 = pitzer_data.beta0_a4(row);
    
    params.beta1_a1 = pitzer_data.beta1_a1(row);
    params.beta1_a2 = pitzer_data.beta1_a2(row);
    params.beta1_a3 = pitzer_data.beta1_a3(row);
    params.beta1_a4 = pitzer_data.beta1_a4(row);
    
    params.beta2_a1 = pitzer_data.beta2_a1(row);
    params.beta2_a2 = pitzer_data.beta2_a2(row);
    params.beta2_a3 = pitzer_data.beta2_a3(row);
    params.beta2_a4 = pitzer_data.beta2_a4(row);
    
    params.cphi_a1 = pitzer_data.cphi_a1(row);
    params.cphi_a2 = pitzer_data.cphi_a2(row);
    params.cphi_a3 = pitzer_data.cphi_a3(row);
    params.cphi_a4 = pitzer_data.cphi_a4(row);
    
    params.alpha1 = pitzer_data.alpha1(row);
    params.alpha2 = pitzer_data.alpha2(row);
    
    % Calculate water uptake for all Atacama points
    water_uptake_atacama = nan(size(atacama_RH));
    n_below_drh = 0;
    
    for i = 1:length(atacama_RH)
        RH = atacama_RH(i);
        T_C = atacama_T(i);
        T_K = T_C + 273.15;
        
        % Check if this point is above DRH
        if RH < DRH
            n_below_drh = n_below_drh + 1;
            continue;  % Below deliquescence point, no solution forms
        end
        
        % Calculate temperature-dependent Pitzer parameters
        beta0 = calc_temp_param(params.beta0_a1, params.beta0_a2, params.beta0_a3, params.beta0_a4, T_K, Tr);
        beta1 = calc_temp_param(params.beta1_a1, params.beta1_a2, params.beta1_a3, params.beta1_a4, T_K, Tr);
        beta2 = calc_temp_param(params.beta2_a1, params.beta2_a2, params.beta2_a3, params.beta2_a4, T_K, Tr);
        cphi  = calc_temp_param(params.cphi_a1, params.cphi_a2, params.cphi_a3, params.cphi_a4, T_K, Tr);
        
        % Solve for molality
        m_solution = solve_for_molality(RH, beta0, beta1, beta2, cphi, params.alpha1, params.alpha2, nu, T_C);
        
        if ~isnan(m_solution) && m_solution > 0
            % Calculate water uptake = mass fraction of water
            mass_salt_g = m_solution * MW;
            mass_water_g = 1000;
            mf_water = mass_water_g / (mass_salt_g + mass_water_g);
            water_uptake_atacama(i) = mf_water;
        end
    end
    
    if n_below_drh > 0
        fprintf('  Points below DRH: %d/%d (%.1f%%)\n', n_below_drh, length(atacama_RH), ...
            n_below_drh/length(atacama_RH)*100);
    end
    
    % Calculate swing
    valid_uptake = water_uptake_atacama(~isnan(water_uptake_atacama));
    n_valid = length(valid_uptake);
    n_above_drh = length(atacama_RH) - n_below_drh;
    
    if n_valid >= 5  % Need at least 5 points
        uptake_min = min(valid_uptake);
        uptake_max = max(valid_uptake);
        swing = uptake_max - uptake_min;
        
        fprintf('  Valid points: %d/%d above DRH (%.1f%% convergence)\n', ...
            n_valid, n_above_drh, n_valid/n_above_drh*100);
        fprintf('  Water uptake range: %.4f - %.4f\n', uptake_min, uptake_max);
        fprintf('  SWING: %.4f (%.1f%% range)\n', swing, swing*100);
        
        % Store results
        results(end+1).salt_name = strrep(salt_key, '_', ' → ');
        results(end).cation = cation;
        results(end).anion = anion;
        results(end).water_uptake_min = uptake_min;
        results(end).water_uptake_max = uptake_max;
        results(end).swing = swing;
        results(end).n_valid_points = n_valid;
        results(end).MW = MW;
        results(end).nu = nu;
    else
        fprintf('  REJECTED: Insufficient valid points (%d/%d above DRH)\n', n_valid, n_above_drh);
        n_screened_nodata = n_screened_nodata + 1;
    end
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Evaluation complete!\n');
fprintf('  Salts meeting criteria: %d\n', length(results));
fprintf('  Salts screened out:\n');
fprintf('    - High DRH (cannot deliquesce at Atacama): %d\n', n_screened_drh);
fprintf('    - Insufficient valid data points: %d\n', n_screened_nodata);
fprintf('    - No temperature dependence: %d\n', n_screened_notemp);
fprintf('  Total salts processed: %d\n\n', salt_count);

%% 5. Sort and Display Results
if isempty(results)
    error('No salts with sufficient data found!');
end

% Sort by swing (descending)
[~, sort_idx] = sort([results.swing], 'descend');
results = results(sort_idx);

% Display top salts
fprintf('\n=== TOP SALTS FOR ATACAMA AWH (by water uptake swing) ===\n\n');
fprintf('Rank  Salt               Swing   Min     Max     Points  MW      Nu\n');
fprintf('----  ----------------  ------  ------  ------  ------  ------  ---\n');

n_display = min(20, length(results));
for i = 1:n_display
    fprintf('%2d    %-16s  %.4f  %.4f  %.4f  %3d/%2d  %6.1f  %d\n', ...
        i, results(i).salt_name, results(i).swing, ...
        results(i).water_uptake_min, results(i).water_uptake_max, ...
        results(i).n_valid_points, length(atacama_RH), ...
        results(i).MW, results(i).nu);
end

%% 6. Create Visualizations

% Figure 1: Bar chart of water uptake swing
figure('Position', [100, 100, 1200, 800]);
n_plot = min(15, length(results));
swings = [results(1:n_plot).swing];
salt_names = {results(1:n_plot).salt_name};

% Clean up salt names for plotting
salt_names_clean = cellfun(@(x) strrep(x, ' → ', '-'), salt_names, 'UniformOutput', false);

barh(swings);
set(gca, 'YTick', 1:n_plot, 'YTickLabel', salt_names_clean);
xlabel('Water Uptake Swing (kg H_2O / total kg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Salt', 'FontSize', 12, 'FontWeight', 'bold');
title('Top Salts for AWH in Atacama Desert (by Daily Water Uptake Swing)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10);
xlim([0, max(swings)*1.1]);

% Add value labels
for i = 1:n_plot
    text(swings(i) + max(swings)*0.02, i, sprintf('%.3f', swings(i)), ...
        'FontSize', 9, 'VerticalAlignment', 'middle');
end

print(fullfile(fig_out_dir, 'Atacama_Salt_Comparison_Swing'), '-dpng', '-r300');

% Figure 2: Scatter plot - Swing vs MW
figure('Position', [150, 150, 1000, 700]);
MWs = [results.MW];
swings_all = [results.swing];
n_points_all = [results.n_valid_points];

scatter(MWs, swings_all, 100, n_points_all, 'filled', 'MarkerEdgeColor', 'k');
colormap(jet);
cb = colorbar;
cb.Label.String = 'Valid Atacama Points';
cb.Label.FontSize = 11;

xlabel('Molecular Weight (g/mol)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Water Uptake Swing', 'FontSize', 12, 'FontWeight', 'bold');
title('Salt Performance vs Molecular Weight', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Label top 5 salts
for i = 1:min(5, length(results))
    text(results(i).MW, results(i).swing, ['  ' salt_names_clean{i}], ...
        'FontSize', 9, 'FontWeight', 'bold');
end

print(fullfile(fig_out_dir, 'Atacama_Salt_Swing_vs_MW'), '-dpng', '-r300');

% Figure 3: Range comparison (min-max bars)
figure('Position', [200, 200, 1200, 800]);
n_plot = min(15, length(results));

mins = [results(1:n_plot).water_uptake_min];
maxs = [results(1:n_plot).water_uptake_max];
midpoints = (mins + maxs) / 2;
ranges = maxs - mins;

barh(midpoints, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
hold on;

% Add error bars showing range
for i = 1:n_plot
    plot([mins(i), maxs(i)], [i, i], 'k-', 'LineWidth', 2);
    plot(mins(i), i, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(maxs(i), i, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
end

set(gca, 'YTick', 1:n_plot, 'YTickLabel', salt_names_clean);
xlabel('Water Uptake (kg H_2O / total kg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Salt', 'FontSize', 12, 'FontWeight', 'bold');
title('Water Uptake Range: Night (Blue) to Day (Red) in Atacama', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10);
xlim([0, 1]);

legend({'Midpoint', 'Daily Range', 'Minimum (Day)', 'Maximum (Night)'}, 'Location', 'southeast');

print(fullfile(fig_out_dir, 'Atacama_Salt_Range_Comparison'), '-dpng', '-r300');

% Figure 4: Coverage map - valid points
figure('Position', [250, 250, 1000, 700]);
coverage = [results.n_valid_points] / length(atacama_RH) * 100;

bar(coverage);
set(gca, 'XTick', 1:length(results), 'XTickLabel', salt_names_clean, 'XTickLabelRotation', 45);
ylabel('Coverage (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Atacama Condition Coverage (% of valid points)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([0, 100]);
set(gca, 'FontSize', 9);

print(fullfile(fig_out_dir, 'Atacama_Salt_Coverage'), '-dpng', '-r300');

fprintf('\nAll figures saved to: %s\n', fig_out_dir);

% Save results to CSV
csv_out = fullfile(fig_out_dir, 'atacama_salt_performance.csv');
T = struct2table(results);
writetable(T, csv_out);
fprintf('Results saved to: %s\n', csv_out);

%% Helper Functions

function param_val = calc_temp_param(a1, a2, a3, a4, T, Tr)
    if isnan(a1), a1 = 0; end
    if isnan(a2), a2 = 0; end
    if isnan(a3), a3 = 0; end
    if isnan(a4), a4 = 0; end
    param_val = a1 + a2*(1/T - 1/Tr) + a3*log(T/Tr) + a4*(T - Tr);
end

function m_solution = solve_for_molality(RH, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T_C)
    % Solve for molality with robust solver
    m_solution = NaN;
    
    % Initial guess based on RH
    if RH < 0.2
        m_guess = 15;
    elseif RH < 0.4
        m_guess = 8;
    elseif RH < 0.6
        m_guess = 4;
    elseif RH < 0.8
        m_guess = 2;
    else
        m_guess = 0.5;
    end
    
    % Try bracket search first
    try
        m_solution = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T_C), ...
            [0.01, 25], optimset('Display', 'off', 'TolX', 1e-5));
        if m_solution < 0, m_solution = NaN; end
    catch
        % Try from initial guess
        try
            m_solution = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T_C), ...
                m_guess, optimset('Display', 'off', 'TolX', 1e-5));
            if m_solution < 0 || m_solution > 30, m_solution = NaN; end
        catch
            % Last resort: multiple starting points
            for m_try = [0.1, 1, 5, 10, 15, 20]
                try
                    m_test = fzero(@(m) pitzer_residual(m, RH, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T_C), ...
                        m_try, optimset('Display', 'off', 'TolX', 1e-5));
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
end

function residual = pitzer_residual(m, RH_target, beta0, beta1, beta2, cphi, alpha1, alpha2, nu, T)
    if m <= 0
        residual = 1e6;
        return;
    end
    aw_calc = pitzer_water_activity(m, nu, beta0, beta1, beta2, cphi, alpha1, alpha2, T);
    residual = aw_calc - RH_target;
end
