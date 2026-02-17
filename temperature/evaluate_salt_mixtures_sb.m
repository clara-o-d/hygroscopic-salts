close all
clear
clc

% Script: evaluate_salt_mixtures_sb.m
% Identify salt mixtures with Pitzer mixing parameters that can deliquesce 
% at Santa Barbara conditions (RH: 42-88%)
%
% This script screens for viable mixtures based on:
% 1. Component salts must deliquesce at max Santa Barbara RH
% 2. Ternary (psi) mixing parameters must exist
% 3. Binary Pitzer parameters must exist
%
% Temperature dependence is not required - uses 25°C values

% Add necessary paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

%% 1. Define Santa Barbara Climate Conditions
fprintf('Santa Barbara climate conditions...\n');

% Santa Barbara RH typically ranges from 42% to 88%
sb_min_RH = 0.42;
sb_max_RH = 0.88;
sb_avg_T = 16;  % Average ~61°F = 16°C

fprintf('  RH range: %.1f%% - %.1f%%\n', sb_min_RH*100, sb_max_RH*100);
fprintf('  Average temp: %.1f°C\n', sb_avg_T);
fprintf('  Using max RH = %.1f%% as DRH threshold\n\n', sb_max_RH*100);

%% 2. Load Pitzer Binary Parameters
fprintf('Loading binary Pitzer parameters...\n');
csv_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_binary.csv');
pitzer_binary = readtable(csv_file);
fprintf('  Total binary parameters: %d\n', height(pitzer_binary));

%% 3. Load Pitzer Ternary (Psi) Parameters - UPDATED VERSION
fprintf('Loading ternary (psi) Pitzer parameters (updated)...\n');
psi_file = fullfile(filepath, '..', 'data', 'pitzer_psi_updated.csv');
pitzer_psi_base = readtable(psi_file);

% Load Lassin 2018 supplemental data
psi_lassin_file = fullfile(filepath, '..', 'data', 'pitzer_psi_lassin2018.csv');
pitzer_psi_lassin = readtable(psi_lassin_file);

% Add missing columns to lassin table to match base table structure
if ~ismember('psi_a2', pitzer_psi_lassin.Properties.VariableNames)
    pitzer_psi_lassin.psi_a2 = NaN(height(pitzer_psi_lassin), 1);
end
if ~ismember('psi_a3', pitzer_psi_lassin.Properties.VariableNames)
    pitzer_psi_lassin.psi_a3 = NaN(height(pitzer_psi_lassin), 1);
end
if ~ismember('psi_a4', pitzer_psi_lassin.Properties.VariableNames)
    pitzer_psi_lassin.psi_a4 = NaN(height(pitzer_psi_lassin), 1);
end

% Reorder columns to match base table
pitzer_psi_lassin = pitzer_psi_lassin(:, pitzer_psi_base.Properties.VariableNames);

% Merge tables (Lassin supplements existing)
pitzer_psi = [pitzer_psi_base; pitzer_psi_lassin];
fprintf('  Total psi parameters: %d (base) + %d (Lassin 2018) = %d\n', ...
    height(pitzer_psi_base), height(pitzer_psi_lassin), height(pitzer_psi));

%% 4. Load Pitzer Theta (Like-Ion) Parameters - UPDATED VERSION
fprintf('Loading theta (like-ion) Pitzer parameters (updated)...\n');
theta_file = fullfile(filepath, '..', 'data', 'pitzer_theta_updated.csv');
pitzer_theta_base = readtable(theta_file);

% Load Lassin 2018 supplemental data
theta_lassin_file = fullfile(filepath, '..', 'data', 'pitzer_theta_lassin2018.csv');
pitzer_theta_lassin = readtable(theta_lassin_file);

% Add missing columns to lassin table to match base table structure if needed
base_vars = pitzer_theta_base.Properties.VariableNames;
for i = 1:length(base_vars)
    if ~ismember(base_vars{i}, pitzer_theta_lassin.Properties.VariableNames)
        pitzer_theta_lassin.(base_vars{i}) = NaN(height(pitzer_theta_lassin), 1);
    end
end

% Reorder columns to match base table
pitzer_theta_lassin = pitzer_theta_lassin(:, pitzer_theta_base.Properties.VariableNames);

% Merge tables (Lassin supplements existing)
pitzer_theta = [pitzer_theta_base; pitzer_theta_lassin];
fprintf('  Total theta parameters: %d (base) + %d (Lassin 2018) = %d\n\n', ...
    height(pitzer_theta_base), height(pitzer_theta_lassin), height(pitzer_theta));

%% 5. Define Salt Properties (MW, stoichiometry, DRH)
salt_props = containers.Map();

% Format: {MW, nu, DRH, full_name}
% Hygroscopic salts (low to moderate DRH - good for Santa Barbara)
salt_props('Li+_Cl-') = {42.4, 2, 0.12, 'LiCl'};
salt_props('Li+_Br-') = {86.85, 2, 0.07, 'LiBr'};
salt_props('Li+_I-') = {133.85, 2, 0.18, 'LiI'};
salt_props('Ca++_Cl-') = {110.98, 3, 0.31, 'CaCl2'};
salt_props('Mg++_Cl-') = {95.2, 3, 0.33, 'MgCl2'};
salt_props('Zn++_Cl-') = {136.3, 3, 0.07, 'ZnCl2'};
salt_props('H+_Cl-') = {36.46, 2, 0.17, 'HCl'};
salt_props('Na+_OH-') = {40.00, 2, 0.23, 'NaOH'};

% Moderately hygroscopic (work well at Santa Barbara RH)
salt_props('Na+_Br-') = {102.89, 2, 0.614, 'NaBr'};
salt_props('Na+_I-') = {149.89, 2, 0.581, 'NaI'};
salt_props('Na+_Cl-') = {58.443, 2, 0.765, 'NaCl'};
salt_props('K+_Cl-') = {74.551, 2, 0.855, 'KCl'};
salt_props('K+_Br-') = {119.00, 2, 0.833, 'KBr'};
salt_props('Cs+_Cl-') = {168.36, 2, 0.82, 'CsCl'};
salt_props('Mg++_SO4--') = {120.37, 2, 0.86, 'MgSO4'};

% Less hygroscopic (may still work at Santa Barbara's high RH)
salt_props('Na+_SO4--') = {142.04, 3, 0.93, 'Na2SO4'};
salt_props('K+_NO3-') = {101.10, 2, 0.932, 'KNO3'};
salt_props('K+_SO4--') = {174.26, 3, 0.97, 'K2SO4'};

%% 6. Screen Psi and Theta Parameters
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Screening ternary mixtures with psi and theta parameters...\n');
fprintf('%s\n', repmat('=', 1, 80));

viable_mixtures = struct('species1', {}, 'species2', {}, 'species3', {}, ...
    'salt1', {}, 'salt2', {}, 'drh1', {}, 'drh2', {}, ...
    'psi_a1', {}, 'theta_a1', {}, 'has_theta', {}, 'mixture_type', {}, 'reason_rejected', {});

n_total = 0;
n_valid_structure = 0;
n_sb_compatible = 0;
n_binary_params = 0;
n_has_theta = 0;

for row = 1:height(pitzer_psi)
    n_total = n_total + 1;
    
    species1 = pitzer_psi.species1{row};
    species2 = pitzer_psi.species2{row};
    species3 = pitzer_psi.species3{row};
    
    % Get psi value
    psi_a1 = pitzer_psi.psi_a1(row);
    
    if isnan(psi_a1)
        continue;
    end
    
    % Identify which are cations and which are anions
    is_cation_1 = contains(species1, '+') && ~contains(species1, '-');
    is_anion_1 = contains(species1, '-') && ~contains(species1, '+');
    is_cation_2 = contains(species2, '+') && ~contains(species2, '-');
    is_anion_2 = contains(species2, '-') && ~contains(species2, '+');
    is_cation_3 = contains(species3, '+') && ~contains(species3, '-');
    is_anion_3 = contains(species3, '-') && ~contains(species3, '+');
    
    % Determine mixture type
    cations = {};
    anions = {};
    cation = '';
    anion = '';
    mixture_type = '';
    
    if is_cation_1 && is_cation_2 && is_anion_3
        % Two cations, one anion: M1-M2-X
        cations = {species1, species2};
        anion = species3;
        mixture_type = sprintf('%s-%s-%s', species1, species2, species3);
    elseif is_cation_1 && is_anion_2 && is_cation_3
        % M1-X-M2 format
        cations = {species1, species3};
        anion = species2;
        mixture_type = sprintf('%s-%s-%s', species1, anion, species3);
    elseif is_anion_1 && is_cation_2 && is_anion_3
        % X1-M-X2 format (two anions, one cation)
        anions = {species1, species3};
        cation = species2;
        mixture_type = sprintf('%s-%s-%s (anion mix)', species1, cation, species3);
    elseif is_anion_1 && is_anion_2 && is_cation_3
        % X1-X2-M format (two anions, one cation)
        anions = {species1, species2};
        cation = species3;
        mixture_type = sprintf('%s-%s-%s (anion mix)', species1, species2, species3);
    else
        % Other configurations, skip
        continue;
    end
    
    % Check if we have valid cation-cation or anion-anion mixture
    if ~(~isempty(cations) && length(cations) == 2 && ~isempty(anion)) && ...
       ~(~isempty(anions) && length(anions) == 2 && ~isempty(cation))
        continue;
    end
    
    n_valid_structure = n_valid_structure + 1;
    
    % Create salt keys based on mixture type
    if ~isempty(cations)
        % Cation-cation mixture (two cations, one anion)
        salt1_key = sprintf('%s_%s', cations{1}, anion);
        salt2_key = sprintf('%s_%s', cations{2}, anion);
        like_ions = cations;
    else
        % Anion-anion mixture (two anions, one cation)
        salt1_key = sprintf('%s_%s', cation, anions{1});
        salt2_key = sprintf('%s_%s', cation, anions{2});
        like_ions = anions;
    end
    
    % Check if both salts are in our database
    if ~isKey(salt_props, salt1_key) || ~isKey(salt_props, salt2_key)
        continue;
    end
    
    props1 = salt_props(salt1_key);
    props2 = salt_props(salt2_key);
    
    drh1 = props1{3};
    drh2 = props2{3};
    salt1_name = props1{4};
    salt2_name = props2{4};
    
    % Check if both salts can deliquesce at Santa Barbara conditions
    reason_rejected = '';
    
    if drh1 > sb_max_RH || drh2 > sb_max_RH
        if drh1 > sb_max_RH && drh2 > sb_max_RH
            reason_rejected = sprintf('Both salts exceed max RH: DRH1=%.1f%%, DRH2=%.1f%%', ...
                drh1*100, drh2*100);
        elseif drh1 > sb_max_RH
            reason_rejected = sprintf('%s DRH=%.1f%% exceeds max RH=%.1f%%', ...
                salt1_name, drh1*100, sb_max_RH*100);
        else
            reason_rejected = sprintf('%s DRH=%.1f%% exceeds max RH=%.1f%%', ...
                salt2_name, drh2*100, sb_max_RH*100);
        end
    else
        n_sb_compatible = n_sb_compatible + 1;
    end
    
    % Check if both salts have binary Pitzer parameters
    salt1_has_binary = false;
    salt2_has_binary = false;
    
    if ~isempty(cations)
        % Cation-cation mixture: check cat1-anion and cat2-anion
        salt1_has_binary = check_salt_binary_params(pitzer_binary, cations{1}, anion);
        salt2_has_binary = check_salt_binary_params(pitzer_binary, cations{2}, anion);
    else
        % Anion-anion mixture: check cation-an1 and cation-an2
        salt1_has_binary = check_salt_binary_params(pitzer_binary, cation, anions{1});
        salt2_has_binary = check_salt_binary_params(pitzer_binary, cation, anions{2});
    end
    
    if isempty(reason_rejected) && ((~salt1_has_binary) || (~salt2_has_binary))
        if ~salt1_has_binary && ~salt2_has_binary
            reason_rejected = 'Both salts lack binary Pitzer parameters';
        elseif ~salt1_has_binary
            reason_rejected = sprintf('%s lacks binary Pitzer parameters', salt1_name);
        else
            reason_rejected = sprintf('%s lacks binary Pitzer parameters', salt2_name);
        end
    end
    
    if isempty(reason_rejected)
        n_binary_params = n_binary_params + 1;
    end
    
    % Check if theta parameter exists for like-ion interaction
    [has_theta, theta_value] = check_theta_params(pitzer_theta, like_ions{1}, like_ions{2});
    
    if isempty(reason_rejected) && ~has_theta
        reason_rejected = sprintf('Missing theta parameter for %s-%s interaction', like_ions{1}, like_ions{2});
    end
    
    if has_theta
        n_has_theta = n_has_theta + 1;
    end
    
    % Store result
    viable_mixtures(end+1).species1 = species1;
    viable_mixtures(end).species2 = species2;
    viable_mixtures(end).species3 = species3;
    viable_mixtures(end).salt1 = salt1_name;
    viable_mixtures(end).salt2 = salt2_name;
    viable_mixtures(end).drh1 = drh1;
    viable_mixtures(end).drh2 = drh2;
    viable_mixtures(end).psi_a1 = psi_a1;
    viable_mixtures(end).theta_a1 = theta_value;
    viable_mixtures(end).has_theta = has_theta;
    viable_mixtures(end).mixture_type = mixture_type;
    viable_mixtures(end).reason_rejected = reason_rejected;
end

%% 7. Display Results
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('SCREENING SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Total psi parameters examined: %d\n', n_total);
fprintf('  Valid mixture structures (2 cations + 1 anion OR 2 anions + 1 cation): %d\n', n_valid_structure);
fprintf('  Santa Barbara-compatible (both salts DRH < %.1f%%): %d\n', sb_max_RH*100, n_sb_compatible);
fprintf('  With complete binary parameters: %d\n', n_binary_params);
fprintf('  With theta (like-ion) parameters: %d\n', n_has_theta);
fprintf('\n');

% Separate viable from rejected
viable = viable_mixtures(arrayfun(@(x) isempty(x.reason_rejected), viable_mixtures));
rejected = viable_mixtures(arrayfun(@(x) ~isempty(x.reason_rejected), viable_mixtures));

fprintf('%s\n', repmat('=', 1, 80));
fprintf('VIABLE MIXTURES FOR SANTA BARBARA AWH (%d found)\n', length(viable));
fprintf('%s\n', repmat('=', 1, 80));

if isempty(viable)
    fprintf('No mixtures meet all criteria.\n');
else
    fprintf('\n%-12s  %-12s  DRH1   DRH2   Psi      Theta    Mixture Type\n', 'Salt 1', 'Salt 2');
    fprintf('%s\n', repmat('-', 1, 85));
    
    for i = 1:length(viable)
        fprintf('%-12s  %-12s  %5.1f%% %5.1f%% %7.4f  %7.4f  %s\n', ...
            viable(i).salt1, viable(i).salt2, ...
            viable(i).drh1*100, viable(i).drh2*100, ...
            viable(i).psi_a1, viable(i).theta_a1, viable(i).mixture_type);
    end
    
    fprintf('\n');
    fprintf('NOTE: These mixtures have complete Pitzer parameters (binary, psi, and theta)\n');
    fprintf('      and can deliquesce at Santa Barbara conditions. Multi-component Pitzer\n');
    fprintf('      calculations would be needed to evaluate actual water uptake swing.\n');
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('SANTA BARBARA-COMPATIBLE MIXTURES (DRH check passed, but missing parameters)\n');
fprintf('%s\n', repmat('=', 1, 80));

% Show mixtures that pass DRH check but fail for other reasons
sb_ok = rejected(arrayfun(@(x) ~contains(x.reason_rejected, 'exceed'), rejected));
if ~isempty(sb_ok)
    fprintf('\n%-12s  %-12s  DRH1   DRH2   Issue\n', 'Salt 1', 'Salt 2');
    fprintf('%s\n', repmat('-', 1, 80));
    for i = 1:min(20, length(sb_ok))
        fprintf('%-12s  %-12s  %5.1f%% %5.1f%% %s\n', ...
            sb_ok(i).salt1, sb_ok(i).salt2, ...
            sb_ok(i).drh1*100, sb_ok(i).drh2*100, ...
            sb_ok(i).reason_rejected);
    end
    if length(sb_ok) > 20
        fprintf('... and %d more\n', length(sb_ok) - 20);
    end
else
    fprintf('None found.\n');
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('REJECTED: DRH TOO HIGH (showing up to 15)\n');
fprintf('%s\n', repmat('=', 1, 80));

% Show mixtures rejected due to high DRH
drh_rejected = rejected(arrayfun(@(x) contains(x.reason_rejected, 'exceed'), rejected));
n_show = min(15, length(drh_rejected));
if n_show > 0
    for i = 1:n_show
        fprintf('%-12s + %-12s (DRH: %.1f%%, %.1f%%) - %s\n', ...
            drh_rejected(i).salt1, drh_rejected(i).salt2, ...
            drh_rejected(i).drh1*100, drh_rejected(i).drh2*100, ...
            drh_rejected(i).reason_rejected);
    end
    
    if length(drh_rejected) > 15
        fprintf('... and %d more\n', length(drh_rejected) - 15);
    end
else
    fprintf('None - all salts can deliquesce at Santa Barbara RH levels!\n');
end

%% 8. Save Results
if ~isempty(viable)
    output_file = fullfile(filepath, 'santa_barbara_mixture_candidates.csv');
    
    % Create table for export
    salt1_list = {viable.salt1}';
    salt2_list = {viable.salt2}';
    drh1_list = [viable.drh1]';
    drh2_list = [viable.drh2]';
    psi_list = [viable.psi_a1]';
    theta_list = [viable.theta_a1]';
    mixture_type_list = {viable.mixture_type}';
    
    result_table = table(salt1_list, salt2_list, drh1_list, drh2_list, psi_list, ...
        theta_list, mixture_type_list, ...
        'VariableNames', {'Salt1', 'Salt2', 'DRH1', 'DRH2', 'Psi_a1', 'Theta_a1', 'MixtureType'});
    
    writetable(result_table, output_file);
    fprintf('\nViable mixture candidates saved to: %s\n', output_file);
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('COMMONLY REQUESTED MIXTURES\n');
fprintf('%s\n', repmat('=', 1, 80));

% Check specific mixtures people might be interested in
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Li+', 'Ca++', 'Cl-', 'LiCl', 'CaCl2');
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Li+', 'Mg++', 'Cl-', 'LiCl', 'MgCl2');
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Ca++', 'Mg++', 'Cl-', 'CaCl2', 'MgCl2');
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Li+', 'Na+', 'Cl-', 'LiCl', 'NaCl');
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Na+', 'K+', 'Cl-', 'NaCl', 'KCl');
check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, 'Na+', 'Mg++', 'Cl-', 'NaCl', 'MgCl2');

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Screening complete!\n');
fprintf('NOTE: Santa Barbara''s higher RH (42-88%%) allows many more salts than Atacama!\n');
fprintf('%s\n', repmat('=', 1, 80));

%% Helper Functions

function has_params = check_salt_binary_params(pitzer_binary, cation, anion)
    % Check if a salt has binary Pitzer parameters
    
    % Find the row for this cation-anion pair
    idx = find((strcmp(pitzer_binary.species1, cation) & strcmp(pitzer_binary.species2, anion)) | ...
               (strcmp(pitzer_binary.species1, anion) & strcmp(pitzer_binary.species2, cation)));
    
    if isempty(idx)
        has_params = false;
        return;
    end
    
    % Take first match and check if beta0 parameter exists (at minimum)
    has_params = ~isnan(pitzer_binary.beta0_a1(idx(1)));
end

function [has_theta, theta_value] = check_theta_params(pitzer_theta, ion1, ion2)
    % Check if theta parameter exists for like-ion interaction (cation-cation or anion-anion)
    
    % Find the row for this ion-ion pair (either order)
    idx = find((strcmp(pitzer_theta.species1, ion1) & strcmp(pitzer_theta.species2, ion2)) | ...
               (strcmp(pitzer_theta.species1, ion2) & strcmp(pitzer_theta.species2, ion1)));
    
    if isempty(idx)
        has_theta = false;
        theta_value = NaN;
        return;
    end
    
    % Get theta value (take first match if multiple)
    theta_value = pitzer_theta.theta_a1(idx(1));
    has_theta = ~isnan(theta_value);
end

function check_specific_mixture(pitzer_binary, pitzer_theta, pitzer_psi, cation1, cation2, anion, salt1_name, salt2_name)
    % Check if a specific mixture has all required parameters
    
    fprintf('\n%s + %s:\n', salt1_name, salt2_name);
    
    % Check binary parameters
    has_salt1 = check_salt_binary_params(pitzer_binary, cation1, anion);
    has_salt2 = check_salt_binary_params(pitzer_binary, cation2, anion);
    fprintf('  Binary params: %s (%s), %s (%s)\n', ...
        salt1_name, ifelse(has_salt1, '✓', '✗'), ...
        salt2_name, ifelse(has_salt2, '✓', '✗'));
    
    % Check theta parameter
    [has_theta, theta_val] = check_theta_params(pitzer_theta, cation1, cation2);
    if has_theta
        fprintf('  Theta %s-%s: ✓ (%.4f)\n', cation1, cation2, theta_val);
    else
        fprintf('  Theta %s-%s: ✗ MISSING\n', cation1, cation2);
    end
    
    % Check psi parameter
    psi_idx = [];
    for i = 1:height(pitzer_psi)
        s = sort({pitzer_psi.species1{i}, pitzer_psi.species2{i}, pitzer_psi.species3{i}});
        target = sort({cation1, cation2, anion});
        if isequal(s, target)
            psi_idx = i;
            break;
        end
    end
    if ~isempty(psi_idx)
        fprintf('  Psi %s-%s-%s: ✓ (%.4f)\n', cation1, cation2, anion, pitzer_psi.psi_a1(psi_idx));
    else
        fprintf('  Psi %s-%s-%s: ✗ MISSING\n', cation1, cation2, anion);
    end
    
    if has_salt1 && has_salt2 && has_theta && ~isempty(psi_idx)
        fprintf('  → STATUS: ✓ VIABLE for AWH evaluation\n');
    else
        fprintf('  → STATUS: ✗ Missing parameters\n');
    end
end

function result = ifelse(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
