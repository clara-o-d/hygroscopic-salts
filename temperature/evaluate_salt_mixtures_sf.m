close all
clear
clc

% Script: evaluate_salt_mixtures_atacama.m
% Identify salt mixtures with temperature-dependent Pitzer mixing parameters
% that can deliquesce at Atacama conditions
%
% This script screens for viable mixtures based on:
% 1. Component salts must deliquesce at max Atacama RH
% 2. Binary Pitzer parameters must have temperature dependence
% 3. Ternary (psi) parameters must exist and have temperature dependence
%
% Note: Full multi-component Pitzer calculations are complex and would require
% a separate implementation. This script identifies candidates for evaluation.

% Add necessary paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'pitzer'));

%% 1. Load Atacama Climate Data to Get Max RH
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

atacama_max_RH = max(atacama_RH);
atacama_min_RH = min(atacama_RH);

fprintf('  RH range: %.1f%% - %.1f%%\n', atacama_min_RH*100, atacama_max_RH*100);
fprintf('  Temp range: %.1f°C - %.1f°C\n\n', min(atacama_T), max(atacama_T));

%% 2. Load Pitzer Binary Parameters
fprintf('Loading binary Pitzer parameters...\n');
csv_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_binary.csv');
pitzer_binary = readtable(csv_file);

%% 3. Load Pitzer Ternary (Psi) Parameters
fprintf('Loading ternary (psi) Pitzer parameters...\n');
psi_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_psi.csv');
pitzer_psi_base = readtable(psi_file);

% Load Lassin 2018 supplemental data
psi_lassin_file = fullfile(filepath, '..', 'data', 'pitzer_psi_lassin2018.csv');
pitzer_psi_lassin = readtable(psi_lassin_file);

% Merge tables (Lassin supplements existing)
pitzer_psi = [pitzer_psi_base; pitzer_psi_lassin];
fprintf('  Total psi parameters: %d (base) + %d (Lassin 2018) = %d\n\n', ...
    height(pitzer_psi_base), height(pitzer_psi_lassin), height(pitzer_psi));

%% 4. Define Salt Properties (MW, stoichiometry, DRH)
salt_props = containers.Map();

% Format: {MW, nu, DRH, full_name}
salt_props('Li+_Cl-') = {42.4, 2, 0.12, 'LiCl'};
salt_props('Na+_Cl-') = {58.443, 2, 0.765, 'NaCl'};
salt_props('K+_Cl-') = {74.551, 2, 0.855, 'KCl'};
salt_props('Mg++_Cl-') = {95.2, 3, 0.33, 'MgCl2'};
salt_props('Ca++_Cl-') = {110.98, 3, 0.31, 'CaCl2'};
salt_props('Li+_Br-') = {86.85, 2, 0.07, 'LiBr'};
salt_props('Na+_Br-') = {102.89, 2, 0.614, 'NaBr'};
salt_props('K+_Br-') = {119.00, 2, 0.833, 'KBr'};
salt_props('Cs+_Cl-') = {168.36, 2, 0.82, 'CsCl'};
salt_props('Zn++_Cl-') = {136.3, 3, 0.07, 'ZnCl2'};
salt_props('Na+_SO4--') = {142.04, 3, 0.93, 'Na2SO4'};
salt_props('K+_SO4--') = {174.26, 3, 0.97, 'K2SO4'};
salt_props('Mg++_SO4--') = {120.37, 2, 0.86, 'MgSO4'};
salt_props('Na+_OH-') = {40.00, 2, 0.23, 'NaOH'};
salt_props('H+_Cl-') = {36.46, 2, 0.17, 'HCl'};

% Constants
MW_water = 18.015;
Tr = 298.15;

%% 5. Screen Psi Parameters with Temperature Dependence
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Screening ternary mixtures with temperature-dependent psi parameters...\n');
fprintf('%s\n', repmat('=', 1, 80));

viable_mixtures = struct('species1', {}, 'species2', {}, 'species3', {}, ...
    'salt1', {}, 'salt2', {}, 'drh1', {}, 'drh2', {}, ...
    'psi_temp_dep', {}, 'salt1_temp_dep', {}, 'salt2_temp_dep', {}, ...
    'reason_rejected', {});

n_total = 0;
n_temp_dep_psi = 0;
n_atacama_compatible = 0;
n_binary_temp_dep = 0;

for row = 1:height(pitzer_psi)
    n_total = n_total + 1;
    
    species1 = pitzer_psi.species1{row};
    species2 = pitzer_psi.species2{row};
    species3 = pitzer_psi.species3{row};
    
    % Check if psi has temperature dependence
    has_psi_temp = false;
    if ~isnan(pitzer_psi.psi_a2(row)) && abs(pitzer_psi.psi_a2(row)) > 1e-15
        has_psi_temp = true;
    end
    if ~isnan(pitzer_psi.psi_a3(row)) && abs(pitzer_psi.psi_a3(row)) > 1e-15
        has_psi_temp = true;
    end
    if ~isnan(pitzer_psi.psi_a4(row)) && abs(pitzer_psi.psi_a4(row)) > 1e-15
        has_psi_temp = true;
    end
    
    if ~has_psi_temp
        continue;
    end
    
    n_temp_dep_psi = n_temp_dep_psi + 1;
    
    % Identify which are cations and which are anions
    % Psi parameters are for: cation-cation-anion or anion-anion-cation interactions
    is_cation_1 = contains(species1, '+') && ~contains(species1, '-');
    is_anion_1 = contains(species1, '-') && ~contains(species1, '+');
    is_cation_2 = contains(species2, '+') && ~contains(species2, '-');
    is_anion_2 = contains(species2, '-') && ~contains(species2, '+');
    is_cation_3 = contains(species3, '+') && ~contains(species3, '-');
    is_anion_3 = contains(species3, '-') && ~contains(species3, '+');
    
    % Determine mixture type
    cations = {};
    anion = '';
    
    if is_cation_1 && is_cation_2 && is_anion_3
        % Two cations, one anion: M1-M2-X
        cations = {species1, species2};
        anion = species3;
    elseif is_cation_1 && is_anion_2 && is_cation_3
        % M1-X-M2 format
        cations = {species1, species3};
        anion = species2;
    elseif is_anion_1 && is_cation_2 && is_anion_3
        % X1-M-X2 format (two anions, one cation) - not typical for simple salts
        continue;
    elseif is_anion_1 && is_anion_2 && is_cation_3
        % X1-X2-M format (two anions, one cation)
        continue;
    else
        % Other configurations, skip
        continue;
    end
    
    if isempty(cations) || isempty(anion)
        continue;
    end
    
    % Create salt keys
    salt1_key = sprintf('%s_%s', cations{1}, anion);
    salt2_key = sprintf('%s_%s', cations{2}, anion);
    
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
    
    % Check if both salts can deliquesce at Atacama conditions
    reason_rejected = '';
    
    if drh1 > atacama_max_RH || drh2 > atacama_max_RH
        if drh1 > atacama_max_RH && drh2 > atacama_max_RH
            reason_rejected = sprintf('Both salts: DRH1=%.1f%%, DRH2=%.1f%% > max RH=%.1f%%', ...
                drh1*100, drh2*100, atacama_max_RH*100);
        elseif drh1 > atacama_max_RH
            reason_rejected = sprintf('%s: DRH=%.1f%% > max RH=%.1f%%', ...
                salt1_name, drh1*100, atacama_max_RH*100);
        else
            reason_rejected = sprintf('%s: DRH=%.1f%% > max RH=%.1f%%', ...
                salt2_name, drh2*100, atacama_max_RH*100);
        end
    else
        n_atacama_compatible = n_atacama_compatible + 1;
    end
    
    % Check if both salts have temperature-dependent binary parameters
    salt1_has_temp = check_salt_temp_dependence(pitzer_binary, cations{1}, anion);
    salt2_has_temp = check_salt_temp_dependence(pitzer_binary, cations{2}, anion);
    
    if isempty(reason_rejected) && (~salt1_has_temp || ~salt2_has_temp)
        if ~salt1_has_temp && ~salt2_has_temp
            reason_rejected = 'Both salts lack temp-dependent binary parameters';
        elseif ~salt1_has_temp
            reason_rejected = sprintf('%s lacks temp-dependent binary parameters', salt1_name);
        else
            reason_rejected = sprintf('%s lacks temp-dependent binary parameters', salt2_name);
        end
    end
    
    if isempty(reason_rejected)
        n_binary_temp_dep = n_binary_temp_dep + 1;
    end
    
    % Store result
    viable_mixtures(end+1).species1 = species1;
    viable_mixtures(end).species2 = species2;
    viable_mixtures(end).species3 = species3;
    viable_mixtures(end).salt1 = salt1_name;
    viable_mixtures(end).salt2 = salt2_name;
    viable_mixtures(end).drh1 = drh1;
    viable_mixtures(end).drh2 = drh2;
    viable_mixtures(end).psi_temp_dep = has_psi_temp;
    viable_mixtures(end).salt1_temp_dep = salt1_has_temp;
    viable_mixtures(end).salt2_temp_dep = salt2_has_temp;
    viable_mixtures(end).reason_rejected = reason_rejected;
end

%% 6. Display Results
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('SCREENING SUMMARY\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('Total psi parameters examined: %d\n', n_total);
fprintf('  With temperature dependence: %d\n', n_temp_dep_psi);
fprintf('  Atacama-compatible (DRH check): %d\n', n_atacama_compatible);
fprintf('  With all temp-dependent binary params: %d\n', n_binary_temp_dep);
fprintf('\n');

% Separate viable from rejected
viable = viable_mixtures(arrayfun(@(x) isempty(x.reason_rejected), viable_mixtures));
rejected = viable_mixtures(arrayfun(@(x) ~isempty(x.reason_rejected), viable_mixtures));

fprintf('%s\n', repmat('=', 1, 80));
fprintf('VIABLE MIXTURES FOR ATACAMA AWH (%d found)\n', length(viable));
fprintf('%s\n', repmat('=', 1, 80));

if isempty(viable)
    fprintf('No mixtures meet all criteria.\n');
else
    fprintf('\n%-15s  %-15s  DRH1    DRH2    Psi Species\n', 'Salt 1', 'Salt 2');
    fprintf('%s\n', repmat('-', 1, 80));
    
    for i = 1:length(viable)
        fprintf('%-15s  %-15s  %5.1f%%  %5.1f%%  %s-%s-%s\n', ...
            viable(i).salt1, viable(i).salt2, ...
            viable(i).drh1*100, viable(i).drh2*100, ...
            viable(i).species1, viable(i).species2, viable(i).species3);
    end
    
    fprintf('\n');
    fprintf('NOTE: These mixtures have all required temperature-dependent parameters\n');
    fprintf('      and can deliquesce at Atacama conditions. Full multi-component\n');
    fprintf('      Pitzer calculations would be needed to evaluate water uptake swing.\n');
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('REJECTED MIXTURES (showing up to 20)\n');
fprintf('%s\n', repmat('=', 1, 80));

n_show = min(20, length(rejected));
for i = 1:n_show
    fprintf('\n%-12s + %-12s (DRH: %.1f%%, %.1f%%)\n', ...
        rejected(i).salt1, rejected(i).salt2, ...
        rejected(i).drh1*100, rejected(i).drh2*100);
    fprintf('  Reason: %s\n', rejected(i).reason_rejected);
end

if length(rejected) > 20
    fprintf('\n... and %d more rejected mixtures\n', length(rejected) - 20);
end

%% 7. Save Results
output_file = fullfile(filepath, 'atacama_mixture_candidates.csv');
if ~isempty(viable)
    % Create table for export
    salt1_list = {viable.salt1}';
    salt2_list = {viable.salt2}';
    drh1_list = [viable.drh1]';
    drh2_list = [viable.drh2]';
    species1_list = {viable.species1}';
    species2_list = {viable.species2}';
    species3_list = {viable.species3}';
    
    result_table = table(salt1_list, salt2_list, drh1_list, drh2_list, ...
        species1_list, species2_list, species3_list, ...
        'VariableNames', {'Salt1', 'Salt2', 'DRH1', 'DRH2', ...
        'Species1', 'Species2', 'Species3'});
    
    writetable(result_table, output_file);
    fprintf('\nViable mixture candidates saved to: %s\n', output_file);
end

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Screening complete!\n');
fprintf('%s\n', repmat('=', 1, 80));

%% Helper Function

function has_temp = check_salt_temp_dependence(pitzer_binary, cation, anion)
    % Check if a salt has temperature-dependent binary Pitzer parameters
    
    % Find the row for this cation-anion pair
    idx = find((strcmp(pitzer_binary.species1, cation) & strcmp(pitzer_binary.species2, anion)) | ...
               (strcmp(pitzer_binary.species1, anion) & strcmp(pitzer_binary.species2, cation)));
    
    if isempty(idx)
        has_temp = false;
        return;
    end
    
    % Check if any beta0, beta1, or cphi has temperature dependence (a2, a3, or a4)
    has_temp = false;
    
    % Check beta0
    if ~isnan(pitzer_binary.beta0_a2(idx)) && abs(pitzer_binary.beta0_a2(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.beta0_a3(idx)) && abs(pitzer_binary.beta0_a3(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.beta0_a4(idx)) && abs(pitzer_binary.beta0_a4(idx)) > 1e-15
        has_temp = true;
    end
    
    % Check beta1
    if ~isnan(pitzer_binary.beta1_a2(idx)) && abs(pitzer_binary.beta1_a2(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.beta1_a3(idx)) && abs(pitzer_binary.beta1_a3(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.beta1_a4(idx)) && abs(pitzer_binary.beta1_a4(idx)) > 1e-15
        has_temp = true;
    end
    
    % Check cphi
    if ~isnan(pitzer_binary.cphi_a2(idx)) && abs(pitzer_binary.cphi_a2(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.cphi_a3(idx)) && abs(pitzer_binary.cphi_a3(idx)) > 1e-15
        has_temp = true;
    end
    if ~isnan(pitzer_binary.cphi_a4(idx)) && abs(pitzer_binary.cphi_a4(idx)) > 1e-15
        has_temp = true;
    end
end
