close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'classification_analysis');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Define all salts with their molecular weights, valid RH ranges, and ION COUNTS
% Format: {salt_name, MW, RH_min, RH_max, function_name, function_args, num_cations, num_anions}
% function_args: 0 = no args, 1 = needs T

salt_data = {
    % Endothermic salts
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1};
    {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1};
    {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1};
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};
    
    % Exothermic salts
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2};
    {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
    {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1};
    {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2};
    {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
    {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1};
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1};
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2};
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2};
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2};
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2};
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2};
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2};
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2};
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
};

% Process each salt
num_points = 100;
all_salts_data = struct();

if isempty(salt_data)
    error('salt_data is empty!');
end

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    
    % Get Ion counts
    n_cat = salt_data{s}{7};
    n_an  = salt_data{s}{8};
    nu = n_cat + n_an; % Total number of dissociated ions
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    
    x_water_mol = zeros(size(RH_vec)); % Molecular basis
    x_water_ion = zeros(size(RH_vec)); % Ionic basis
    
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
            
            % 1. Molecular Definition: x_w = n_w / (n_w + n_s)
            x_water_mol(i) = n_w / (n_w + n_s);
            
            % 2. Ionic Definition: x_w = n_w / (n_w + nu * n_s)
            x_water_ion(i) = n_w / (n_w + (nu * n_s));
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        gamma_w_mol = RH_vec ./ x_water_mol;
        gamma_w_ion = RH_vec ./ x_water_ion;
        
        % Calculate water activity (aw = RH)
        aw = RH_vec;
        
        % Calculate logarithmic quantities
        ln_gamma_w_mol = log(gamma_w_mol);
        ln_gamma_w_ion = log(gamma_w_ion);
        ln_aw = log(aw);
        
        % Store data
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        
        all_salts_data.(salt_name).x_water_mol = x_water_mol;
        all_salts_data.(salt_name).gamma_w_mol = gamma_w_mol;
        all_salts_data.(salt_name).ln_gamma_w_mol = ln_gamma_w_mol;
        
        all_salts_data.(salt_name).x_water_ion = x_water_ion;
        all_salts_data.(salt_name).gamma_w_ion = gamma_w_ion;
        all_salts_data.(salt_name).ln_gamma_w_ion = ln_gamma_w_ion;
        
        all_salts_data.(salt_name).aw = aw;
        all_salts_data.(salt_name).ln_aw = ln_aw;
        
        all_salts_data.(salt_name).mf_water = mf_water;
        all_salts_data.(salt_name).mf_salt = mf_salt;
        all_salts_data.(salt_name).MW = MW;
        all_salts_data.(salt_name).nu = nu;
        all_salts_data.(salt_name).n_cat = n_cat;
        all_salts_data.(salt_name).n_an = n_an;
        all_salts_data.(salt_name).display_name = salt_name;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

%% ========================================================================
%% CLASSIFICATION DEFINITIONS
%% ========================================================================

% --- Classification 1: Endothermic vs Exothermic ---
endothermic_salts = {'NaCl', 'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', ...
                     'LiNO3', 'KNO3', 'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', ...
                     'RbCl', 'CsBr', 'CsI', 'Na2SO4', 'K2SO4', 'NH42SO4', ...
                     'MgSO4', 'MnSO4', 'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4', ...
                     'BaNO3', 'BaCl2'};

% --- Classification 2: Cation Kosmotrope vs Chaotrope ---
% Kosmotropes: small, highly charged ions that structure water (Li+, Na+, Mg2+, Ca2+, Al3+)
% Chaotropes: large, weakly charged ions that disrupt water structure (K+, Rb+, Cs+, NH4+)
cation_kosmotropes = {'Li', 'Na', 'Mg', 'Ca', 'Sr', 'Ba', 'Zn', 'Cu', 'Ni', 'Mn', 'Ag', 'H'};
cation_chaotropes = {'K', 'Rb', 'Cs', 'NH4'};

% --- Classification 3: Anion Kosmotrope vs Chaotrope ---
% Kosmotropes: small, highly charged anions (SO4^2-, OH-)
% Chaotropes: large, weakly charged anions (I-, Br-, NO3-, ClO4-)
anion_kosmotropes = {'SO4', 'OH'};
anion_chaotropes = {'I', 'Br', 'NO3', 'ClO4', 'ClO3'};
anion_intermediate = {'Cl'}; % Cl- is often considered intermediate

%% Helper function to extract cation from salt name
function cation = extract_cation(salt_name)
    if startsWith(salt_name, 'NH4')
        cation = 'NH4';
    elseif startsWith(salt_name, 'Li')
        cation = 'Li';
    elseif startsWith(salt_name, 'Na')
        cation = 'Na';
    elseif startsWith(salt_name, 'K')
        cation = 'K';
    elseif startsWith(salt_name, 'Rb')
        cation = 'Rb';
    elseif startsWith(salt_name, 'Cs')
        cation = 'Cs';
    elseif startsWith(salt_name, 'Mg')
        cation = 'Mg';
    elseif startsWith(salt_name, 'Ca')
        cation = 'Ca';
    elseif startsWith(salt_name, 'Sr')
        cation = 'Sr';
    elseif startsWith(salt_name, 'Ba')
        cation = 'Ba';
    elseif startsWith(salt_name, 'Zn')
        cation = 'Zn';
    elseif startsWith(salt_name, 'Cu')
        cation = 'Cu';
    elseif startsWith(salt_name, 'Ni')
        cation = 'Ni';
    elseif startsWith(salt_name, 'Mn')
        cation = 'Mn';
    elseif startsWith(salt_name, 'Ag')
        cation = 'Ag';
    elseif startsWith(salt_name, 'H')
        cation = 'H';
    else
        cation = 'Unknown';
    end
end

%% Helper function to extract anion from salt name
function anion = extract_anion(salt_name)
    if contains(salt_name, 'SO4')
        anion = 'SO4';
    elseif contains(salt_name, 'NO3')
        anion = 'NO3';
    elseif contains(salt_name, 'ClO4')
        anion = 'ClO4';
    elseif contains(salt_name, 'ClO3')
        anion = 'ClO3';
    elseif contains(salt_name, 'OH')
        anion = 'OH';
    elseif contains(salt_name, 'Cl')
        anion = 'Cl';
    elseif contains(salt_name, 'Br')
        anion = 'Br';
    elseif contains(salt_name, 'I')
        anion = 'I';
    else
        anion = 'Unknown';
    end
end

% Apply classifications to each salt
classifications = struct();

for s = 1:num_salts
    salt_name = salt_names{s};
    
    % Classification 1: Endothermic vs Exothermic
    classifications.(salt_name).thermal = 'Exothermic';
    if any(strcmp(salt_name, endothermic_salts))
        classifications.(salt_name).thermal = 'Endothermic';
    end
    
    % Extract cation and anion
    cation = extract_cation(salt_name);
    anion = extract_anion(salt_name);
    classifications.(salt_name).cation = cation;
    classifications.(salt_name).anion = anion;
    
    % Classification 2: Cation Kosmotrope vs Chaotrope
    classifications.(salt_name).cation_class = 'Unknown';
    if any(strcmp(cation, cation_kosmotropes))
        classifications.(salt_name).cation_class = 'Kosmotrope';
    elseif any(strcmp(cation, cation_chaotropes))
        classifications.(salt_name).cation_class = 'Chaotrope';
    end
    
    % Classification 3: Anion Kosmotrope vs Chaotrope
    classifications.(salt_name).anion_class = 'Unknown';
    if any(strcmp(anion, anion_kosmotropes))
        classifications.(salt_name).anion_class = 'Kosmotrope';
    elseif any(strcmp(anion, anion_chaotropes))
        classifications.(salt_name).anion_class = 'Chaotrope';
    elseif any(strcmp(anion, anion_intermediate))
        classifications.(salt_name).anion_class = 'Intermediate';
    end
    
    % Classification 4: Electrolyte Type
    n_cat = all_salts_data.(salt_name).n_cat;
    n_an = all_salts_data.(salt_name).n_an;
    classifications.(salt_name).electrolyte_type = sprintf('%d:%d', n_cat, n_an);
    
    % Classification 5: Anion Family
    if strcmp(anion, 'SO4')
        classifications.(salt_name).anion_family = 'Sulfate';
    elseif strcmp(anion, 'NO3')
        classifications.(salt_name).anion_family = 'Nitrate';
    elseif any(strcmp(anion, {'Cl', 'Br', 'I'}))
        classifications.(salt_name).anion_family = 'Halide';
    elseif any(strcmp(anion, {'ClO4', 'ClO3'}))
        classifications.(salt_name).anion_family = 'Perchlorate/Chlorate';
    elseif strcmp(anion, 'OH')
        classifications.(salt_name).anion_family = 'Hydroxide';
    else
        classifications.(salt_name).anion_family = 'Other';
    end
    
    % Classification 6: Combined Ion Classes
    cat_class = classifications.(salt_name).cation_class;
    an_class = classifications.(salt_name).anion_class;
    classifications.(salt_name).combined_ion_class = sprintf('%s-Cat/%s-An', cat_class, an_class);
end

%% ========================================================================
%% QUANTITATIVE SEPARATION ANALYSIS
%% ========================================================================

fprintf('\n========================================\n');
fprintf('QUANTITATIVE CLASSIFICATION ANALYSIS\n');
fprintf('========================================\n');

% Test multiple RH values to get robust metrics
test_RH_values = [0.50, 0.60, 0.70, 0.80, 0.90];

class_fields = {'thermal', 'cation_class', 'anion_class', 'electrolyte_type', ...
               'anion_family', 'combined_ion_class'};
class_titles = {'Endothermic vs Exothermic', 'Cation Class', 'Anion Class', ...
               'Electrolyte Type', 'Anion Family', 'Combined Ion Classes'};

% Analyze for both molecular and ionic bases
basis_types = {'molecular', 'ionic'};
basis_fields = {'ln_gamma_w_mol', 'ln_gamma_w_ion'};
basis_titles = {'Molecular Basis', 'Ionic Basis'};

% Store separation metrics for each classification and basis
separation_results_mol = struct();
separation_results_ion = struct();

for basis_idx = 1:length(basis_types)
    basis_type = basis_types{basis_idx};
    metric_field = basis_fields{basis_idx};
    basis_title = basis_titles{basis_idx};
    
    fprintf('\n========================================\n');
    fprintf('ANALYZING: %s\n', basis_title);
    fprintf('========================================\n');
    
    % Store results for this basis
    if strcmp(basis_type, 'molecular')
        separation_results = separation_results_mol;
    else
        separation_results = separation_results_ion;
    end
    
    for cf = 1:length(class_fields)
        class_field = class_fields{cf};
        class_title = class_titles{cf};
        
        % Extract unique classes
        classes = {};
        for s = 1:length(salt_names)
            class_val = classifications.(salt_names{s}).(class_field);
            if ~any(strcmp(classes, class_val)) && ~strcmp(class_val, 'Unknown')
                classes{end+1} = class_val;
            end
        end
        
        % Skip if only one class or no valid classes
        if length(classes) < 2
            fprintf('\n--- %s ---\n', class_title);
            fprintf('  Insufficient classes for comparison (only %d class)\n', length(classes));
            continue;
        end
        
        % Calculate separation metrics across all test RH values
        all_f_statistics = [];
        all_mean_separations = [];
        all_coefficient_variations = [];
        valid_rh_count = 0;
        
        for rh_idx = 1:length(test_RH_values)
            target_RH = test_RH_values(rh_idx);
            
            class_values = containers.Map();
            
            % Collect values for each class
            for c = 1:length(classes)
                class_val = classes{c};
                values = [];
                
                for s = 1:length(salt_names)
                    salt_name = salt_names{s};
                    if strcmp(classifications.(salt_name).(class_field), class_val)
                        data = all_salts_data.(salt_name);
                        
                        % Check if target RH is in valid range
                        if target_RH >= min(data.RH) && target_RH <= max(data.RH)
                            % Find closest RH value
                            [~, idx] = min(abs(data.RH - target_RH));
                            values = [values; data.(metric_field)(idx)];
                        end
                    end
                end
                
                if length(values) >= 2  % Need at least 2 samples per class
                    class_values(class_val) = values;
                end
            end
            
            % Calculate metrics if we have at least 2 classes with data
            if length(class_values) >= 2
                valid_rh_count = valid_rh_count + 1;
                
                % Extract all values and group labels
                all_vals = [];
                group_labels = {};
                for c = 1:length(classes)
                    class_val = classes{c};
                    if isKey(class_values, class_val)
                        vals = class_values(class_val);
                        all_vals = [all_vals; vals];
                        group_labels = [group_labels; repmat({class_val}, length(vals), 1)];
                    end
                end
                
                % Calculate F-statistic (between-group variance / within-group variance)
                group_means = [];
                group_stds = [];
                group_counts = [];
                overall_mean = mean(all_vals);
                
                for c = 1:length(classes)
                    class_val = classes{c};
                    if isKey(class_values, class_val)
                        vals = class_values(class_val);
                        group_means = [group_means; mean(vals)];
                        group_stds = [group_stds; std(vals)];
                        group_counts = [group_counts; length(vals)];
                    end
                end
                
                % Between-group sum of squares
                SS_between = sum(group_counts .* (group_means - overall_mean).^2);
                df_between = length(group_means) - 1;
                MS_between = SS_between / df_between;
                
                % Within-group sum of squares
                SS_within = sum((group_counts - 1) .* group_stds.^2);
                df_within = length(all_vals) - length(group_means);
                MS_within = SS_within / df_within;
                
                % F-statistic
                if MS_within > 0
                    f_stat = MS_between / MS_within;
                    all_f_statistics = [all_f_statistics; f_stat];
                end
                
                % Mean separation (distance between class means)
                if length(group_means) >= 2
                    mean_sep = max(group_means) - min(group_means);
                    all_mean_separations = [all_mean_separations; mean_sep];
                end
                
                % Coefficient of variation (std/mean) - lower is better separation
                cv = std(group_means) / abs(mean(group_means));
                if ~isnan(cv) && ~isinf(cv)
                    all_coefficient_variations = [all_coefficient_variations; cv];
                end
            end
        end
        
        % Store results
        if valid_rh_count > 0
            separation_results.(class_field).mean_f_stat = mean(all_f_statistics);
            separation_results.(class_field).mean_separation = mean(all_mean_separations);
            separation_results.(class_field).mean_cv = mean(all_coefficient_variations);
            separation_results.(class_field).num_classes = length(classes);
            separation_results.(class_field).valid_rh_count = valid_rh_count;
            separation_results.(class_field).title = class_title;
        end
    end
    
    % Store back to the appropriate structure
    if strcmp(basis_type, 'molecular')
        separation_results_mol = separation_results;
    else
        separation_results_ion = separation_results;
    end
end

%% ========================================================================
%% DETERMINE MOST EFFECTIVE CLASSIFICATION FOR EACH BASIS
%% ========================================================================

% Find best for molecular basis
fprintf('\n========================================\n');
fprintf('SEPARATION METRICS SUMMARY - Molecular Basis\n');
fprintf('========================================\n');

result_fields_mol = fieldnames(separation_results_mol);
best_classification_mol = '';
best_score_mol = -inf;

for rf = 1:length(result_fields_mol)
    field = result_fields_mol{rf};
    result = separation_results_mol.(field);
    
    fprintf('\n--- %s ---\n', result.title);
    fprintf('  Number of classes: %d\n', result.num_classes);
    fprintf('  Valid RH points: %d\n', result.valid_rh_count);
    fprintf('  Mean F-statistic: %.4f (higher = better separation)\n', result.mean_f_stat);
    fprintf('  Mean separation: %.4f (higher = better separation)\n', result.mean_separation);
    fprintf('  Mean CV: %.4f (lower = better separation)\n', result.mean_cv);
    
    % Combined score: F-statistic * mean_separation / CV
    combined_score = result.mean_f_stat * result.mean_separation / (result.mean_cv + 0.001);
    
    if combined_score > best_score_mol
        best_score_mol = combined_score;
        best_classification_mol = field;
    end
end

% Also check individual metrics for molecular
[~, idx_f_mol] = max(cellfun(@(x) separation_results_mol.(x).mean_f_stat, result_fields_mol));
[~, idx_sep_mol] = max(cellfun(@(x) separation_results_mol.(x).mean_separation, result_fields_mol));
[~, idx_cv_mol] = min(cellfun(@(x) separation_results_mol.(x).mean_cv, result_fields_mol));

fprintf('\n========================================\n');
fprintf('BEST CLASSIFICATIONS BY METRIC - Molecular Basis\n');
fprintf('========================================\n');
fprintf('Best by F-statistic: %s (F=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_f_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_f_mol}).mean_f_stat);
fprintf('Best by mean separation: %s (sep=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_sep_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_sep_mol}).mean_separation);
fprintf('Best by CV (lowest): %s (CV=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_cv_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_cv_mol}).mean_cv);
fprintf('Best by combined score: %s (score=%.4f)\n', ...
        separation_results_mol.(best_classification_mol).title, best_score_mol);

% Find best for ionic basis
fprintf('\n========================================\n');
fprintf('SEPARATION METRICS SUMMARY - Ionic Basis\n');
fprintf('========================================\n');

result_fields_ion = fieldnames(separation_results_ion);
best_classification_ion = '';
best_score_ion = -inf;

for rf = 1:length(result_fields_ion)
    field = result_fields_ion{rf};
    result = separation_results_ion.(field);
    
    fprintf('\n--- %s ---\n', result.title);
    fprintf('  Number of classes: %d\n', result.num_classes);
    fprintf('  Valid RH points: %d\n', result.valid_rh_count);
    fprintf('  Mean F-statistic: %.4f (higher = better separation)\n', result.mean_f_stat);
    fprintf('  Mean separation: %.4f (higher = better separation)\n', result.mean_separation);
    fprintf('  Mean CV: %.4f (lower = better separation)\n', result.mean_cv);
    
    % Combined score: F-statistic * mean_separation / CV
    combined_score = result.mean_f_stat * result.mean_separation / (result.mean_cv + 0.001);
    
    if combined_score > best_score_ion
        best_score_ion = combined_score;
        best_classification_ion = field;
    end
end

% Also check individual metrics for ionic
[~, idx_f_ion] = max(cellfun(@(x) separation_results_ion.(x).mean_f_stat, result_fields_ion));
[~, idx_sep_ion] = max(cellfun(@(x) separation_results_ion.(x).mean_separation, result_fields_ion));
[~, idx_cv_ion] = min(cellfun(@(x) separation_results_ion.(x).mean_cv, result_fields_ion));

fprintf('\n========================================\n');
fprintf('BEST CLASSIFICATIONS BY METRIC - Ionic Basis\n');
fprintf('========================================\n');
fprintf('Best by F-statistic: %s (F=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_f_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_f_ion}).mean_f_stat);
fprintf('Best by mean separation: %s (sep=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_sep_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_sep_ion}).mean_separation);
fprintf('Best by CV (lowest): %s (CV=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_cv_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_cv_ion}).mean_cv);
fprintf('Best by combined score: %s (score=%.4f)\n', ...
        separation_results_ion.(best_classification_ion).title, best_score_ion);

% Compare results
fprintf('\n========================================\n');
fprintf('COMPARISON: MOLECULAR vs IONIC BASIS\n');
fprintf('========================================\n');
fprintf('Molecular Basis - Best Classification: %s (score=%.4f)\n', ...
        separation_results_mol.(best_classification_mol).title, best_score_mol);
fprintf('Ionic Basis - Best Classification: %s (score=%.4f)\n', ...
        separation_results_ion.(best_classification_ion).title, best_score_ion);

% Use best classification for each basis in plotting
best_classification_mol_final = best_classification_mol;
best_classification_ion_final = best_classification_ion;

%% ========================================================================
%% PLOT ln(gamma) vs RH COLORED BY BEST CLASSIFICATION
%% ========================================================================

fprintf('\n========================================\n');
fprintf('GENERATING PLOTS\n');
fprintf('========================================\n');

% Plot for Molecular Basis using its best classification
fprintf('Using classification for Molecular Basis: %s\n', ...
        get_classification_title(best_classification_mol_final));

% Extract unique classes for molecular basis best classification
classes_mol = {};
for s = 1:length(salt_names)
    class_val = classifications.(salt_names{s}).(best_classification_mol_final);
    if ~any(strcmp(classes_mol, class_val)) && ~strcmp(class_val, 'Unknown')
        classes_mol{end+1} = class_val;
    end
end

% Define colors for molecular plot
num_classes_mol = length(classes_mol);
if num_classes_mol <= 2
    colors_map_mol = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410];
elseif num_classes_mol <= 4
    colors_map_mol = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; ...
                      0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560];
else
    colors_map_mol = hsv(num_classes_mol);
end

alpha = 0.6;

figure('Position', [300, 300, 1200, 800]);
hold on; grid on; box on;

for c = 1:length(classes_mol)
    class_val = classes_mol{c};
    if c <= size(colors_map_mol, 1)
        color = colors_map_mol(c, :);
    else
        color = [0.5, 0.5, 0.5];
    end
    count = 0;
    
    for s = 1:length(salt_names)
        salt_name = salt_names{s};
        if strcmp(classifications.(salt_name).(best_classification_mol_final), class_val)
            data = all_salts_data.(salt_name);
            h = plot(data.RH*100, data.ln_gamma_w_mol, 'LineWidth', 2, ...
                     'color', color, 'HandleVisibility', 'off');
            try
                h.Color(4) = alpha;
            catch
                h.Color = color * alpha + [1, 1, 1] * (1-alpha);
            end
            count = count + 1;
        end
    end
    
    plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', color, ...
               'DisplayName', sprintf('%s (n=%d)', class_val, count));
end

plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('ln(\\gamma_w) vs RH - Molecular Basis - Colored by %s', ...
      get_classification_title(best_classification_mol_final)), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 12)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

filename = sprintf('ln_gamma_vs_RH_molecular_best_classification_%s', best_classification_mol_final);
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
fprintf('Plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));

% Plot for Ionic Basis using its best classification
fprintf('Using classification for Ionic Basis: %s\n', ...
        get_classification_title(best_classification_ion_final));

% Extract unique classes for ionic basis best classification
classes_ion = {};
for s = 1:length(salt_names)
    class_val = classifications.(salt_names{s}).(best_classification_ion_final);
    if ~any(strcmp(classes_ion, class_val)) && ~strcmp(class_val, 'Unknown')
        classes_ion{end+1} = class_val;
    end
end

% Define colors for ionic plot
num_classes_ion = length(classes_ion);
if num_classes_ion <= 2
    colors_map_ion = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410];
elseif num_classes_ion <= 4
    colors_map_ion = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; ...
                      0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560];
else
    colors_map_ion = hsv(num_classes_ion);
end

figure('Position', [300, 300, 1200, 800]);
hold on; grid on; box on;

for c = 1:length(classes_ion)
    class_val = classes_ion{c};
    if c <= size(colors_map_ion, 1)
        color = colors_map_ion(c, :);
    else
        color = [0.5, 0.5, 0.5];
    end
    count = 0;
    
    for s = 1:length(salt_names)
        salt_name = salt_names{s};
        if strcmp(classifications.(salt_name).(best_classification_ion_final), class_val)
            data = all_salts_data.(salt_name);
            h = plot(data.RH*100, data.ln_gamma_w_ion, 'LineWidth', 2, ...
                     'color', color, 'HandleVisibility', 'off');
            try
                h.Color(4) = alpha;
            catch
                h.Color = color * alpha + [1, 1, 1] * (1-alpha);
            end
            count = count + 1;
        end
    end
    
    plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', color, ...
               'DisplayName', sprintf('%s (n=%d)', class_val, count));
end

plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('ln(\\gamma_w) vs RH - Ionic Basis - Colored by %s', ...
      get_classification_title(best_classification_ion_final)), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 12)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

filename = sprintf('ln_gamma_vs_RH_ionic_best_classification_%s', best_classification_ion_final);
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
fprintf('Plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));

%% ========================================================================
%% SAVE RESULTS TO FILE
%% ========================================================================

% Write results summary to text file
results_file = fullfile(fig_out_dir, 'classification_analysis_results.txt');
fid = fopen(results_file, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'CLASSIFICATION EFFECTIVENESS ANALYSIS\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Analysis Date: %s\n\n', datestr(now));

fprintf(fid, 'METRICS EXPLANATION:\n');
fprintf(fid, '  - F-statistic: Ratio of between-group to within-group variance\n');
fprintf(fid, '    Higher values indicate better separation between classes\n');
fprintf(fid, '  - Mean Separation: Average distance between class means\n');
fprintf(fid, '    Higher values indicate classes are further apart\n');
fprintf(fid, '  - Coefficient of Variation: Standard deviation / mean of class means\n');
fprintf(fid, '    Lower values indicate more consistent separation\n\n');

fprintf(fid, '========================================\n');
fprintf(fid, 'SEPARATION METRICS BY CLASSIFICATION - MOLECULAR BASIS\n');
fprintf(fid, '========================================\n\n');

for rf = 1:length(result_fields_mol)
    field = result_fields_mol{rf};
    result = separation_results_mol.(field);
    
    fprintf(fid, '--- %s ---\n', result.title);
    fprintf(fid, '  Number of classes: %d\n', result.num_classes);
    fprintf(fid, '  Valid RH points: %d\n', result.valid_rh_count);
    fprintf(fid, '  Mean F-statistic: %.4f\n', result.mean_f_stat);
    fprintf(fid, '  Mean separation: %.4f\n', result.mean_separation);
    fprintf(fid, '  Mean CV: %.4f\n\n', result.mean_cv);
end

fprintf(fid, '========================================\n');
fprintf(fid, 'BEST CLASSIFICATIONS - MOLECULAR BASIS\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Best by F-statistic: %s (F=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_f_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_f_mol}).mean_f_stat);
fprintf(fid, 'Best by mean separation: %s (sep=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_sep_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_sep_mol}).mean_separation);
fprintf(fid, 'Best by CV (lowest): %s (CV=%.4f)\n', ...
        separation_results_mol.(result_fields_mol{idx_cv_mol}).title, ...
        separation_results_mol.(result_fields_mol{idx_cv_mol}).mean_cv);
fprintf(fid, 'Best by combined score: %s (score=%.4f)\n\n', ...
        separation_results_mol.(best_classification_mol_final).title, best_score_mol);

fprintf(fid, '========================================\n');
fprintf(fid, 'SEPARATION METRICS BY CLASSIFICATION - IONIC BASIS\n');
fprintf(fid, '========================================\n\n');

for rf = 1:length(result_fields_ion)
    field = result_fields_ion{rf};
    result = separation_results_ion.(field);
    
    fprintf(fid, '--- %s ---\n', result.title);
    fprintf(fid, '  Number of classes: %d\n', result.num_classes);
    fprintf(fid, '  Valid RH points: %d\n', result.valid_rh_count);
    fprintf(fid, '  Mean F-statistic: %.4f\n', result.mean_f_stat);
    fprintf(fid, '  Mean separation: %.4f\n', result.mean_separation);
    fprintf(fid, '  Mean CV: %.4f\n\n', result.mean_cv);
end

fprintf(fid, '========================================\n');
fprintf(fid, 'BEST CLASSIFICATIONS - IONIC BASIS\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Best by F-statistic: %s (F=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_f_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_f_ion}).mean_f_stat);
fprintf(fid, 'Best by mean separation: %s (sep=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_sep_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_sep_ion}).mean_separation);
fprintf(fid, 'Best by CV (lowest): %s (CV=%.4f)\n', ...
        separation_results_ion.(result_fields_ion{idx_cv_ion}).title, ...
        separation_results_ion.(result_fields_ion{idx_cv_ion}).mean_cv);
fprintf(fid, 'Best by combined score: %s (score=%.4f)\n\n', ...
        separation_results_ion.(best_classification_ion_final).title, best_score_ion);

fprintf(fid, '========================================\n');
fprintf(fid, 'COMPARISON: MOLECULAR vs IONIC BASIS\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Molecular Basis - Best Classification: %s (score=%.4f)\n', ...
        separation_results_mol.(best_classification_mol_final).title, best_score_mol);
fprintf(fid, 'Ionic Basis - Best Classification: %s (score=%.4f)\n\n', ...
        separation_results_ion.(best_classification_ion_final).title, best_score_ion);

fprintf(fid, '========================================\n');
fprintf(fid, 'CONCLUSION\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'The most effective classification for separating salts based on\n');
fprintf(fid, 'ln(gamma_w) vs RH relationships:\n\n');
fprintf(fid, '  Molecular Basis: %s (score=%.4f)\n', ...
        separation_results_mol.(best_classification_mol_final).title, best_score_mol);
fprintf(fid, '  Ionic Basis: %s (score=%.4f)\n\n', ...
        separation_results_ion.(best_classification_ion_final).title, best_score_ion);
fprintf(fid, 'These classifications were selected based on a combined score that\n');
fprintf(fid, 'considers F-statistic, mean separation, and coefficient of variation.\n');

fclose(fid);
fprintf('\nResults summary saved to: %s\n', results_file);

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');

%% ========================================================================
%% HELPER FUNCTION
%% ========================================================================

function title_str = get_classification_title(class_field)
    if strcmp(class_field, 'thermal')
        title_str = 'Endothermic vs Exothermic';
    elseif strcmp(class_field, 'cation_class')
        title_str = 'Cation Class';
    elseif strcmp(class_field, 'anion_class')
        title_str = 'Anion Class';
    elseif strcmp(class_field, 'electrolyte_type')
        title_str = 'Electrolyte Type';
    elseif strcmp(class_field, 'anion_family')
        title_str = 'Anion Family';
    elseif strcmp(class_field, 'combined_ion_class')
        title_str = 'Combined Ion Classes';
    else
        title_str = 'Unknown';
    end
end
