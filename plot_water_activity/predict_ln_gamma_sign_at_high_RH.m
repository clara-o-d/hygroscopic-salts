close all 
clear
clc 

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'ln_gamma_sign_prediction');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015;

% Define all salts with their molecular weights, valid RH ranges, and ION COUNTS
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

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    n_cat = salt_data{s}{7};
    n_an  = salt_data{s}{8};
    nu = n_cat + n_an;
    
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    mf_salt = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water_mol = zeros(size(RH_vec));
    x_water_ion = zeros(size(RH_vec));
    
    success = true;
    for i = 1:length(RH_vec)
        try
            if func_args == 1
                mf_salt(i) = feval(func_name, RH_vec(i), T);
            else
                mf_salt(i) = feval(func_name, RH_vec(i));
            end
            mf_water(i) = 1 - mf_salt(i);
            n_w = mf_water(i) / MWw;  % Moles of water
            n_s = mf_salt(i) / MW;    % Moles of salt (molecular)
            
            % Molecular basis: x_w = n_w / (n_w + n_s)
            % (1 mole NaCl counts as 1 mole)
            x_water_mol(i) = n_w / (n_w + n_s);
            
            % Ionic basis: x_w = n_w / (n_w + nu * n_s)
            % where nu = n_cat + n_an (total ions per salt molecule)
            % For NaCl: nu = 2 (1 Na+ + 1 Cl-), so 1 mole NaCl counts as 2 moles
            x_water_ion(i) = n_w / (n_w + (nu * n_s));
        catch ME
            warning('Error processing %s at RH=%.4f: %s', salt_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        gamma_w_mol = RH_vec ./ x_water_mol;
        gamma_w_ion = RH_vec ./ x_water_ion;
        ln_gamma_w_mol = log(gamma_w_mol);
        ln_gamma_w_ion = log(gamma_w_ion);
        
        all_salts_data.(salt_name) = struct();
        all_salts_data.(salt_name).RH = RH_vec;
        all_salts_data.(salt_name).ln_gamma_w_mol = ln_gamma_w_mol;
        all_salts_data.(salt_name).ln_gamma_w_ion = ln_gamma_w_ion;
        all_salts_data.(salt_name).MW = MW;
        all_salts_data.(salt_name).nu = nu;
        all_salts_data.(salt_name).n_cat = n_cat;
        all_salts_data.(salt_name).n_an = n_an;
    end
end

salt_names = fieldnames(all_salts_data);
num_salts = length(salt_names);
fprintf('Successfully processed %d salts\n', num_salts);

%% ========================================================================
%% DETERMINE ln(gamma) SIGN AT HIGH RH
%% ========================================================================

% Define high RH threshold
high_RH_threshold = 0.85;  % Consider RH >= 0.85 as "high RH"
fprintf('\nAnalyzing ln(gamma) sign at RH >= %.2f\n', high_RH_threshold);

% Use ionic basis only
ln_gamma_field = 'ln_gamma_w_ion';
basis_title = 'Ionic Basis';

fprintf('\n--- %s ---\n', basis_title);

% For each salt, determine if ln(gamma) > 0 at high RH
salt_labels = struct();
salt_high_RH_values = struct();

for s = 1:num_salts
    salt_name = salt_names{s};
    data = all_salts_data.(salt_name);
    
    % Find indices where RH >= threshold
    high_RH_idx = data.RH >= high_RH_threshold;
    
    if any(high_RH_idx)
        ln_gamma_high_RH = data.(ln_gamma_field)(high_RH_idx);
        % Use mean value at high RH to determine label
        mean_ln_gamma_high = mean(ln_gamma_high_RH);
        salt_labels.(salt_name) = double(mean_ln_gamma_high > 0);
        salt_high_RH_values.(salt_name) = mean_ln_gamma_high;
    else
        % If salt doesn't reach high RH, use maximum available RH
        [max_RH, max_idx] = max(data.RH);
        salt_labels.(salt_name) = double(data.(ln_gamma_field)(max_idx) > 0);
        salt_high_RH_values.(salt_name) = data.(ln_gamma_field)(max_idx);
    end
end

% Count positive and negative
num_positive = sum(cellfun(@(x) salt_labels.(x), salt_names));
num_negative = num_salts - num_positive;
fprintf('Salts with ln(gamma) > 0 at high RH: %d\n', num_positive);
fprintf('Salts with ln(gamma) <= 0 at high RH: %d\n', num_negative);

%% ========================================================================
%% CLASSIFICATION DEFINITIONS
%% ========================================================================

endothermic_salts = {'NaCl', 'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', ...
                     'LiNO3', 'KNO3', 'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', ...
                     'RbCl', 'CsBr', 'CsI', 'Na2SO4', 'K2SO4', 'NH42SO4', ...
                     'MgSO4', 'MnSO4', 'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4', ...
                     'BaNO3', 'BaCl2'};

cation_kosmotropes = {'Li', 'Na', 'Mg', 'Ca', 'Sr', 'Ba', 'Zn', 'Cu', 'Ni', 'Mn', 'Ag', 'H'};
cation_chaotropes = {'K', 'Rb', 'Cs', 'NH4'};
anion_kosmotropes = {'SO4', 'OH'};
anion_chaotropes = {'I', 'Br', 'NO3', 'ClO4', 'ClO3'};
anion_intermediate = {'Cl'};

% Helper functions
function cation = extract_cation(salt_name)
    if startsWith(salt_name, 'NH4'), cation = 'NH4';
    elseif startsWith(salt_name, 'Li'), cation = 'Li';
    elseif startsWith(salt_name, 'Na'), cation = 'Na';
    elseif startsWith(salt_name, 'K'), cation = 'K';
    elseif startsWith(salt_name, 'Rb'), cation = 'Rb';
    elseif startsWith(salt_name, 'Cs'), cation = 'Cs';
    elseif startsWith(salt_name, 'Mg'), cation = 'Mg';
    elseif startsWith(salt_name, 'Ca'), cation = 'Ca';
    elseif startsWith(salt_name, 'Sr'), cation = 'Sr';
    elseif startsWith(salt_name, 'Ba'), cation = 'Ba';
    elseif startsWith(salt_name, 'Zn'), cation = 'Zn';
    elseif startsWith(salt_name, 'Cu'), cation = 'Cu';
    elseif startsWith(salt_name, 'Ni'), cation = 'Ni';
    elseif startsWith(salt_name, 'Mn'), cation = 'Mn';
    elseif startsWith(salt_name, 'Ag'), cation = 'Ag';
    elseif startsWith(salt_name, 'H'), cation = 'H';
    else, cation = 'Unknown';
    end
end

function anion = extract_anion(salt_name)
    if contains(salt_name, 'SO4'), anion = 'SO4';
    elseif contains(salt_name, 'NO3'), anion = 'NO3';
    elseif contains(salt_name, 'ClO4'), anion = 'ClO4';
    elseif contains(salt_name, 'ClO3'), anion = 'ClO3';
    elseif contains(salt_name, 'OH'), anion = 'OH';
    elseif contains(salt_name, 'Cl'), anion = 'Cl';
    elseif contains(salt_name, 'Br'), anion = 'Br';
    elseif contains(salt_name, 'I'), anion = 'I';
    else, anion = 'Unknown';
    end
end

% Apply classifications
classifications = struct();
for s = 1:num_salts
    salt_name = salt_names{s};
    
    classifications.(salt_name).thermal = 'Exothermic';
    if any(strcmp(salt_name, endothermic_salts))
        classifications.(salt_name).thermal = 'Endothermic';
    end
    
    cation = extract_cation(salt_name);
    anion = extract_anion(salt_name);
    classifications.(salt_name).cation = cation;
    classifications.(salt_name).anion = anion;
    
    classifications.(salt_name).cation_class = 'Unknown';
    if any(strcmp(cation, cation_kosmotropes))
        classifications.(salt_name).cation_class = 'Kosmotrope';
    elseif any(strcmp(cation, cation_chaotropes))
        classifications.(salt_name).cation_class = 'Chaotrope';
    end
    
    classifications.(salt_name).anion_class = 'Unknown';
    if any(strcmp(anion, anion_kosmotropes))
        classifications.(salt_name).anion_class = 'Kosmotrope';
    elseif any(strcmp(anion, anion_chaotropes))
        classifications.(salt_name).anion_class = 'Chaotrope';
    elseif any(strcmp(anion, anion_intermediate))
        classifications.(salt_name).anion_class = 'Intermediate';
    end
    
    n_cat = all_salts_data.(salt_name).n_cat;
    n_an = all_salts_data.(salt_name).n_an;
    classifications.(salt_name).electrolyte_type = sprintf('%d:%d', n_cat, n_an);
    
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
    
    cat_class = classifications.(salt_name).cation_class;
    an_class = classifications.(salt_name).anion_class;
    classifications.(salt_name).combined_ion_class = sprintf('%s-Cat/%s-An', cat_class, an_class);
end

%% ========================================================================
%% EVALUATE CLASSIFICATION SCHEMES FOR PREDICTING ln(gamma) SIGN
%% ========================================================================

fprintf('\n========================================\n');
fprintf('EVALUATING CLASSIFICATION SCHEMES\n');
fprintf('========================================\n');

class_fields = {'thermal', 'cation_class', 'anion_class', 'electrolyte_type', ...
               'anion_family', 'combined_ion_class'};
class_titles = {'Endothermic vs Exothermic', 'Cation Class', 'Anion Class', ...
               'Electrolyte Type', 'Anion Family', 'Combined Ion Classes'};

% Evaluate classification schemes for ionic basis
classification_scores = struct();

fprintf('\n========================================\n');
fprintf('EVALUATING CLASSIFICATION SCHEMES - Ionic Basis\n');
fprintf('========================================\n');

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
        
        if length(classes) < 2
            continue;
        end
        
        % Calculate confusion matrix for each class
        true_positives = 0;
        true_negatives = 0;
        false_positives = 0;
        false_negatives = 0;
        
        % For each class, determine if it predicts positive or negative
        class_predictions = containers.Map();
        for c = 1:length(classes)
            class_val = classes{c};
            positive_count = 0;
            total_count = 0;
            
        for s = 1:length(salt_names)
            salt_name = salt_names{s};
            if strcmp(classifications.(salt_name).(class_field), class_val)
                total_count = total_count + 1;
                if salt_labels.(salt_name) == 1
                    positive_count = positive_count + 1;
                end
            end
        end
        
        % Class predicts positive if majority are positive
        class_predictions(class_val) = double(positive_count > total_count / 2);
    end
    
    % Calculate accuracy
    correct = 0;
    total = 0;
    for s = 1:length(salt_names)
        salt_name = salt_names{s};
        class_val = classifications.(salt_name).(class_field);
        if isKey(class_predictions, class_val)
            predicted = class_predictions(class_val);
            actual = salt_labels.(salt_name);
            if predicted == actual
                correct = correct + 1;
            end
            total = total + 1;
            
            % Update confusion matrix
            if predicted == 1 && actual == 1
                true_positives = true_positives + 1;
            elseif predicted == 0 && actual == 0
                true_negatives = true_negatives + 1;
            elseif predicted == 1 && actual == 0
                false_positives = false_positives + 1;
            elseif predicted == 0 && actual == 1
                false_negatives = false_negatives + 1;
            end
        end
    end
    
    accuracy = correct / total;
    precision = true_positives / (true_positives + false_positives + eps);
    recall = true_positives / (true_positives + false_negatives + eps);
    f1_score = 2 * (precision * recall) / (precision + recall + eps);
    
    classification_scores.(class_field).accuracy = accuracy;
    classification_scores.(class_field).precision = precision;
    classification_scores.(class_field).recall = recall;
    classification_scores.(class_field).f1_score = f1_score;
    classification_scores.(class_field).title = class_title;
    classification_scores.(class_field).num_classes = length(classes);
    
    fprintf('\n--- %s ---\n', class_title);
    fprintf('  Accuracy: %.4f\n', accuracy);
    fprintf('  Precision: %.4f\n', precision);
    fprintf('  Recall: %.4f\n', recall);
    fprintf('  F1 Score: %.4f\n', f1_score);
end

% Find best classification
result_fields = fieldnames(classification_scores);
[~, best_idx] = max(cellfun(@(x) classification_scores.(x).f1_score, result_fields));
best_classification_field = result_fields{best_idx};

fprintf('\n========================================\n');
fprintf('BEST CLASSIFICATION SCHEME - Ionic Basis\n');
fprintf('========================================\n');
fprintf('Best by F1 Score: %s\n', classification_scores.(best_classification_field).title);
fprintf('  Accuracy: %.4f\n', classification_scores.(best_classification_field).accuracy);
fprintf('  F1 Score: %.4f\n', classification_scores.(best_classification_field).f1_score);

%% ========================================================================
%% PLOT: SPLIT BY BEST CLASSIFICATION
%% ========================================================================

fprintf('\n========================================\n');
fprintf('GENERATING PLOTS\n');
fprintf('========================================\n');

fprintf('Creating plot for Ionic Basis using classification: %s\n', ...
        classification_scores.(best_classification_field).title);

% Extract classes for best classification
classes_best = {};
for s = 1:length(salt_names)
    class_val = classifications.(salt_names{s}).(best_classification_field);
    if ~any(strcmp(classes_best, class_val)) && ~strcmp(class_val, 'Unknown')
        classes_best{end+1} = class_val;
    end
end

% Determine prediction for each class
class_predictions_best = containers.Map();
for c = 1:length(classes_best)
    class_val = classes_best{c};
    positive_count = 0;
    total_count = 0;
    
    for s = 1:length(salt_names)
        salt_name = salt_names{s};
        if strcmp(classifications.(salt_name).(best_classification_field), class_val)
            total_count = total_count + 1;
            if salt_labels.(salt_name) == 1
                positive_count = positive_count + 1;
            end
        end
    end
    
    class_predictions_best(class_val) = double(positive_count > total_count / 2);
end

% Create plot
figure('Position', [300, 300, 1400, 800]);
hold on; grid on; box on;

colors_map = hsv(length(classes_best));
alpha = 0.6;

for c = 1:length(classes_best)
    class_val = classes_best{c};
    color = colors_map(c, :);
    count = 0;
    
    for s = 1:length(salt_names)
        salt_name = salt_names{s};
        if strcmp(classifications.(salt_name).(best_classification_field), class_val)
            data = all_salts_data.(salt_name);
            predicted_positive = class_predictions_best(class_val);
            actual_positive = salt_labels.(salt_name);
            
            % Use different line style based on correctness
            if predicted_positive == actual_positive
                linestyle = '-';
            else
                linestyle = '--';  % Dashed for misclassified
            end
            
            h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', 2, ...
                     'color', color, 'LineStyle', linestyle, 'HandleVisibility', 'off');
            try
                h.Color(4) = alpha;
            catch
                h.Color = color * alpha + [1, 1, 1] * (1-alpha);
            end
            count = count + 1;
        end
    end
    
    % Add legend entry
    label = sprintf('%s (n=%d, pred=%s)', class_val, count, ...
                   char('+' * class_predictions_best(class_val) + '-' * (1-class_predictions_best(class_val))));
    plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', color, 'DisplayName', label);
end

plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'ln(\gamma_w) = 0')

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('ln(\\gamma_w) vs RH - Ionic Basis - Split by %s (Best Classification)', ...
      classification_scores.(best_classification_field).title), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 12)
xlim([0 100])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

filename = sprintf('ln_gamma_vs_RH_ionic_split_by_best_classification_%s', best_classification_field);
print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
fprintf('Plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));

%% ========================================================================
%% MACHINE LEARNING APPROACH
%% ========================================================================

fprintf('\n========================================\n');
fprintf('MACHINE LEARNING PREDICTION\n');
fprintf('========================================\n');

% Prepare feature matrix
% Features: thermal (binary), cation_class (categorical), anion_class (categorical),
%           electrolyte_type (categorical), anion_family (categorical),
%           MW, n_cat, n_an, nu

% Convert categorical to numeric
thermal_numeric = zeros(num_salts, 1);
cation_class_numeric = zeros(num_salts, 1);
anion_class_numeric = zeros(num_salts, 1);
electrolyte_type_numeric = zeros(num_salts, 1);
anion_family_numeric = zeros(num_salts, 1);

for s = 1:num_salts
    salt_name = salt_names{s};
    thermal_numeric(s) = strcmp(classifications.(salt_name).thermal, 'Endothermic');
    
    % Cation class: 1=Kosmotrope, 2=Chaotrope, 0=Unknown
    if strcmp(classifications.(salt_name).cation_class, 'Kosmotrope')
        cation_class_numeric(s) = 1;
    elseif strcmp(classifications.(salt_name).cation_class, 'Chaotrope')
        cation_class_numeric(s) = 2;
    end
    
    % Anion class: 1=Kosmotrope, 2=Chaotrope, 3=Intermediate, 0=Unknown
    if strcmp(classifications.(salt_name).anion_class, 'Kosmotrope')
        anion_class_numeric(s) = 1;
    elseif strcmp(classifications.(salt_name).anion_class, 'Chaotrope')
        anion_class_numeric(s) = 2;
    elseif strcmp(classifications.(salt_name).anion_class, 'Intermediate')
        anion_class_numeric(s) = 3;
    end
    
    % Electrolyte type: encode as number
    if strcmp(classifications.(salt_name).electrolyte_type, '1:1')
        electrolyte_type_numeric(s) = 1;
    elseif strcmp(classifications.(salt_name).electrolyte_type, '1:2')
        electrolyte_type_numeric(s) = 2;
    elseif strcmp(classifications.(salt_name).electrolyte_type, '2:1')
        electrolyte_type_numeric(s) = 3;
    end
    
    % Anion family: encode as number
    an_fam = classifications.(salt_name).anion_family;
    if strcmp(an_fam, 'Halide'), anion_family_numeric(s) = 1;
    elseif strcmp(an_fam, 'Sulfate'), anion_family_numeric(s) = 2;
    elseif strcmp(an_fam, 'Nitrate'), anion_family_numeric(s) = 3;
    elseif strcmp(an_fam, 'Perchlorate/Chlorate'), anion_family_numeric(s) = 4;
    elseif strcmp(an_fam, 'Hydroxide'), anion_family_numeric(s) = 5;
    end
end

% Extract numeric features and ensure they are column vectors
MW_values = cell2mat(cellfun(@(x) all_salts_data.(x).MW, salt_names, 'UniformOutput', false));
n_cat_values = cell2mat(cellfun(@(x) all_salts_data.(x).n_cat, salt_names, 'UniformOutput', false));
n_an_values = cell2mat(cellfun(@(x) all_salts_data.(x).n_an, salt_names, 'UniformOutput', false));
nu_values = cell2mat(cellfun(@(x) all_salts_data.(x).nu, salt_names, 'UniformOutput', false));

% Ensure all are column vectors
MW_values = MW_values(:);
n_cat_values = n_cat_values(:);
n_an_values = n_an_values(:);
nu_values = nu_values(:);

% Create feature matrix (all columns must have same number of rows)
X = [thermal_numeric, cation_class_numeric, anion_class_numeric, ...
     electrolyte_type_numeric, anion_family_numeric, ...
     MW_values, n_cat_values, n_an_values, nu_values];

% Run ML for ionic basis only
ml_results = struct();

% Split into train and test sets (70/30 split)
rng(42);  % For reproducibility
train_ratio = 0.7;
num_train = round(num_salts * train_ratio);
indices = randperm(num_salts);
train_idx = indices(1:num_train);
test_idx = indices(num_train+1:end);

fprintf('\n========================================\n');
fprintf('MACHINE LEARNING PREDICTION - Ionic Basis\n');
fprintf('========================================\n');

% Create label vector (ensure it's a column vector)
y = cell2mat(cellfun(@(x) salt_labels.(x), salt_names, 'UniformOutput', false));
y = y(:);
    
    X_train = X(train_idx, :);
    y_train = y(train_idx);
    X_test = X(test_idx, :);
    y_test = y(test_idx);
    
    fprintf('Training set: %d samples\n', length(train_idx));
    fprintf('Test set: %d samples\n', length(test_idx));

% Try multiple ML algorithms
% Note: MATLAB's Statistics and Machine Learning Toolbox is required for some functions

% 1. Decision Tree
try
    tree = fitctree(X_train, y_train, 'PredictorNames', ...
        {'Thermal', 'CationClass', 'AnionClass', 'ElectrolyteType', ...
         'AnionFamily', 'MW', 'n_cat', 'n_an', 'nu'});
    y_pred_tree_train = predict(tree, X_train);
    y_pred_tree_test = predict(tree, X_test);
    
    acc_tree_train = sum(y_pred_tree_train == y_train) / length(y_train);
    acc_tree_test = sum(y_pred_tree_test == y_test) / length(y_test);
    
    fprintf('\nDecision Tree:\n');
    fprintf('  Training Accuracy: %.4f\n', acc_tree_train);
    fprintf('  Test Accuracy: %.4f\n', acc_tree_test);
    
    ml_results.tree = struct();
    ml_results.tree.accuracy_train = acc_tree_train;
    ml_results.tree.accuracy_test = acc_tree_test;
    ml_results.tree.model = tree;
    ml_results.tree.y_pred_train = y_pred_tree_train;
    ml_results.tree.y_pred_test = y_pred_tree_test;
catch ME
    fprintf('Decision Tree not available: %s\n', ME.message);
    ml_results.tree = [];
end

% 2. Random Forest (using TreeBagger)
try
    num_trees = 100;
    rf = TreeBagger(num_trees, X_train, y_train, 'Method', 'classification', ...
                    'OOBPrediction', 'on');
    y_pred_rf_train = str2double(predict(rf, X_train));
    y_pred_rf_test = str2double(predict(rf, X_test));
    
    acc_rf_train = sum(y_pred_rf_train == y_train) / length(y_train);
    acc_rf_test = sum(y_pred_rf_test == y_test) / length(y_test);
    
    fprintf('\nRandom Forest:\n');
    fprintf('  Training Accuracy: %.4f\n', acc_rf_train);
    fprintf('  Test Accuracy: %.4f\n', acc_rf_test);
    
    ml_results.rf = struct();
    ml_results.rf.accuracy_train = acc_rf_train;
    ml_results.rf.accuracy_test = acc_rf_test;
    ml_results.rf.model = rf;
    ml_results.rf.y_pred_train = y_pred_rf_train;
    ml_results.rf.y_pred_test = y_pred_rf_test;
catch ME
    fprintf('Random Forest not available: %s\n', ME.message);
    ml_results.rf = [];
end

% 3. Logistic Regression
try
    % Use glmfit for logistic regression
    [b, dev, stats] = glmfit(X_train, y_train, 'binomial', 'link', 'logit');
    y_pred_lr_train = glmval(b, X_train, 'logit') > 0.5;
    y_pred_lr_test = glmval(b, X_test, 'logit') > 0.5;
    
    acc_lr_train = sum(y_pred_lr_train == y_train) / length(y_train);
    acc_lr_test = sum(y_pred_lr_test == y_test) / length(y_test);
    
    fprintf('\nLogistic Regression:\n');
    fprintf('  Training Accuracy: %.4f\n', acc_lr_train);
    fprintf('  Test Accuracy: %.4f\n', acc_lr_test);
    
    ml_results.lr = struct();
    ml_results.lr.accuracy_train = acc_lr_train;
    ml_results.lr.accuracy_test = acc_lr_test;
    ml_results.lr.coefficients = b;
    ml_results.lr.y_pred_train = y_pred_lr_train;
    ml_results.lr.y_pred_test = y_pred_lr_test;
catch ME
    fprintf('Logistic Regression not available: %s\n', ME.message);
    ml_results.lr = [];
end

% Find best ML model
ml_fields = fieldnames(ml_results);
best_ml_model = '';
best_ml_test_acc = 0;

for mf = 1:length(ml_fields)
    field = ml_fields{mf};
    if ~isempty(ml_results.(field)) && isfield(ml_results.(field), 'accuracy_test')
        if ml_results.(field).accuracy_test > best_ml_test_acc
            best_ml_test_acc = ml_results.(field).accuracy_test;
            best_ml_model = field;
        end
    end
end

if ~isempty(best_ml_model)
    fprintf('\nBest ML Model for Ionic Basis: %s (Test Accuracy: %.4f)\n', ...
            best_ml_model, best_ml_test_acc);
    
    %% ========================================================================
    %% CONFUSION MATRIX ANALYSIS
    %% ========================================================================
    
    fprintf('\n========================================\n');
    fprintf('CONFUSION MATRIX ANALYSIS\n');
    fprintf('========================================\n');
    
    % Get predictions and actual labels
    y_pred_train = ml_results.(best_ml_model).y_pred_train;
    y_pred_test = ml_results.(best_ml_model).y_pred_test;
    
    % Calculate confusion matrix for training set
    % TP: predicted=1, actual=1
    % TN: predicted=0, actual=0
    % FP: predicted=1, actual=0
    % FN: predicted=0, actual=1
    TP_train = sum((y_pred_train == 1) & (y_train == 1));
    TN_train = sum((y_pred_train == 0) & (y_train == 0));
    FP_train = sum((y_pred_train == 1) & (y_train == 0));
    FN_train = sum((y_pred_train == 0) & (y_train == 1));
    
    % Calculate confusion matrix for test set
    TP_test = sum((y_pred_test == 1) & (y_test == 1));
    TN_test = sum((y_pred_test == 0) & (y_test == 0));
    FP_test = sum((y_pred_test == 1) & (y_test == 0));
    FN_test = sum((y_pred_test == 0) & (y_test == 1));
    
    % Store confusion matrices
    confusion_train = [TP_train, FP_train; FN_train, TN_train];
    confusion_test = [TP_test, FP_test; FN_test, TN_test];
    
    ml_results.(best_ml_model).confusion_train = confusion_train;
    ml_results.(best_ml_model).confusion_test = confusion_test;
    ml_results.(best_ml_model).TP_train = TP_train;
    ml_results.(best_ml_model).TN_train = TN_train;
    ml_results.(best_ml_model).FP_train = FP_train;
    ml_results.(best_ml_model).FN_train = FN_train;
    ml_results.(best_ml_model).TP_test = TP_test;
    ml_results.(best_ml_model).TN_test = TN_test;
    ml_results.(best_ml_model).FP_test = FP_test;
    ml_results.(best_ml_model).FN_test = FN_test;
    
    % Print confusion matrices
    fprintf('\nTraining Set Confusion Matrix:\n');
    fprintf('                              Predicted\n');
    fprintf('                    ln(γ_w)>0 (γ_w>1)  ln(γ_w)≤0 (γ_w≤1)\n');
    fprintf('Actual ln(γ_w)>0 (γ_w>1)      %3d              %3d\n', TP_train, FN_train);
    fprintf('Actual ln(γ_w)≤0 (γ_w≤1)     %3d              %3d\n', FP_train, TN_train);
    fprintf('\n  True Positives (TP):  %d (correctly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', TP_train);
    fprintf('  True Negatives (TN):  %d (correctly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', TN_train);
    fprintf('  False Positives (FP): %d (incorrectly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', FP_train);
    fprintf('  False Negatives (FN): %d (incorrectly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', FN_train);
    
    fprintf('\nTest Set Confusion Matrix:\n');
    fprintf('                              Predicted\n');
    fprintf('                    ln(γ_w)>0 (γ_w>1)  ln(γ_w)≤0 (γ_w≤1)\n');
    fprintf('Actual ln(γ_w)>0 (γ_w>1)      %3d              %3d\n', TP_test, FN_test);
    fprintf('Actual ln(γ_w)≤0 (γ_w≤1)     %3d              %3d\n', FP_test, TN_test);
    fprintf('\n  True Positives (TP):  %d (correctly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', TP_test);
    fprintf('  True Negatives (TN):  %d (correctly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', TN_test);
    fprintf('  False Positives (FP): %d (incorrectly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', FP_test);
    fprintf('  False Negatives (FN): %d (incorrectly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', FN_test);
    
    % Calculate additional metrics
    total_train = TP_train + TN_train + FP_train + FN_train;
    total_test = TP_test + TN_test + FP_test + FN_test;
    
    % Precision = TP / (TP + FP)
    precision_train = TP_train / (TP_train + FP_train + eps);
    precision_test = TP_test / (TP_test + FP_test + eps);
    
    % Recall (Sensitivity) = TP / (TP + FN)
    recall_train = TP_train / (TP_train + FN_train + eps);
    recall_test = TP_test / (TP_test + FN_test + eps);
    
    % Specificity = TN / (TN + FP)
    specificity_train = TN_train / (TN_train + FP_train + eps);
    specificity_test = TN_test / (TN_test + FP_test + eps);
    
    % F1 Score = 2 * (Precision * Recall) / (Precision + Recall)
    f1_train = 2 * (precision_train * recall_train) / (precision_train + recall_train + eps);
    f1_test = 2 * (precision_test * recall_test) / (precision_test + recall_test + eps);
    
    fprintf('\nTraining Set Metrics:\n');
    fprintf('  Accuracy:  %.4f (%.1f%%)\n', (TP_train + TN_train) / total_train, ...
            (TP_train + TN_train) / total_train * 100);
    fprintf('  Precision: %.4f\n', precision_train);
    fprintf('  Recall:    %.4f\n', recall_train);
    fprintf('  Specificity: %.4f\n', specificity_train);
    fprintf('  F1 Score:  %.4f\n', f1_train);
    
    fprintf('\nTest Set Metrics:\n');
    fprintf('  Accuracy:  %.4f (%.1f%%)\n', (TP_test + TN_test) / total_test, ...
            (TP_test + TN_test) / total_test * 100);
    fprintf('  Precision: %.4f\n', precision_test);
    fprintf('  Recall:    %.4f\n', recall_test);
    fprintf('  Specificity: %.4f\n', specificity_test);
    fprintf('  F1 Score:  %.4f\n', f1_test);
    
    % Create confusion matrix visualization
    figure('Position', [100, 100, 1200, 500]);
    
    % Training set confusion matrix
    subplot(1, 2, 1);
    axis([0.5 2.5 0.5 2.5]);
    axis square;
    hold on;
    % Draw grid lines
    plot([1.5 1.5], [0.5 2.5], 'k-', 'LineWidth', 1);
    plot([0.5 2.5], [1.5 1.5], 'k-', 'LineWidth', 1);
    % Add text for confusion matrix values
    % With YDir reverse: y=1 is top (Actual Positive), y=2 is bottom (Actual Negative)
    % x=1 is left (Predicted Positive), x=2 is right (Predicted Negative)
    text(1, 1, num2str(TP_train), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(2, 1, num2str(FN_train), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(1, 2, num2str(FP_train), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(2, 2, num2str(TN_train), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    % Add labels for TP, TN, FP, FN
    text(1, 0.7, 'TP', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(2, 0.7, 'FN', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(1, 2.3, 'FP', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(2, 2.3, 'TN', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'γ_w>1', 'γ_w≤1'}, ...
        'YTick', [1, 2], 'YTickLabel', {'γ_w>1', 'γ_w≤1'}, ...
        'TickLength', [0 0], 'YDir', 'reverse');
    xlabel('Predicted', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Actual', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Training Set Confusion Matrix\n(Accuracy: %.2f%%)', ...
          (TP_train + TN_train) / total_train * 100), ...
          'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12);
    hold off;
    
    % Test set confusion matrix
    subplot(1, 2, 2);
    axis([0.5 2.5 0.5 2.5]);
    axis square;
    hold on;
    % Draw grid lines
    plot([1.5 1.5], [0.5 2.5], 'k-', 'LineWidth', 1);
    plot([0.5 2.5], [1.5 1.5], 'k-', 'LineWidth', 1);
    % Add text for confusion matrix values
    % With YDir reverse: y=1 is top (Actual Positive), y=2 is bottom (Actual Negative)
    % x=1 is left (Predicted Positive), x=2 is right (Predicted Negative)
    text(1, 1, num2str(TP_test), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(2, 1, num2str(FN_test), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(1, 2, num2str(FP_test), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    text(2, 2, num2str(TN_test), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 20, 'FontWeight', 'bold');
    % Add labels for TP, TN, FP, FN
    text(1, 0.7, 'TP', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(2, 0.7, 'FN', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(1, 2.3, 'FP', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    text(2, 2.3, 'TN', 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 10);
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'γ_w>1', 'γ_w≤1'}, ...
        'YTick', [1, 2], 'YTickLabel', {'γ_w>1', 'γ_w≤1'}, ...
        'TickLength', [0 0], 'YDir', 'reverse');
    xlabel('Predicted', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Actual', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Test Set Confusion Matrix\n(Accuracy: %.2f%%)', ...
          (TP_test + TN_test) / total_test * 100), ...
          'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12);
    hold off;
    
    set(gcf, 'color', 'w');
    
    filename = sprintf('confusion_matrix_%s_ionic', best_ml_model);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    fprintf('\nConfusion matrix plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));
    
    %% ========================================================================
    %% FEATURE IMPORTANCE ANALYSIS - GAIN METRIC
    %% ========================================================================
    
    fprintf('\n========================================\n');
    fprintf('FEATURE IMPORTANCE ANALYSIS - GAIN METRIC\n');
    fprintf('========================================\n');
    
    feature_names = {'Thermal', 'CationClass', 'AnionClass', 'ElectrolyteType', ...
                     'AnionFamily', 'MW', 'n_cat', 'n_an', 'nu'};
    
    % Calculate feature gain: how much each feature contributes to model performance
    fprintf('\nCalculating feature gain (performance contribution)...\n');
    fprintf('This measures how much model performance drops when each feature is removed.\n\n');
    
    % Get baseline performance
    baseline_acc = ml_results.(best_ml_model).accuracy_test;
    fprintf('Baseline test accuracy (all features): %.4f\n\n', baseline_acc);
    
    % Calculate gain for each feature by removing it and retraining
    feature_gain = zeros(length(feature_names), 1);
    
    for feat_idx = 1:length(feature_names)
        fprintf('  Calculating gain for %s (%d/%d)...', feature_names{feat_idx}, feat_idx, length(feature_names));
        
        % Create feature matrix without this feature
        X_train_reduced = X_train;
        X_test_reduced = X_test;
        X_train_reduced(:, feat_idx) = [];
        X_test_reduced(:, feat_idx) = [];
        
        % Retrain model without this feature
        try
            if strcmp(best_ml_model, 'tree')
                tree_reduced = fitctree(X_train_reduced, y_train, 'PredictorNames', ...
                    feature_names([1:feat_idx-1, feat_idx+1:end]));
                y_pred_reduced = predict(tree_reduced, X_test_reduced);
                acc_reduced = sum(y_pred_reduced == y_test) / length(y_test);
            elseif strcmp(best_ml_model, 'rf')
                rf_reduced = TreeBagger(100, X_train_reduced, y_train, 'Method', 'classification', ...
                    'OOBPrediction', 'off');
                y_pred_reduced = str2double(predict(rf_reduced, X_test_reduced));
                acc_reduced = sum(y_pred_reduced == y_test) / length(y_test);
            elseif strcmp(best_ml_model, 'lr')
                [b_reduced, ~, ~] = glmfit(X_train_reduced, y_train, 'binomial', 'link', 'logit');
                y_pred_reduced = glmval(b_reduced, X_test_reduced, 'logit') > 0.5;
                acc_reduced = sum(y_pred_reduced == y_test) / length(y_test);
            else
                acc_reduced = baseline_acc;  % Fallback
            end
            
            % Gain is the drop in accuracy (how much we lose by removing the feature)
            % Higher drop = higher gain (more important feature)
            gain = baseline_acc - acc_reduced;
            feature_gain(feat_idx) = gain;
            
            fprintf(' Gain: %.4f (acc without feature: %.4f)\n', gain, acc_reduced);
        catch ME
            fprintf(' Error: %s\n', ME.message);
            feature_gain(feat_idx) = 0;  % If we can't calculate, assume no gain
        end
    end
    
    % Normalize gains (make them sum to 1, but keep negative gains as negative)
    % Shift to make all positive for normalization
    min_gain = min(feature_gain);
    if min_gain < 0
        feature_gain_shifted = feature_gain - min_gain;
    else
        feature_gain_shifted = feature_gain;
    end
    
    if sum(feature_gain_shifted) > 0
        feature_importance = feature_gain_shifted / sum(feature_gain_shifted);
    else
        feature_importance = ones(length(feature_names), 1) / length(feature_names);
    end
    
    fprintf('\n========================================\n');
    fprintf('FEATURE GAIN RESULTS\n');
    fprintf('========================================\n');
    fprintf('Gain = Baseline Accuracy - Accuracy Without Feature\n');
    fprintf('Higher gain = more important feature\n\n');
    
    % Sort features by gain
    [sorted_gain, sort_idx] = sort(feature_gain, 'descend');
    sorted_features = feature_names(sort_idx);
    
    % Print feature gain
    for i = 1:length(feature_names)
        fprintf('  %d. %s: Gain = %.4f (%.2f%% of baseline)\n', i, sorted_features{i}, ...
                sorted_gain(i), (sorted_gain(i) / baseline_acc) * 100);
    end
    
    % Use gain-based importance (already calculated above)
    % Sort features by gain for plotting
    [sorted_gain_plot, sort_idx_plot] = sort(feature_gain, 'descend');
    sorted_features_plot = feature_names(sort_idx_plot);
    
    % Create feature gain plot
    figure('Position', [300, 300, 1000, 600]);
    barh(sorted_gain_plot);
    set(gca, 'YTickLabel', sorted_features_plot);
    xlabel('Feature Gain (Accuracy Drop When Removed)', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel('Feature', 'FontSize', 14, 'FontWeight', 'bold')
    title(sprintf('Feature Gain Analysis - %s Model (Baseline Acc: %.2f%%)', ...
          upper(best_ml_model), baseline_acc * 100), ...
          'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    grid on; box on;
    
    % Add baseline reference line
    hold on;
    plot([0 0], ylim, 'r--', 'LineWidth', 2, 'DisplayName', 'No Gain (Baseline)');
    legend('Location', 'best');
    hold off;
    
    filename = sprintf('feature_gain_%s_ionic', best_ml_model);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    fprintf('\nFeature gain plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));
    
    % Update sorted_features and sorted_importance for use in rest of code
    sorted_features = sorted_features_plot;
    sorted_importance = sorted_gain_plot;  % Use actual gain values
    
    % Store feature gain in results
    ml_results.(best_ml_model).feature_gain = feature_gain;
    ml_results.(best_ml_model).feature_importance = sorted_gain_plot / sum(sorted_gain_plot);  % Normalized for compatibility
    ml_results.(best_ml_model).feature_names = feature_names;
    ml_results.(best_ml_model).sorted_features = sorted_features;
    ml_results.(best_ml_model).sorted_importance = sorted_importance;
    ml_results.(best_ml_model).baseline_accuracy = baseline_acc;
    
    % Describe most important features
    fprintf('\n========================================\n');
    fprintf('MOST IMPORTANT FEATURES (BY GAIN)\n');
    fprintf('========================================\n');
    fprintf('Top 3 most important features for predicting ln(gamma) sign:\n');
    for i = 1:min(3, length(sorted_features))
        gain_val = sorted_importance(i);
        gain_pct = (gain_val / baseline_acc) * 100;
        fprintf('  %d. %s (gain: %.4f, %.2f%% of baseline accuracy)\n', ...
                i, sorted_features{i}, gain_val, gain_pct);
    end
    
    % Provide interpretation
    fprintf('\nInterpretation:\n');
    top_feature = sorted_features{1};
    top_gain = sorted_importance(1);
    top_gain_pct = (top_gain / baseline_acc) * 100;
    fprintf('  The most important feature is "%s" with a gain of %.4f.\n', ...
            top_feature, top_gain);
    fprintf('  Removing this feature reduces model accuracy by %.2f%% (from %.2f%% to %.2f%%).\n', ...
            top_gain_pct, baseline_acc * 100, (baseline_acc - top_gain) * 100);
    
    if top_gain > baseline_acc * 0.1
        fprintf('  This feature is critical for model performance.\n');
    elseif top_gain > baseline_acc * 0.05
        fprintf('  This feature is highly important for model performance.\n');
    else
        fprintf('  Multiple features contribute to model performance.\n');
    end
    
    % Analyze top features in detail
    fprintf('\n========================================\n');
    fprintf('DETAILED FEATURE ANALYSIS\n');
    fprintf('========================================\n');
    
    % Get feature values for positive and negative classes
    for i = 1:min(3, length(sorted_features))
        feat_name = sorted_features{i};
        feat_idx = find(strcmp(feature_names, feat_name));
        
        if ~isempty(feat_idx)
            feat_values = X(:, feat_idx);
            feat_values_positive = feat_values(y == 1);
            feat_values_negative = feat_values(y == 0);
            
            % Find the gain for this feature
            feat_gain_idx = find(strcmp(sorted_features, feat_name));
            if ~isempty(feat_gain_idx)
                feat_gain = sorted_importance(feat_gain_idx);
            else
                feat_gain = 0;
            end
            fprintf('\n%s (Gain: %.4f, %.2f%% of baseline):\n', feat_name, feat_gain, (feat_gain / baseline_acc) * 100);
            if ~isempty(feat_values_positive)
                fprintf('  ln(gamma) > 0: mean=%.4f, std=%.4f, range=[%.4f, %.4f], n=%d\n', ...
                        mean(feat_values_positive), std(feat_values_positive), ...
                        min(feat_values_positive), max(feat_values_positive), ...
                        length(feat_values_positive));
            end
            if ~isempty(feat_values_negative)
                fprintf('  ln(gamma) <= 0: mean=%.4f, std=%.4f, range=[%.4f, %.4f], n=%d\n', ...
                        mean(feat_values_negative), std(feat_values_negative), ...
                        min(feat_values_negative), max(feat_values_negative), ...
                        length(feat_values_negative));
            end
            
            % Calculate separation
            if ~isempty(feat_values_positive) && ~isempty(feat_values_negative)
                mean_sep = abs(mean(feat_values_positive) - mean(feat_values_negative));
                pooled_std = sqrt((std(feat_values_positive)^2 + std(feat_values_negative)^2) / 2);
                if pooled_std > 0
                    cohens_d = mean_sep / pooled_std;
                    fprintf('  Effect size (Cohen''s d): %.4f', cohens_d);
                    if cohens_d > 0.8
                        fprintf(' (large effect)\n');
                    elseif cohens_d > 0.5
                        fprintf(' (medium effect)\n');
                    elseif cohens_d > 0.2
                        fprintf(' (small effect)\n');
                    else
                        fprintf(' (negligible effect)\n');
                    end
                end
            end
        end
    end
    
    % Model-specific insights
    fprintf('\n========================================\n');
    fprintf('MODEL-SPECIFIC INSIGHTS\n');
    fprintf('========================================\n');
    
    if strcmp(best_ml_model, 'tree')
        fprintf('Decision Tree Structure:\n');
        fprintf('  Number of nodes: %d\n', ml_results.tree.model.NumNodes);
        % Calculate tree depth
        try
            % Simple depth calculation: find maximum path length
            children = ml_results.tree.model.Children;
            depth = 0;
            stack = [1, 0];  % [node_index, current_depth]
            while ~isempty(stack)
                node = stack(1, 1);
                curr_depth = stack(1, 2);
                stack(1, :) = [];
                depth = max(depth, curr_depth);
                if children(node, 1) > 0  % Not a leaf
                    stack = [stack; children(node, 1), curr_depth + 1];
                    stack = [stack; children(node, 2), curr_depth + 1];
                end
            end
            fprintf('  Tree depth: %d\n', depth);
        catch
            fprintf('  Tree depth: N/A\n');
        end
        fprintf('  The tree makes decisions by splitting on features at each node.\n');
        fprintf('  Features used more frequently in splits are more important.\n');
    elseif strcmp(best_ml_model, 'rf')
        fprintf('Random Forest Structure:\n');
        fprintf('  Number of trees: %d\n', ml_results.rf.model.NumTrees);
        fprintf('  Features are ranked by their average importance across all trees.\n');
        fprintf('  Importance is measured by how much prediction error increases\n');
        fprintf('  when that feature is randomly permuted (OOB error).\n');
    elseif strcmp(best_ml_model, 'lr')
        fprintf('Logistic Regression Coefficients:\n');
        try
            coefficients = ml_results.lr.coefficients(2:end);  % Skip intercept
            intercept = ml_results.lr.coefficients(1);
            fprintf('  Intercept: %.4f\n', intercept);
            for i = 1:length(feature_names)
                fprintf('  %s: %.4f', feature_names{i}, coefficients(i));
                if coefficients(i) > 0
                    fprintf(' (increases probability of ln(gamma) > 0)\n');
                else
                    fprintf(' (decreases probability of ln(gamma) > 0)\n');
                end
            end
        catch
            fprintf('  Coefficient information not available\n');
        end
    end
end

%% ========================================================================
%% PLOT: ML PREDICTIONS
%% ========================================================================

if ~isempty(best_ml_model) && isfield(ml_results, best_ml_model)
    fprintf('\nGenerating ML prediction plot for Ionic Basis...\n');
    
    % Get predictions for all salts in original order
    y_pred_all = zeros(num_salts, 1);
    
    % Store train predictions in correct order
    for i = 1:length(train_idx)
        y_pred_all(train_idx(i)) = ml_results.(best_ml_model).y_pred_train(i);
    end
    
    % Store test predictions in correct order
    for i = 1:length(test_idx)
        y_pred_all(test_idx(i)) = ml_results.(best_ml_model).y_pred_test(i);
    end
    
    % Get test accuracy
    test_acc = ml_results.(best_ml_model).accuracy_test;
        
        figure('Position', [300, 300, 1400, 800]);
        hold on; grid on; box on;
        
    % Plot all salts with predictions
    for s = 1:num_salts
        salt_name = salt_names{s};
        data = all_salts_data.(salt_name);
        predicted = y_pred_all(s);
        actual = salt_labels.(salt_name);
        is_test = any(test_idx == s);
        
        if predicted == 1
            color = [0.8500, 0.3250, 0.0980];  % Orange for predicted positive
        else
            color = [0, 0.4470, 0.7410];  % Blue for predicted negative
        end
        
        if predicted == actual
            linestyle = '-';
            linewidth = 2;
        else
            linestyle = '--';
            linewidth = 3;  % Thicker for misclassified
        end
        
        % Make test set slightly more visible
        if is_test
            linewidth = linewidth + 0.5;
        end
        
        h = plot(data.RH*100, data.(ln_gamma_field), 'LineWidth', linewidth, ...
                 'color', color, 'LineStyle', linestyle, 'HandleVisibility', 'off');
        try
            if is_test
                h.Color(4) = 0.9;
            else
                h.Color(4) = 0.6;
            end
        catch
            if is_test
                h.Color = color * 0.9 + [1, 1, 1] * 0.1;
            else
                h.Color = color * 0.6 + [1, 1, 1] * 0.4;
            end
        end
    end
    
    plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'ln(\gamma_w) = 0')
    plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', [0.8500, 0.3250, 0.0980], ...
         'DisplayName', sprintf('Predicted ln(\\gamma) > 0 (Train: %d, Test: %d)', ...
         sum(y_pred_all(train_idx) == 1), sum(y_pred_all(test_idx) == 1)));
    plot(NaN, NaN, '-', 'LineWidth', 2.5, 'color', [0, 0.4470, 0.7410], ...
         'DisplayName', sprintf('Predicted ln(\\gamma) <= 0 (Train: %d, Test: %d)', ...
         sum(y_pred_all(train_idx) == 0), sum(y_pred_all(test_idx) == 0)));
    plot(NaN, NaN, '--', 'LineWidth', 3, 'color', 'r', ...
         'DisplayName', sprintf('Misclassified (Train: %d, Test: %d)', ...
         sum(y_pred_all(train_idx) ~= y_train), sum(y_pred_all(test_idx) ~= y_test)));
    
    xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold')
    ylabel('ln(\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
    title(sprintf('ln(\\gamma_w) vs RH - Ionic Basis - %s Predictions (Test Accuracy: %.2f%%)', ...
          upper(best_ml_model), test_acc * 100), ...
          'FontSize', 16, 'FontWeight', 'bold')
    legend('Location', 'best', 'FontSize', 12)
    xlim([0 100])
    set(gca, 'FontSize', 12)
    set(gcf, 'color', 'w');
    
    filename = sprintf('ln_gamma_vs_RH_ionic_ML_predictions_%s', best_ml_model);
    print(fullfile(fig_out_dir, filename), '-dpng', '-r600')
    fprintf('Plot saved to: %s\n', fullfile(fig_out_dir, [filename '.png']));
end

%% ========================================================================
%% SAVE RESULTS
%% ========================================================================

results_file = fullfile(fig_out_dir, 'prediction_analysis_results.txt');
fid = fopen(results_file, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'PREDICTION OF ln(gamma) SIGN AT HIGH RH\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Analysis Date: %s\n\n', datestr(now));
fprintf(fid, 'High RH Threshold: %.2f\n\n', high_RH_threshold);

fprintf(fid, 'Class Distribution - Ionic Basis:\n');
fprintf(fid, '  ln(gamma) > 0: %d salts\n', num_positive);
fprintf(fid, '  ln(gamma) <= 0: %d salts\n\n', num_negative);

fprintf(fid, '========================================\n');
fprintf(fid, 'CLASSIFICATION SCHEME PERFORMANCE - IONIC BASIS\n');
fprintf(fid, '========================================\n\n');

for rf = 1:length(result_fields)
    field = result_fields{rf};
    fprintf(fid, '--- %s ---\n', classification_scores.(field).title);
    fprintf(fid, '  Accuracy: %.4f\n', classification_scores.(field).accuracy);
    fprintf(fid, '  Precision: %.4f\n', classification_scores.(field).precision);
    fprintf(fid, '  Recall: %.4f\n', classification_scores.(field).recall);
    fprintf(fid, '  F1 Score: %.4f\n\n', classification_scores.(field).f1_score);
end

fprintf(fid, 'Best Classification - Ionic: %s\n', classification_scores.(best_classification_field).title);
fprintf(fid, '  Accuracy: %.4f\n', classification_scores.(best_classification_field).accuracy);
fprintf(fid, '  F1 Score: %.4f\n\n', classification_scores.(best_classification_field).f1_score);

if ~isempty(best_ml_model) && isfield(ml_results, best_ml_model)
    fprintf(fid, '========================================\n');
    fprintf(fid, 'MACHINE LEARNING RESULTS - IONIC BASIS\n');
    fprintf(fid, '========================================\n\n');
    fprintf(fid, 'Best Model: %s\n', best_ml_model);
    fprintf(fid, '  Training Accuracy: %.4f\n', ml_results.(best_ml_model).accuracy_train);
    fprintf(fid, '  Test Accuracy: %.4f\n\n', ml_results.(best_ml_model).accuracy_test);
    
    % Add confusion matrix if available
    if isfield(ml_results.(best_ml_model), 'confusion_train')
        fprintf(fid, '========================================\n');
        fprintf(fid, 'CONFUSION MATRIX ANALYSIS\n');
        fprintf(fid, '========================================\n\n');
        
        TP_train = ml_results.(best_ml_model).TP_train;
        TN_train = ml_results.(best_ml_model).TN_train;
        FP_train = ml_results.(best_ml_model).FP_train;
        FN_train = ml_results.(best_ml_model).FN_train;
        TP_test = ml_results.(best_ml_model).TP_test;
        TN_test = ml_results.(best_ml_model).TN_test;
        FP_test = ml_results.(best_ml_model).FP_test;
        FN_test = ml_results.(best_ml_model).FN_test;
        
        fprintf(fid, 'Training Set Confusion Matrix:\n');
        fprintf(fid, '                              Predicted\n');
        fprintf(fid, '                    ln(γ_w)>0 (γ_w>1)  ln(γ_w)≤0 (γ_w≤1)\n');
        fprintf(fid, 'Actual ln(γ_w)>0 (γ_w>1)      %3d              %3d\n', TP_train, FN_train);
        fprintf(fid, 'Actual ln(γ_w)≤0 (γ_w≤1)     %3d              %3d\n', FP_train, TN_train);
        fprintf(fid, '\n  True Positives (TP):  %d (correctly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', TP_train);
        fprintf(fid, '  True Negatives (TN):  %d (correctly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', TN_train);
        fprintf(fid, '  False Positives (FP): %d (incorrectly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', FP_train);
        fprintf(fid, '  False Negatives (FN): %d (incorrectly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n\n', FN_train);
        
        fprintf(fid, 'Test Set Confusion Matrix:\n');
        fprintf(fid, '                              Predicted\n');
        fprintf(fid, '                    ln(γ_w)>0 (γ_w>1)  ln(γ_w)≤0 (γ_w≤1)\n');
        fprintf(fid, 'Actual ln(γ_w)>0 (γ_w>1)      %3d              %3d\n', TP_test, FN_test);
        fprintf(fid, 'Actual ln(γ_w)≤0 (γ_w≤1)     %3d              %3d\n', FP_test, TN_test);
        fprintf(fid, '\n  True Positives (TP):  %d (correctly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', TP_test);
        fprintf(fid, '  True Negatives (TN):  %d (correctly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n', TN_test);
        fprintf(fid, '  False Positives (FP): %d (incorrectly predicted ln(gamma) > 0, i.e., gamma_w > 1)\n', FP_test);
        fprintf(fid, '  False Negatives (FN): %d (incorrectly predicted ln(gamma) <= 0, i.e., gamma_w <= 1)\n\n', FN_test);
        
        % Calculate metrics
        total_train = TP_train + TN_train + FP_train + FN_train;
        total_test = TP_test + TN_test + FP_test + FN_test;
        precision_train = TP_train / (TP_train + FP_train + eps);
        precision_test = TP_test / (TP_test + FP_test + eps);
        recall_train = TP_train / (TP_train + FN_train + eps);
        recall_test = TP_test / (TP_test + FN_test + eps);
        specificity_train = TN_train / (TN_train + FP_train + eps);
        specificity_test = TN_test / (TN_test + FP_test + eps);
        f1_train = 2 * (precision_train * recall_train) / (precision_train + recall_train + eps);
        f1_test = 2 * (precision_test * recall_test) / (precision_test + recall_test + eps);
        
        fprintf(fid, 'Training Set Metrics:\n');
        fprintf(fid, '  Accuracy:    %.4f (%.1f%%)\n', (TP_train + TN_train) / total_train, ...
                (TP_train + TN_train) / total_train * 100);
        fprintf(fid, '  Precision:   %.4f\n', precision_train);
        fprintf(fid, '  Recall:      %.4f\n', recall_train);
        fprintf(fid, '  Specificity: %.4f\n', specificity_train);
        fprintf(fid, '  F1 Score:    %.4f\n\n', f1_train);
        
        fprintf(fid, 'Test Set Metrics:\n');
        fprintf(fid, '  Accuracy:    %.4f (%.1f%%)\n', (TP_test + TN_test) / total_test, ...
                (TP_test + TN_test) / total_test * 100);
        fprintf(fid, '  Precision:   %.4f\n', precision_test);
        fprintf(fid, '  Recall:      %.4f\n', recall_test);
        fprintf(fid, '  Specificity: %.4f\n', specificity_test);
        fprintf(fid, '  F1 Score:    %.4f\n\n', f1_test);
    end
    
    % Add feature gain analysis if available
    if isfield(ml_results.(best_ml_model), 'feature_gain') || isfield(ml_results.(best_ml_model), 'feature_importance')
        if isfield(ml_results.(best_ml_model), 'baseline_accuracy')
            baseline_acc_file = ml_results.(best_ml_model).baseline_accuracy;
        else
            baseline_acc_file = ml_results.(best_ml_model).accuracy_test;
        end
        fprintf(fid, '========================================\n');
        fprintf(fid, 'FEATURE GAIN ANALYSIS\n');
        fprintf(fid, '========================================\n\n');
        fprintf(fid, 'Baseline Test Accuracy (all features): %.4f\n\n', baseline_acc_file);
        fprintf(fid, 'Gain = Baseline Accuracy - Accuracy Without Feature\n');
        fprintf(fid, 'Higher gain = more important feature\n\n');
        
        sorted_features = ml_results.(best_ml_model).sorted_features;
        sorted_gain = ml_results.(best_ml_model).sorted_importance;  % This contains gain values
        
        fprintf(fid, 'Feature Gain Rankings:\n');
        for i = 1:length(sorted_features)
            gain_val = sorted_gain(i);
            gain_pct = (gain_val / baseline_acc_file) * 100;
            fprintf(fid, '  %d. %s: Gain = %.4f (%.2f%% of baseline, acc without: %.4f)\n', ...
                    i, sorted_features{i}, gain_val, gain_pct, baseline_acc_file - gain_val);
        end
        
        fprintf(fid, '\nMost Important Feature: %s (Gain: %.4f, %.2f%% of baseline)\n\n', ...
                sorted_features{1}, sorted_gain(1), (sorted_gain(1) / baseline_acc_file) * 100);
        
        % Add detailed statistics for top 3 features
        fprintf(fid, 'Detailed Statistics for Top 3 Features:\n\n');
        feature_names = ml_results.(best_ml_model).feature_names;
        
        for i = 1:min(3, length(sorted_features))
            feat_name = sorted_features{i};
            feat_idx = find(strcmp(feature_names, feat_name));
            
            if ~isempty(feat_idx)
                feat_values = X(:, feat_idx);
                feat_values_positive = feat_values(y == 1);
                feat_values_negative = feat_values(y == 0);
                
                gain_val = sorted_gain(i);
                gain_pct = (gain_val / baseline_acc_file) * 100;
                fprintf(fid, '%s (Gain: %.4f, %.2f%% of baseline):\n', feat_name, gain_val, gain_pct);
                if ~isempty(feat_values_positive)
                    fprintf(fid, '  ln(gamma) > 0: mean=%.4f, std=%.4f, range=[%.4f, %.4f], n=%d\n', ...
                            mean(feat_values_positive), std(feat_values_positive), ...
                            min(feat_values_positive), max(feat_values_positive), ...
                            length(feat_values_positive));
                end
                if ~isempty(feat_values_negative)
                    fprintf(fid, '  ln(gamma) <= 0: mean=%.4f, std=%.4f, range=[%.4f, %.4f], n=%d\n', ...
                            mean(feat_values_negative), std(feat_values_negative), ...
                            min(feat_values_negative), max(feat_values_negative), ...
                            length(feat_values_negative));
                end
                
                % Calculate effect size
                if ~isempty(feat_values_positive) && ~isempty(feat_values_negative)
                    mean_sep = abs(mean(feat_values_positive) - mean(feat_values_negative));
                    pooled_std = sqrt((std(feat_values_positive)^2 + std(feat_values_negative)^2) / 2);
                    if pooled_std > 0
                        cohens_d = mean_sep / pooled_std;
                        fprintf(fid, '  Effect size (Cohen''s d): %.4f', cohens_d);
                        if cohens_d > 0.8
                            fprintf(fid, ' (large effect)\n');
                        elseif cohens_d > 0.5
                            fprintf(fid, ' (medium effect)\n');
                        elseif cohens_d > 0.2
                            fprintf(fid, ' (small effect)\n');
                        else
                            fprintf(fid, ' (negligible effect)\n');
                        end
                    end
                end
                fprintf(fid, '\n');
            end
        end
        
        % Model-specific information
        fprintf(fid, 'Model-Specific Information:\n');
        if strcmp(best_ml_model, 'tree')
            fprintf(fid, '  Model Type: Decision Tree\n');
            fprintf(fid, '  Number of nodes: %d\n', ml_results.tree.model.NumNodes);
            fprintf(fid, '  The tree splits on features to make decisions.\n');
            fprintf(fid, '  Features used more frequently in splits are more important.\n');
        elseif strcmp(best_ml_model, 'rf')
            fprintf(fid, '  Model Type: Random Forest\n');
            fprintf(fid, '  Number of trees: %d\n', ml_results.rf.model.NumTrees);
            fprintf(fid, '  Importance measured by OOB permuted predictor error.\n');
        elseif strcmp(best_ml_model, 'lr')
            fprintf(fid, '  Model Type: Logistic Regression\n');
            try
                coefficients = ml_results.lr.coefficients(2:end);
                intercept = ml_results.lr.coefficients(1);
                fprintf(fid, '  Intercept: %.4f\n', intercept);
                fprintf(fid, '  Coefficients (positive = increases prob of ln(gamma) > 0):\n');
                for i = 1:length(feature_names)
                    fprintf(fid, '    %s: %.4f\n', feature_names{i}, coefficients(i));
                end
            catch
                fprintf(fid, '  Coefficient information not available\n');
            end
        end
        
        % Add summary section
        fprintf(fid, '\n========================================\n');
        fprintf(fid, 'SUMMARY: KEY FINDINGS\n');
        fprintf(fid, '========================================\n\n');
        fprintf(fid, 'Best ML Model: %s (Test Accuracy: %.2f%%)\n\n', ...
                best_ml_model, ml_results.(best_ml_model).accuracy_test * 100);
        fprintf(fid, 'Most Important Features for Prediction (by Gain):\n');
        for i = 1:min(3, length(sorted_features))
            gain_val = sorted_gain(i);
            gain_pct = (gain_val / baseline_acc_file) * 100;
            fprintf(fid, '  %d. %s (Gain: %.4f, %.2f%% of baseline)\n', ...
                    i, sorted_features{i}, gain_val, gain_pct);
        end
        fprintf(fid, '\nThese features are the primary drivers of whether a salt will have\n');
        fprintf(fid, 'ln(gamma) > 0 or <= 0 at high relative humidity (RH >= %.2f).\n', high_RH_threshold);
        fprintf(fid, 'Gain represents how much model accuracy drops when each feature is removed.\n');
    end
end

fclose(fid);
fprintf('\nResults saved to: %s\n', results_file);

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');
