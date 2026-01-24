close all 
clear
clc 

% Script to save all salt water activity data to CSV files
% Saves: molality, water activity, mass fraction (salt & water), 
%        mole fraction, activity coefficient, molecular weight, temperature

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T_default = 25; % Default temperature in Celsius
MWw = 18.015; % Molecular weight of water (g/mol)

% Define all salts with their molecular weights and valid RH ranges
% Format: {salt_name, MW, RH_min, RH_max, function_name, function_args, salt_type, valence, temperature, nu, ionic_strength_factor}
% function_args: 0 = no args, 1 = needs T, 2 = other
% temperature: temperature in Celsius (if not specified, uses T_default)
% nu: number of ions the salt dissociates into (e.g., NaCl->2, CaCl2->3, Na2SO4->3)
% ionic_strength_factor: I = factor * molality (e.g., 1:1 salts->1, 2:1 or 1:2 salts->3)

salt_data = {
    % Endothermic salts (1:1 salts dissociate into 2 ions, I = m)
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1, 25, 2, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1, 25, 2, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1, 25, 2, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1, 25, 2, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1, 25, 2, 1};
    {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1, 25, 2, 1};
    {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1, 25, 2, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1, 25, 2, 1};
    % {'NH4NO3', 80.043, 0.118, 0.732, 'calculate_mf_NH4NO3', 0, 1, 1, 25, 2, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1, 25, 2, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1, 25, 2, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1, 25, 2, 1};
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1, 30, 2, 1};
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1, 30, 2, 1};
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1, 30, 2, 1};
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1, 30, 2, 1};
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1, 30, 2, 1};
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1, 30, 2, 1};
    
    % Exothermic salts (No source data provided to verify)
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1, 25, 2, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1, 25, 2, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1, 25, 2, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1, 25, 2, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2, 25, 3, 3};  % Ca2+ + 2Cl-, I = 3m
    {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2, 25, 3, 3};  % Mg2+ + 2Cl-, I = 3m
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2, 25, 3, 3};  % Mg2+ + 2NO3-, I = 3m
    {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1, 25, 2, 1};
    {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2, 25, 3, 3};  % Zn2+ + 2Cl-, I = 3m
    {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2, 25, 3, 3};  % Zn2+ + 2I-, I = 3m
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2, 25, 3, 3};  % Zn2+ + 2Br-, I = 3m
    {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1, 25, 2, 1};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1, 25, 3, 3};  % 2Na+ + SO4^2-, I = 3m
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1, 25, 3, 3};  % 2K+ + SO4^2-, I = 3m
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1, 25, 3, 3};  % 2NH4+ + SO4^2-, I = 3m
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1, 25, 2, 4};  % Mg2+ + SO4^2-, I = 4m
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1, 25, 2, 4};  % Mn2+ + SO4^2-, I = 4m
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1, 25, 3, 3};  % 2Li+ + SO4^2-, I = 3m
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1, 25, 2, 4};  % Ni2+ + SO4^2-, I = 4m
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1, 25, 2, 4};  % Cu2+ + SO4^2-, I = 4m
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1, 25, 2, 4};  % Zn2+ + SO4^2-, I = 4m
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2, 25, 3, 3};  % Ba2+ + 2NO3-, I = 3m
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2, 25, 3, 3};  % Ca2+ + 2NO3-, I = 3m
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2, 30, 3, 3}  % Ca2+ + 2Br-, I = 3m
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2, 30, 3, 3};  % Ca2+ + 2I-, I = 3m
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2, 30, 3, 3};  % Sr2+ + 2Cl-, I = 3m
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2, 30, 3, 3};  % Sr2+ + 2Br-, I = 3m
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2, 30, 3, 3};  % Sr2+ + 2I-, I = 3m
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2, 30, 3, 3};  % Ba2+ + 2Cl-, I = 3m
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2, 30, 3, 3};  % Ba2+ + 2Br-, I = 3m
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1, 25, 2, 1};
};

% Process each salt
num_points = 100;
data_output_dir = filepath; % Script is already in data/ folder

% Create output directory if it doesn't exist
if ~exist(data_output_dir, 'dir')
    mkdir(data_output_dir);
end

% Initialize cell array for combined data
all_data_combined = {};
all_data_combined{1, 1} = 'Salt';
all_data_combined{1, 2} = 'Temperature_C';
all_data_combined{1, 3} = 'MW_Salt_g_per_mol';
all_data_combined{1, 4} = 'RH_Water_Activity';
all_data_combined{1, 5} = 'Mass_Fraction_Salt';
all_data_combined{1, 6} = 'Mass_Fraction_Water';
all_data_combined{1, 7} = 'Mole_Fraction_Water';
all_data_combined{1, 8} = 'Mole_Fraction_Water_Ionic_Basis';
all_data_combined{1, 9} = 'Molality_mol_per_kg';
all_data_combined{1, 10} = 'Activity_Coefficient_Water';
all_data_combined{1, 11} = 'Ionic_Strength_mol_per_kg';
all_data_combined{1, 12} = 'ln_aw';
all_data_combined{1, 13} = 'ln_gammaw';

row_idx = 2; % Start from row 2 (after header)

fprintf('\n=== Saving Water Activity Data to CSV ===\n');
fprintf('Output directory: %s\n\n', data_output_dir);

for s = 1:length(salt_data)
    salt_name = salt_data{s}{1};
    MW = salt_data{s}{2};
    RH_min = salt_data{s}{3};
    RH_max = salt_data{s}{4};
    func_name = salt_data{s}{5};
    func_args = salt_data{s}{6};
    
    % Get temperature for this salt (use default if not specified)
    if length(salt_data{s}) >= 9
        T = salt_data{s}{9};
    else
        T = T_default;
    end
    
    % Get dissociation number (nu) for this salt
    if length(salt_data{s}) >= 10
        nu = salt_data{s}{10};
    else
        nu = 2;  % Default to 2 for simple 1:1 salts
    end
    
    % Get ionic strength factor for this salt
    if length(salt_data{s}) >= 11
        I_factor = salt_data{s}{11};
    else
        I_factor = 1;  % Default to 1 for simple 1:1 salts
    end
    
    fprintf('Processing %s (T=%dC, nu=%d, I_factor=%d)...', salt_name, T, nu, I_factor);
    
    % Create RH vector
    RH_vec_full = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays (use full size first, will trim later)
    RH_vec = zeros(size(RH_vec_full));
    mf_salt = zeros(size(RH_vec_full));
    mf_water = zeros(size(RH_vec_full));
    x_water = zeros(size(RH_vec_full));  % Mole fraction treating salt as one mole
    x_water_ionic = zeros(size(RH_vec_full));  % Mole fraction on ionic basis
    molality = zeros(size(RH_vec_full));
    ionic_strength = zeros(size(RH_vec_full));
    
    % Calculate for each RH value
    valid_count = 0;
    for i = 1:length(RH_vec_full)
        try
            if func_args == 1
                % Function needs temperature
                mf_s = feval(func_name, RH_vec_full(i), T);
            else
                % Function doesn't need temperature
                mf_s = feval(func_name, RH_vec_full(i));
            end
            mf_w = 1 - mf_s;
            
            % Validate physical constraints
            if mf_s < 0 || mf_s > 1 || mf_w < 0 || mf_w > 1
                % Skip this point - outside valid range
                continue;
            end
            
            valid_count = valid_count + 1;
            RH_vec(valid_count) = RH_vec_full(i);
            mf_salt(valid_count) = mf_s;
            mf_water(valid_count) = mf_w;
            
            % Calculate mole fraction of water (treating salt as one mole)
            x_water(valid_count) = (mf_w / MWw) / ...
                ((mf_w / MWw) + (mf_s / MW));
            
            % Calculate mole fraction of water on ionic basis
            % One mole of salt dissociates into nu moles of ions
            x_water_ionic(valid_count) = (mf_w / MWw) / ...
                ((mf_w / MWw) + nu * (mf_s / MW));
            
            % Calculate molality (mol solute / kg water)
            % molality = (mass_salt / MW_salt) / (mass_water / 1000)
            % For 1 kg of solution: molality = (mf_salt / MW) / (mf_water / 1000)
            molality(valid_count) = (mf_s / MW) / (mf_w / 1000);
            
            % Calculate ionic strength: I = I_factor * molality
            ionic_strength(valid_count) = I_factor * molality(valid_count);
            
        catch ME
            % Skip this point if there's an error
            continue;
        end
    end
    
    % Trim arrays to valid data only
    if valid_count > 0
        RH_vec = RH_vec(1:valid_count);
        mf_salt = mf_salt(1:valid_count);
        mf_water = mf_water(1:valid_count);
        x_water = x_water(1:valid_count);
        x_water_ionic = x_water_ionic(1:valid_count);
        molality = molality(1:valid_count);
        ionic_strength = ionic_strength(1:valid_count);
        num_points = valid_count;
        success = true;
    else
        success = false;
    end
    
    if success
        % Calculate Activity Coefficient (gamma_w = a_w / x_w)
        % Use ionic basis for activity coefficient (more physically meaningful)
        gamma_w = RH_vec ./ x_water_ionic;
        
        % Calculate natural logarithms
        ln_aw = log(RH_vec);
        ln_gammaw = log(gamma_w);
        
        % Add to combined data array
        for i = 1:num_points
            all_data_combined{row_idx, 1} = salt_name;
            all_data_combined{row_idx, 2} = T;
            all_data_combined{row_idx, 3} = MW;
            all_data_combined{row_idx, 4} = RH_vec(i);
            all_data_combined{row_idx, 5} = mf_salt(i);
            all_data_combined{row_idx, 6} = mf_water(i);
            all_data_combined{row_idx, 7} = x_water(i);
            all_data_combined{row_idx, 8} = x_water_ionic(i);
            all_data_combined{row_idx, 9} = molality(i);
            all_data_combined{row_idx, 10} = gamma_w(i);
            all_data_combined{row_idx, 11} = ionic_strength(i);
            all_data_combined{row_idx, 12} = ln_aw(i);
            all_data_combined{row_idx, 13} = ln_gammaw(i);
            row_idx = row_idx + 1;
        end
        
        fprintf(' ✓ Saved (%d points)\n', num_points);
    else
        fprintf(' ✗ Failed\n');
    end
end

% Save combined CSV file
combined_filename = fullfile(data_output_dir, 'water_activity_all_salts_combined.csv');
fid = fopen(combined_filename, 'w');

% Write header
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
    all_data_combined{1,1}, all_data_combined{1,2}, all_data_combined{1,3}, ...
    all_data_combined{1,4}, all_data_combined{1,5}, all_data_combined{1,6}, ...
    all_data_combined{1,7}, all_data_combined{1,8}, all_data_combined{1,9}, ...
    all_data_combined{1,10}, all_data_combined{1,11}, all_data_combined{1,12}, ...
    all_data_combined{1,13});

% Write data rows
for i = 2:size(all_data_combined, 1)
    fprintf(fid, '%s,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
        all_data_combined{i,1}, all_data_combined{i,2}, all_data_combined{i,3}, ...
        all_data_combined{i,4}, all_data_combined{i,5}, all_data_combined{i,6}, ...
        all_data_combined{i,7}, all_data_combined{i,8}, all_data_combined{i,9}, ...
        all_data_combined{i,10}, all_data_combined{i,11}, all_data_combined{i,12}, ...
        all_data_combined{i,13});
end

fclose(fid);

fprintf('\n=== Summary ===\n');
fprintf('Salts processed: %d\n', s);
fprintf('Combined file: water_activity_all_salts_combined.csv\n');
fprintf('Total data points: %d\n', row_idx - 2);
fprintf('\nFile saved in: %s\n', data_output_dir);
