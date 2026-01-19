close all 
clear
clc 

% Script to save all salt water activity data to CSV files
% Saves: molality, water activity, mass fraction (salt & water), 
%        mole fraction, activity coefficient, molecular weight, temperature

% Add calculate_mf and util folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, 'calculate_mf'));
addpath(fullfile(filepath, 'util'));

T = 25; % Temperature in Celsius
MWw = 18.015; % Molecular weight of water (g/mol)

% Define all salts with their molecular weights and valid RH ranges
% Format: {salt_name, MW, RH_min, RH_max, function_name, function_args}
% function_args: 0 = no args, 1 = needs T, 2 = other

salt_data = {
    % Endothermic salts (RH ranges matched to Isotherms_screening_endothermic.m)
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl_', 0};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl_', 0};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl_', 0};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl_', 0};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3_', 0};
    {'AgNO3', 169.87, 0.865, 0.99, 'calculate_mf_AgNO3_', 0};
    {'KI', 165.998, 0.975, 0.995, 'calculate_mf_KI_', 0};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0};
    {'KClO3', 122.55, 0.981, 0.995, 'calculate_mf_KClO3', 0};
    {'NaBr', 102.89, 0.614, 0.99, 'calculate_mf_NaBr', 0};
    {'NaI', 149.89, 0.581, 0.99, 'calculate_mf_NaI', 0};
    {'KBr', 119.00, 0.833, 0.99, 'calculate_mf_KBr', 0};
    {'RbCl', 120.92, 0.743, 0.99, 'calculate_mf_RbCl', 0};
    {'CsBr', 212.81, 0.848, 0.99, 'calculate_mf_CsBr', 0};
    {'CsI', 259.81, 0.913, 0.995, 'calculate_mf_CsI', 0};
    
    % Exothermic salts
    {'LiCl', 42.4, 0.12, 0.9, 'calculate_mf_LiCl', 1};
    {'LiOH', 24, 0.85, 0.9, 'calculate_mf_LiOH', 0};
    {'NaOH', 40, 0.23, 0.9, 'calculate_mf_NaOH', 0};
    {'HCl', 36.5, 0.17, 0.9, 'calculate_mf_HCl', 0};
    {'CaCl2', 111, 0.31, 0.9, 'calculate_mf_CaCl', 1};
    {'MgCl2', 95.2, 0.33, 0.9, 'calculate_mf_MgCl', 0};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0};
    {'LiBr', 86.85, 0.07, 0.9, 'calculate_mf_LiBr', 0};
    {'ZnCl2', 136.3, 0.07, 0.8, 'calculate_mf_ZnCl', 0};
    {'ZnI2', 319.18, 0.25, 0.9, 'calculate_mf_ZnI', 0};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0};
    {'LiI', 133.85, 0.18, 0.9, 'calculate_mf_LiI', 0};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.8990, 0.9957, 'calculate_mf_Na2SO4_', 0};
    {'K2SO4', 174.26, 0.9720, 0.9958, 'calculate_mf_K2SO4_', 0};
    {'NH42SO4', 132.14, 0.8310, 0.9959, 'calculate_mf_NH42SO4_', 0};
    {'MgSO4', 120.37, 0.9050, 0.9960, 'calculate_mf_MgSO4_', 0};
    {'MnSO4', 151.00, 0.8620, 0.9961, 'calculate_mf_MnSO4_', 0};
    {'Li2SO4', 109.94, 0.8530, 0.9956, 'calculate_mf_Li2SO4_', 0};
    {'NiSO4', 154.75, 0.9390, 0.9962, 'calculate_mf_NiSO4_', 0};
    {'CuSO4', 159.61, 0.9750, 0.9963, 'calculate_mf_CuSO4_', 0};
    {'ZnSO4', 161.44, 0.9130, 0.9962, 'calculate_mf_ZnSO4_', 0};
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9859, 0.9958, 'calculate_mf_BaNO3', 0};
    {'CaNO3', 164.09, 0.6464, 0.9955, 'calculate_mf_CaNO3', 0};
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6395, 0.9540, 'calculate_mf_CaBr2', 0};
    {'CaI2', 293.89, 0.8321, 0.9524, 'calculate_mf_CaI2', 0};
    {'SrCl2', 158.53, 0.8059, 0.9778, 'calculate_mf_SrCl2', 0};
    {'SrBr2', 247.43, 0.7776, 0.9571, 'calculate_mf_SrBr2', 0};
    {'SrI2', 341.43, 0.6785, 0.9569, 'calculate_mf_SrI2', 0};
    {'BaCl2', 208.23, 0.9375, 0.9731, 'calculate_mf_BaCl2', 0};
    {'BaBr2', 297.14, 0.8221, 0.9587, 'calculate_mf_BaBr2', 0};
    
    % Chlorates
    {'LiClO4', 106.39, 0.7775, 0.9869, 'calculate_mf_LiClO4', 0};
};

% Process each salt
num_points = 100;
data_output_dir = fullfile(filepath, 'data');

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
all_data_combined{1, 8} = 'Molality_mol_per_kg';
all_data_combined{1, 9} = 'Activity_Coefficient_Water';

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
    
    fprintf('Processing %s...', salt_name);
    
    % Create RH vector
    RH_vec_full = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays (use full size first, will trim later)
    RH_vec = zeros(size(RH_vec_full));
    mf_salt = zeros(size(RH_vec_full));
    mf_water = zeros(size(RH_vec_full));
    x_water = zeros(size(RH_vec_full));
    molality = zeros(size(RH_vec_full));
    
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
            
            % Calculate mole fraction of water
            x_water(valid_count) = (mf_w / MWw) / ...
                ((mf_w / MWw) + (mf_s / MW));
            
            % Calculate molality (mol solute / kg water)
            % molality = (mass_salt / MW_salt) / (mass_water / 1000)
            % For 1 kg of solution: molality = (mf_salt / MW) / (mf_water / 1000)
            molality(valid_count) = (mf_s / MW) / (mf_w / 1000);
            
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
        molality = molality(1:valid_count);
        num_points = valid_count;
        success = true;
    else
        success = false;
    end
    
    if success
        % Calculate Activity Coefficient (gamma_w = a_w / x_w)
        gamma_w = RH_vec ./ x_water;
        
        % Save individual CSV file for this salt
        individual_filename = fullfile(data_output_dir, sprintf('water_activity_%s.csv', salt_name));
        
        % Create data matrix
        data_matrix = [repmat(T, num_points, 1), ...
                      repmat(MW, num_points, 1), ...
                      RH_vec', ...
                      mf_salt', ...
                      mf_water', ...
                      x_water', ...
                      molality', ...
                      gamma_w'];
        
        % Write to individual CSV file with header
        fid = fopen(individual_filename, 'w');
        fprintf(fid, 'Temperature_C,MW_Salt_g_per_mol,RH_Water_Activity,Mass_Fraction_Salt,Mass_Fraction_Water,Mole_Fraction_Water,Molality_mol_per_kg,Activity_Coefficient_Water\n');
        fclose(fid);
        dlmwrite(individual_filename, data_matrix, '-append', 'delimiter', ',', 'precision', '%.10f');
        
        % Add to combined data array
        for i = 1:num_points
            all_data_combined{row_idx, 1} = salt_name;
            all_data_combined{row_idx, 2} = T;
            all_data_combined{row_idx, 3} = MW;
            all_data_combined{row_idx, 4} = RH_vec(i);
            all_data_combined{row_idx, 5} = mf_salt(i);
            all_data_combined{row_idx, 6} = mf_water(i);
            all_data_combined{row_idx, 7} = x_water(i);
            all_data_combined{row_idx, 8} = molality(i);
            all_data_combined{row_idx, 9} = gamma_w(i);
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
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
    all_data_combined{1,1}, all_data_combined{1,2}, all_data_combined{1,3}, ...
    all_data_combined{1,4}, all_data_combined{1,5}, all_data_combined{1,6}, ...
    all_data_combined{1,7}, all_data_combined{1,8}, all_data_combined{1,9});

% Write data rows
for i = 2:size(all_data_combined, 1)
    fprintf(fid, '%s,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n', ...
        all_data_combined{i,1}, all_data_combined{i,2}, all_data_combined{i,3}, ...
        all_data_combined{i,4}, all_data_combined{i,5}, all_data_combined{i,6}, ...
        all_data_combined{i,7}, all_data_combined{i,8}, all_data_combined{i,9});
end

fclose(fid);

fprintf('\n=== Summary ===\n');
fprintf('Individual salt files: %d\n', s);
fprintf('Combined file: water_activity_all_salts_combined.csv\n');
fprintf('Total data points: %d\n', row_idx - 2);
fprintf('\nAll files saved in: %s\n', data_output_dir);
