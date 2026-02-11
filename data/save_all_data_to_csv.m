close all 
clear
clc 

% Script to save all salt water activity data to CSV files
% Saves: molality, water activity, mass fraction (salt & water), 
%        mole fraction, activity coefficient, molecular weight, temperature

% Add calculate_mf, util, and data folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

T_default = 25; % Default temperature in Celsius
MWw = 18.015; % Molecular weight of water (g/mol)

% Load canonical salt data (exclude NH4NO3 - comment back in load_salt_data if needed)
salt_data = load_salt_data();
% Optionally exclude specific salts from CSV output (e.g. NH4NO3 was previously commented)
exclude_salts = {'NH4NO3', 'MgNO32'};  % MgNO32 is alias for MgNO3
keep = true(size(salt_data));
for s = 1:length(salt_data)
    if any(strcmp(salt_data{s}{1}, exclude_salts))
        keep(s) = false;
    end
end
salt_data = salt_data(keep);

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
