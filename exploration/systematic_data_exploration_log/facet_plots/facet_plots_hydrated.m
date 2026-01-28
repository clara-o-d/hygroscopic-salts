close all
clear
clc

% =========================================================================
% SYSTEMATIC EXPLORATION OF WATER ACTIVITY DATA (HYDRATED RADII VERSION)
% =========================================================================
% This script performs comprehensive analysis using HYDRATED ionic radii
% from baseline_numeric_only.csv instead of Shannon ionic radii.
% =========================================================================

%% Setup
fprintf('\n========================================\n');
fprintf('SYSTEMATIC WATER ACTIVITY DATA EXPLORATION (LOGARITHMIC)\n');
fprintf('========================================\n\n');

% Add paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'util'));

% Open output file for writing analysis results
output_filename = fullfile(filepath, 'systematic_exploration_results_log_hydrated.txt');
fid_out = fopen(output_filename, 'w');
if fid_out == -1
    warning('Could not open output file for writing. Only printing to console.');
    fid_out = 1;  % Use stdout
end

% Define helper function to write to both console and file
if fid_out == 1
    write_output = @(varargin) fprintf(varargin{:});
else
    write_output = @(varargin) cellfun(@(x) fprintf(x, varargin{:}), {1, fid_out}, 'UniformOutput', false);
end

% Write header to file
write_output('\n========================================\n');
write_output('SYSTEMATIC WATER ACTIVITY DATA EXPLORATION (LOGARITHMIC)\n');
write_output('========================================\n');
write_output('Date: %s\n', datestr(now));
write_output('========================================\n\n');

% Constants
MWw = 18.015; % g/mol
R = 8.314; % J/(mol·K)
T_kelvin = 298.15; % 25°C

%% Load data
write_output('Loading data...\n');
data_file = fullfile(filepath, '../../..', 'data', 'water_activity_all_salts_combined.csv');
data_table = readtable(data_file);

% Extract unique salts
unique_salts = unique(data_table.Salt);
n_salts = length(unique_salts);
write_output('Loaded data for %d salts\n', n_salts);
write_output('Total data points: %d\n\n', height(data_table));

%% Define comprehensive salt properties database
write_output('Building salt properties database...\n');

% Hydrated ionic radii (pm) - from baseline_numeric_only.csv
% Converted from Angstroms to pm (×100)
% Cations
hydrated_radii = struct();
hydrated_radii.Ag = 341.0;  % Ag+
hydrated_radii.Ba = 404.0;  % Ba+
hydrated_radii.Ca = 412.0;  % Ca+
hydrated_radii.Cs = 329.0;  % Cs+
hydrated_radii.Cu = 358.0;  % Cu+
hydrated_radii.H = 358.0;  % H+
hydrated_radii.K = 331.0;  % K+
hydrated_radii.Li = 382.0;  % Li+
hydrated_radii.Mg = 428.0;  % Mg+
hydrated_radii.Mn = 438.0;  % Mn+
hydrated_radii.NH4 = 331.0;  % NH4+
hydrated_radii.Na = 358.0;  % Na+
hydrated_radii.Ni = 404.0;  % Ni+
hydrated_radii.Rb = 329.0;  % Rb+
hydrated_radii.Sr = 412.0;  % Sr+
hydrated_radii.Zn = 430.0;  % Zn+

% Anions
hydrated_radii.Br = 330.0;  % Br-
hydrated_radii.Cl = 332.0;  % Cl-
hydrated_radii.ClO3 = 341.0;  % ClO3-
hydrated_radii.ClO4 = 338.0;  % ClO4-
hydrated_radii.I = 331.0;  % I-
hydrated_radii.NO3 = 335.0;  % NO3-
hydrated_radii.OH = 300.0;  % OH-
hydrated_radii.SO4 = 379.0;  % SO4-
 % ClO4- (effective)

% Salt properties: [cation, cation_charge, cation_stoich, anion, anion_charge, 
%                   anion_stoich, enthalpy_solution_kJ_mol, category, 
%                   cation_type, anion_type]
% Types: 'k'=kosmotrope, 'c'=chaotrope, 'u'=unknown
salt_props = containers.Map();

% Alkali metal chlorides
salt_props('LiCl') = {'Li', 1, 1, 'Cl', 1, 1, -37.0, 'exothermic_halide', 'k', 'c'};
salt_props('NaCl') = {'Na', 1, 1, 'Cl', 1, 1, -3.9, 'exothermic_halide', 'c', 'c'};
salt_props('KCl') = {'K', 1, 1, 'Cl', 1, 1, 17.2, 'endothermic_halide', 'c', 'c'};
salt_props('RbCl') = {'Rb', 1, 1, 'Cl', 1, 1, 17.3, 'endothermic_halide', 'c', 'c'};
salt_props('CsCl') = {'Cs', 1, 1, 'Cl', 1, 1, 17.8, 'endothermic_halide', 'c', 'c'};
salt_props('NH4Cl') = {'NH4', 1, 1, 'Cl', 1, 1, 14.8, 'endothermic_halide', 'c', 'c'};

% Alkali metal bromides
salt_props('LiBr') = {'Li', 1, 1, 'Br', 1, 1, -48.8, 'exothermic_halide', 'k', 'c'};
salt_props('NaBr') = {'Na', 1, 1, 'Br', 1, 1, -0.6, 'endothermic_halide', 'c', 'c'};
salt_props('KBr') = {'K', 1, 1, 'Br', 1, 1, 19.9, 'endothermic_halide', 'c', 'c'};
salt_props('CsBr') = {'Cs', 1, 1, 'Br', 1, 1, 25.0, 'endothermic_halide', 'c', 'c'};

% Alkali metal iodides
salt_props('LiI') = {'Li', 1, 1, 'I', 1, 1, -63.3, 'exothermic_halide', 'k', 'c'};
salt_props('NaI') = {'Na', 1, 1, 'I', 1, 1, -7.5, 'endothermic_halide', 'c', 'c'};
salt_props('KI') = {'K', 1, 1, 'I', 1, 1, 20.3, 'endothermic_halide', 'c', 'c'};
salt_props('CsI') = {'Cs', 1, 1, 'I', 1, 1, 33.3, 'endothermic_halide', 'c', 'c'};

% Alkaline earth halides
salt_props('MgCl2') = {'Mg', 2, 1, 'Cl', 1, 2, -160.0, 'exothermic_halide', 'k', 'c'};
salt_props('CaCl2') = {'Ca', 2, 1, 'Cl', 1, 2, -82.8, 'exothermic_halide', 'k', 'c'};
salt_props('SrCl2') = {'Sr', 2, 1, 'Cl', 1, 2, -47.7, 'exothermic_halide', 'k', 'c'};
salt_props('BaCl2') = {'Ba', 2, 1, 'Cl', 1, 2, -8.7, 'endothermic_halide', 'k', 'c'};

salt_props('CaBr2') = {'Ca', 2, 1, 'Br', 1, 2, -103.1, 'exothermic_halide', 'k', 'c'};
salt_props('SrBr2') = {'Sr', 2, 1, 'Br', 1, 2, -71.1, 'exothermic_halide', 'k', 'c'};
salt_props('BaBr2') = {'Ba', 2, 1, 'Br', 1, 2, -27.6, 'exothermic_halide', 'k', 'c'};

salt_props('CaI2') = {'Ca', 2, 1, 'I', 1, 2, -110.0, 'exothermic_halide', 'k', 'c'};
salt_props('SrI2') = {'Sr', 2, 1, 'I', 1, 2, -89.5, 'exothermic_halide', 'k', 'c'};

% Transition metal halides
salt_props('ZnCl2') = {'Zn', 2, 1, 'Cl', 1, 2, -65.3, 'exothermic_halide', 'k', 'c'};
salt_props('ZnBr2') = {'Zn', 2, 1, 'Br', 1, 2, -72.4, 'exothermic_halide', 'k', 'c'};
salt_props('ZnI2') = {'Zn', 2, 1, 'I', 1, 2, -52.3, 'exothermic_halide', 'k', 'c'};

% Nitrates
salt_props('LiNO3') = {'Li', 1, 1, 'NO3', 1, 1, -2.5, 'endothermic_nitrate', 'k', 'c'};
salt_props('NaNO3') = {'Na', 1, 1, 'NO3', 1, 1, 20.5, 'endothermic_nitrate', 'c', 'c'};
salt_props('KNO3') = {'K', 1, 1, 'NO3', 1, 1, 34.9, 'endothermic_nitrate', 'c', 'c'};
salt_props('NH4NO3') = {'NH4', 1, 1, 'NO3', 1, 1, 25.7, 'endothermic_nitrate', 'c', 'c'};
salt_props('AgNO3') = {'Ag', 1, 1, 'NO3', 1, 1, 22.6, 'endothermic_nitrate', 'c', 'c'};
salt_props('Ca(NO3)2') = {'Ca', 2, 1, 'NO3', 1, 2, -19.2, 'exothermic_nitrate', 'k', 'c'};
salt_props('Mg(NO3)2') = {'Mg', 2, 1, 'NO3', 1, 2, -90.9, 'exothermic_nitrate', 'k', 'c'};
salt_props('Ba(NO3)2') = {'Ba', 2, 1, 'NO3', 1, 2, 42.0, 'endothermic_nitrate', 'k', 'c'};

% Sulfates
salt_props('Li2SO4') = {'Li', 1, 2, 'SO4', 2, 1, -29.8, 'endothermic_sulfate', 'k', 'c'};
salt_props('Na2SO4') = {'Na', 1, 2, 'SO4', 2, 1, -2.4, 'endothermic_sulfate', 'c', 'c'};
salt_props('K2SO4') = {'K', 1, 2, 'SO4', 2, 1, 23.8, 'endothermic_sulfate', 'c', 'c'};
salt_props('(NH4)2SO4') = {'NH4', 1, 2, 'SO4', 2, 1, 6.6, 'endothermic_sulfate', 'c', 'c'};
salt_props('MgSO4') = {'Mg', 2, 1, 'SO4', 2, 1, -91.2, 'endothermic_sulfate', 'k', 'c'};
salt_props('CuSO4') = {'Cu', 2, 1, 'SO4', 2, 1, -66.5, 'endothermic_sulfate', 'u', 'c'}; % Cu+2 not on list
salt_props('ZnSO4') = {'Zn', 2, 1, 'SO4', 2, 1, -78.9, 'endothermic_sulfate', 'k', 'c'};
salt_props('MnSO4') = {'Mn', 2, 1, 'SO4', 2, 1, -54.0, 'endothermic_sulfate', 'k', 'c'};
salt_props('NiSO4') = {'Ni', 2, 1, 'SO4', 2, 1, -73.2, 'endothermic_sulfate', 'k', 'c'};

% Hydroxides and strong acid
salt_props('LiOH') = {'Li', 1, 1, 'OH', 1, 1, -23.6, 'exothermic_base', 'k', 'k'};
salt_props('NaOH') = {'Na', 1, 1, 'OH', 1, 1, -44.5, 'exothermic_base', 'c', 'k'};
salt_props('HCl') = {'H', 1, 1, 'Cl', 1, 1, -74.8, 'exothermic_acid', 'k', 'c'};

% Chlorates/Perchlorates
salt_props('KClO3') = {'K', 1, 1, 'ClO3', 1, 1, 41.4, 'endothermic_oxyanion', 'c', 'c'};
salt_props('NaClO4') = {'Na', 1, 1, 'ClO4', 1, 1, 13.9, 'endothermic_oxyanion', 'c', 'c'};
salt_props('LiClO4') = {'Li', 1, 1, 'ClO4', 1, 1, -32.6, 'exothermic_oxyanion', 'k', 'c'};

write_output('Salt properties database built (%d salts)\n\n', salt_props.Count);

% Define valid RH ranges for each salt (from save_all_data_to_csv.m)
% Format: {salt_name, RH_min, RH_max}
salt_RH_ranges = containers.Map();
salt_RH_ranges('NaCl') = [0.765, 0.99];
salt_RH_ranges('KCl') = [0.855, 0.99];
salt_RH_ranges('NH4Cl') = [0.815, 0.99];
salt_RH_ranges('CsCl') = [0.82, 0.99];
salt_RH_ranges('NaNO3') = [0.971, 0.995];
salt_RH_ranges('AgNO3') = [0.865, 0.99];
salt_RH_ranges('KI') = [0.975, 0.995];
salt_RH_ranges('LiNO3') = [0.736, 0.99];
salt_RH_ranges('KNO3') = [0.932, 0.995];
salt_RH_ranges('NaClO4') = [0.778, 0.99];
salt_RH_ranges('KClO3') = [0.981, 0.995];
salt_RH_ranges('NaBr') = [0.614, 0.99];
salt_RH_ranges('NaI') = [0.581, 0.99];
salt_RH_ranges('KBr') = [0.833, 0.99];
salt_RH_ranges('RbCl') = [0.743, 0.99];
salt_RH_ranges('CsBr') = [0.848, 0.99];
salt_RH_ranges('CsI') = [0.913, 0.995];
salt_RH_ranges('LiCl') = [0.12, 0.9];
salt_RH_ranges('LiOH') = [0.85, 0.9];
salt_RH_ranges('NaOH') = [0.23, 0.9];
salt_RH_ranges('HCl') = [0.17, 0.9];
salt_RH_ranges('CaCl2') = [0.31, 0.9];
salt_RH_ranges('MgCl2') = [0.33, 0.9];
salt_RH_ranges('MgNO3') = [0.55, 0.9];
salt_RH_ranges('LiBr') = [0.07, 0.9];
salt_RH_ranges('ZnCl2') = [0.07, 0.8];
salt_RH_ranges('ZnI2') = [0.25, 0.9];
salt_RH_ranges('ZnBr2') = [0.08, 0.85];
salt_RH_ranges('LiI') = [0.18, 0.9];
salt_RH_ranges('Na2SO4') = [0.8990, 0.9957];
salt_RH_ranges('K2SO4') = [0.9720, 0.9958];
salt_RH_ranges('NH42SO4') = [0.8310, 0.9959];
salt_RH_ranges('MgSO4') = [0.9050, 0.9960];
salt_RH_ranges('MnSO4') = [0.8620, 0.9961];
salt_RH_ranges('Li2SO4') = [0.8530, 0.9956];
salt_RH_ranges('NiSO4') = [0.9390, 0.9962];
salt_RH_ranges('CuSO4') = [0.9750, 0.9963];
salt_RH_ranges('ZnSO4') = [0.9130, 0.9962];
salt_RH_ranges('NH4NO3') = [0.118, 0.732];
salt_RH_ranges('BaNO3') = [0.9859, 0.9958];
salt_RH_ranges('CaNO3') = [0.6464, 0.9955];
salt_RH_ranges('CaBr2') = [0.6395, 0.9540];
salt_RH_ranges('CaI2') = [0.8321, 0.9524];
salt_RH_ranges('SrCl2') = [0.8059, 0.9778];
salt_RH_ranges('SrBr2') = [0.7776, 0.9571];
salt_RH_ranges('SrI2') = [0.6785, 0.9569];
salt_RH_ranges('BaCl2') = [0.9375, 0.9731];
salt_RH_ranges('BaBr2') = [0.8221, 0.9587];
salt_RH_ranges('LiClO4') = [0.7775, 0.9869];

%% Calculate derived properties for each salt
write_output('Calculating derived properties...\n');

% Initialize storage
salt_analysis = struct();

for i = 1:n_salts
    salt_name = unique_salts{i};
    
    % Get data for this salt
    idx = strcmp(data_table.Salt, salt_name);
    salt_data = data_table(idx, :);
    
    % Basic properties
    salt_analysis(i).name = salt_name;
    salt_analysis(i).MW = salt_data.MW_Salt_g_per_mol(1);
    salt_analysis(i).n_points = sum(idx);
    
    % Extract properties if available in database
    if isKey(salt_props, salt_name)
        props = salt_props(salt_name);
        salt_analysis(i).cation = props{1};
        salt_analysis(i).cat_charge = props{2};
        salt_analysis(i).cat_stoich = props{3};
        salt_analysis(i).anion = props{4};
        salt_analysis(i).an_charge = props{5};
        salt_analysis(i).an_stoich = props{6};
        salt_analysis(i).delta_H_soln = props{7};
        salt_analysis(i).category = props{8};
        salt_analysis(i).cat_type = props{9};
        salt_analysis(i).an_type = props{10};
        salt_analysis(i).kc_pair_type = [salt_analysis(i).cat_type '-' salt_analysis(i).an_type];

        % Get ionic radii
        if isfield(hydrated_radii, salt_analysis(i).cation)
            salt_analysis(i).cat_radius = hydrated_radii.(salt_analysis(i).cation);
        else
            salt_analysis(i).cat_radius = NaN;
        end
        
        if isfield(hydrated_radii, salt_analysis(i).anion)
            salt_analysis(i).an_radius = hydrated_radii.(salt_analysis(i).anion);
        else
            salt_analysis(i).an_radius = NaN;
        end
        
        % Calculate charge densities (e/pm^3)
        if ~isnan(salt_analysis(i).cat_radius)
            salt_analysis(i).cat_charge_density = salt_analysis(i).cat_charge / ...
                (salt_analysis(i).cat_radius^3);
        else
            salt_analysis(i).cat_charge_density = NaN;
        end
        
        if ~isnan(salt_analysis(i).an_radius)
            salt_analysis(i).an_charge_density = salt_analysis(i).an_charge / ...
                (salt_analysis(i).an_radius^3);
        else
            salt_analysis(i).an_charge_density = NaN;
        end
        
        % Mean charge density
        if ~isnan(salt_analysis(i).cat_charge_density) && ...
           ~isnan(salt_analysis(i).an_charge_density)
            salt_analysis(i).mean_charge_density = ...
                (salt_analysis(i).cat_charge_density + ...
                 salt_analysis(i).an_charge_density) / 2;
        else
            salt_analysis(i).mean_charge_density = NaN;
        end
        
        % Ionic strength factor (per molality)
        % I = 0.5 * sum(c_i * z_i^2) = 0.5 * m * [nu_+ * z_+^2 + nu_- * z_-^2]
        salt_analysis(i).ionic_strength_factor = 0.5 * ...
            (salt_analysis(i).cat_stoich * salt_analysis(i).cat_charge^2 + ...
             salt_analysis(i).an_stoich * salt_analysis(i).an_charge^2);
        
    else
        % Properties not in database - mark as unknown
        salt_analysis(i).cation = 'Unknown';
        salt_analysis(i).anion = 'Unknown';
        salt_analysis(i).category = 'unknown';
        salt_analysis(i).delta_H_soln = NaN;
        salt_analysis(i).cat_charge = NaN;
        salt_analysis(i).an_charge = NaN;
        salt_analysis(i).cat_radius = NaN;
        salt_analysis(i).an_radius = NaN;
        salt_analysis(i).cat_charge_density = NaN;
        salt_analysis(i).an_charge_density = NaN;
        salt_analysis(i).mean_charge_density = NaN;
        salt_analysis(i).ionic_strength_factor = NaN;
        salt_analysis(i).cat_type = 'u';
        salt_analysis(i).an_type = 'u';
        salt_analysis(i).kc_pair_type = 'u-u';
    end
    
    % Get valid RH range for this salt
    if isKey(salt_RH_ranges, salt_name)
        salt_analysis(i).valid_RH_range = salt_RH_ranges(salt_name);
    else
        % If not in database, use data-derived range
        salt_analysis(i).valid_RH_range = [min(salt_data.RH_Water_Activity), ...
                                           max(salt_data.RH_Water_Activity)];
    end
    
    % Filter data to only include points within valid RH range
    RH_valid_mask = salt_data.RH_Water_Activity >= salt_analysis(i).valid_RH_range(1) & ...
                    salt_data.RH_Water_Activity <= salt_analysis(i).valid_RH_range(2);
    
    salt_data_valid = salt_data(RH_valid_mask, :);
    
    % Statistical properties from valid data only
    if height(salt_data_valid) > 0
        salt_analysis(i).molality_range = [min(salt_data_valid.Molality_mol_per_kg), ...
                                           max(salt_data_valid.Molality_mol_per_kg)];
        salt_analysis(i).RH_range = [min(salt_data_valid.RH_Water_Activity), ...
                                     max(salt_data_valid.RH_Water_Activity)];
        salt_analysis(i).gamma_range = [min(salt_data_valid.Activity_Coefficient_Water), ...
                                        max(salt_data_valid.Activity_Coefficient_Water)];
        
        % Mean values (from valid data only)
        salt_analysis(i).mean_molality = mean(salt_data_valid.Molality_mol_per_kg);
        salt_analysis(i).mean_RH = mean(salt_data_valid.RH_Water_Activity);
        salt_analysis(i).mean_gamma = mean(salt_data_valid.Activity_Coefficient_Water);
        
        % Deviation from ideality
        salt_analysis(i).mean_deviation_from_ideal = salt_analysis(i).mean_gamma - 1.0;
        
        % Store full data arrays (filtered to valid RH range only)
        salt_analysis(i).molality = salt_data_valid.Molality_mol_per_kg;
        salt_analysis(i).RH = salt_data_valid.RH_Water_Activity;
        salt_analysis(i).gamma = salt_data_valid.Activity_Coefficient_Water;
        salt_analysis(i).x_water = salt_data_valid.Mole_Fraction_Water;
        
        % Calculate logarithmic quantities
        salt_analysis(i).ln_gamma = log(salt_analysis(i).gamma);
        salt_analysis(i).ln_aw = log(salt_analysis(i).RH);  % aw = RH
        salt_analysis(i).mean_ln_gamma = mean(salt_analysis(i).ln_gamma);
        salt_analysis(i).mean_ln_aw = mean(salt_analysis(i).ln_aw);
    else
        % No valid data - set to empty/NaN
        salt_analysis(i).molality_range = [NaN, NaN];
        salt_analysis(i).RH_range = [NaN, NaN];
        salt_analysis(i).gamma_range = [NaN, NaN];
        salt_analysis(i).mean_molality = NaN;
        salt_analysis(i).mean_RH = NaN;
        salt_analysis(i).mean_gamma = NaN;
        salt_analysis(i).mean_deviation_from_ideal = NaN;
        salt_analysis(i).molality = [];
        salt_analysis(i).RH = [];
        salt_analysis(i).gamma = [];
        salt_analysis(i).x_water = [];
        salt_analysis(i).ln_gamma = [];
        salt_analysis(i).ln_aw = [];
        salt_analysis(i).mean_ln_gamma = NaN;
        salt_analysis(i).mean_ln_aw = NaN;
    end
end

write_output('Property calculations complete\n\n');

%% Define custom color palette for better distinguishability
% Base color palette with highly distinguishable colors
% Colors chosen to maximize perceptual distance
distinguishable_base_colors = [
    0, 0.4470, 0.7410;      % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.6350, 0.0780, 0.1840;  % Red
    0.3010, 0.7450, 0.9330;  % Cyan
    0.7500, 0.0000, 0.7500;  % Magenta
    0.0000, 0.5000, 0.0000;  % Dark green
    0.8700, 0.4900, 0.0000;  % Dark orange
    0.5000, 0.5000, 0.0000;  % Olive
    0.0000, 0.7500, 0.7500;  % Teal
    0.7500, 0.7500, 0.0000;  % Yellow-green
    0.7500, 0.0000, 0.0000;  % Dark red
    0.0000, 0.0000, 0.7500;  % Dark blue
    0.5000, 0.0000, 0.5000;  % Dark magenta
    0.0000, 0.5000, 0.5000;  % Dark cyan
    0.5500, 0.5500, 0.5500;  % Gray
];

% Function to get distinguishable colors (function defined at end of file)
get_distinguishable_colors = @(n) get_colors_func(n, distinguishable_base_colors);

% Function to get viridis-like color from value between 0 and 1
get_viridis_color = @(t) get_viridis_color_func(t);

%% ANALYSIS: Faceted Plots by Ion Family
write_output('=== Generating Faceted Plots by Ion Family ===\n');

% Prepare data for faceted plotting
% Collect all data points with their ion identities and properties
plot_data = struct();
data_idx = 1;

for i = 1:length(salt_analysis)
    if ~strcmp(salt_analysis(i).category, 'unknown') && ...
       ~isempty(salt_analysis(i).molality) && ...
       ~isnan(salt_analysis(i).cat_radius) && ...
       ~isnan(salt_analysis(i).an_radius)
        
        n_points = length(salt_analysis(i).molality);
        
        for j = 1:n_points
            plot_data(data_idx).salt_name = salt_analysis(i).name;
            plot_data(data_idx).cation = salt_analysis(i).cation;
            plot_data(data_idx).anion = salt_analysis(i).anion;
            plot_data(data_idx).cat_radius = salt_analysis(i).cat_radius;
            plot_data(data_idx).an_radius = salt_analysis(i).an_radius;
            plot_data(data_idx).molality = salt_analysis(i).molality(j);
            plot_data(data_idx).RH = salt_analysis(i).RH(j);
            plot_data(data_idx).gamma = salt_analysis(i).gamma(j);
            plot_data(data_idx).ln_aw = salt_analysis(i).ln_aw(j);
            plot_data(data_idx).ln_gamma = salt_analysis(i).ln_gamma(j);
            plot_data(data_idx).x_water = salt_analysis(i).x_water(j);
            
            % Calculate ionic strength for this point
            I = salt_analysis(i).molality(j) * salt_analysis(i).ionic_strength_factor;
            plot_data(data_idx).ionic_strength = I;
            
            data_idx = data_idx + 1;
            end
        end
    end
    
write_output('Prepared %d data points for faceted plotting\n', length(plot_data));

% Get unique cations and anions
unique_cations = unique({plot_data.cation});
unique_anions = unique({plot_data.anion});

write_output('Found %d unique cations and %d unique anions\n', ...
             length(unique_cations), length(unique_anions));

%% PLOT 1: ln(aw) vs Ionic Strength

% Version A: Faceted by Anion, Colored by Cation Radius
figure('Position', [100, 100, 1800, 1200]);
n_anions = length(unique_anions);
n_cols = min(4, n_anions);
n_rows = ceil(n_anions / n_cols);

for a_idx = 1:n_anions
    subplot(n_rows, n_cols, a_idx);
    hold on; grid on; box on;
    
    anion = unique_anions{a_idx};
    anion_mask = strcmp({plot_data.anion}, anion);
    anion_data = plot_data(anion_mask);
    
    if ~isempty(anion_data)
        % Get unique salts for this anion
        salts_in_facet = unique({anion_data.salt_name});
        
        % Get cation radii for color mapping
        all_cat_radii = [anion_data.cat_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_radius = 329.0;   % Smallest hydrated cation
        global_max_radius = 438.0;  % Largest hydrated cation
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({anion_data.salt_name}, salt);
            salt_data = anion_data(salt_mask);
            
            % Get cation radius for this salt
            cat_radius = salt_data(1).cat_radius;
            
            % Map radius to viridis color (blue to green to yellow) - NO normalization
            t = (cat_radius - global_min_radius) / (global_max_radius - global_min_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by ionic strength for smooth lines
            I_vals = [salt_data.ionic_strength];
            ln_aw_vals = [salt_data.ln_aw];
            [I_sorted, sort_idx] = sort(I_vals);
            ln_aw_sorted = ln_aw_vals(sort_idx);
            
            plot(I_sorted, ln_aw_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at end of line
            text(I_sorted(end), ln_aw_sorted(end), ['  ' strrep(salt, '_', '\_')], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold');
        end
    end
    
    xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
    ylabel('ln(a_w)', 'FontWeight', 'bold');
    title(sprintf('%s^-', strrep(anion, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(a_w) vs Ionic Strength - Faceted by Anion (colored by cation hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
caxis(cb_ax, [329.0, 438.0]);  % Set to actual hydrated radius range in pm
ylabel(cb, 'Cation Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_I_by_anion_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_I_by_anion_hydrated.fig'));

% Version B: Faceted by Cation, Colored by Anion Radius
figure('Position', [100, 100, 1800, 1200]);
n_cations = length(unique_cations);
n_cols = min(4, n_cations);
n_rows = ceil(n_cations / n_cols);

for c_idx = 1:n_cations
    subplot(n_rows, n_cols, c_idx);
hold on; grid on; box on;

    cation = unique_cations{c_idx};
    cation_mask = strcmp({plot_data.cation}, cation);
    cation_data = plot_data(cation_mask);
    
    if ~isempty(cation_data)
        % Get unique salts for this cation
        salts_in_facet = unique({cation_data.salt_name});
        
        % Get anion radii for color mapping
        all_an_radii = [cation_data.an_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_an_radius = 300.0;   % Smallest hydrated anion
        global_max_an_radius = 379.0;  % Largest hydrated anion
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({cation_data.salt_name}, salt);
            salt_data = cation_data(salt_mask);
            
            % Get anion radius for this salt
            an_radius = salt_data(1).an_radius;
            
            % Map radius to viridis color (blue to green to yellow) - NO normalization
            t = (an_radius - global_min_an_radius) / (global_max_an_radius - global_min_an_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by ionic strength for smooth lines
            I_vals = [salt_data.ionic_strength];
            ln_aw_vals = [salt_data.ln_aw];
            [I_sorted, sort_idx] = sort(I_vals);
            ln_aw_sorted = ln_aw_vals(sort_idx);
            
            plot(I_sorted, ln_aw_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at end of line
            text(I_sorted(end), ln_aw_sorted(end), ['  ' strrep(salt, '_', '\_')], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold');
        end
    end
    
    xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
    ylabel('ln(a_w)', 'FontWeight', 'bold');
    title(sprintf('%s^+', strrep(cation, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(a_w) vs Ionic Strength - Faceted by Cation (colored by anion hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
caxis(cb_ax, [300.0, 379.0]);  % Set to actual hydrated anion radius range in pm
ylabel(cb, 'Anion Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_I_by_cation_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_I_by_cation_hydrated.fig'));

write_output('Generated ln(aw) vs Ionic Strength faceted plots\n');

%% PLOT 2: ln(aw) vs Mole Fraction Water

% Version A: Faceted by Anion, Colored by Cation Radius
figure('Position', [100, 100, 1800, 1200]);

for a_idx = 1:n_anions
    subplot(n_rows, n_cols, a_idx);
hold on; grid on; box on;
    
    anion = unique_anions{a_idx};
    anion_mask = strcmp({plot_data.anion}, anion);
    anion_data = plot_data(anion_mask);
    
    if ~isempty(anion_data)
        % Get unique salts for this anion
        salts_in_facet = unique({anion_data.salt_name});
        
        % Get cation radii for color mapping
        all_cat_radii = [anion_data.cat_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_radius = 329.0;   % Smallest hydrated cation
        global_max_radius = 438.0;  % Largest hydrated cation
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({anion_data.salt_name}, salt);
            salt_data = anion_data(salt_mask);
            
            % Get cation radius for this salt
            cat_radius = salt_data(1).cat_radius;
            
            % Map radius to color - NO normalization
            t = (cat_radius - global_min_radius) / (global_max_radius - global_min_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by x_water for smooth lines
            x_water_vals = [salt_data.x_water];
            ln_aw_vals = [salt_data.ln_aw];
            [x_sorted, sort_idx] = sort(x_water_vals);
            ln_aw_sorted = ln_aw_vals(sort_idx);
            
            plot(x_sorted, ln_aw_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(x_sorted(1), ln_aw_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
        end
    end
    
    xlabel('Mole Fraction Water', 'FontWeight', 'bold');
    ylabel('ln(a_w)', 'FontWeight', 'bold');
    title(sprintf('%s^-', strrep(anion, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(a_w) vs Mole Fraction Water - Faceted by Anion (colored by cation hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
caxis(cb_ax, [329.0, 438.0]);  % Set to actual hydrated radius range in pm
ylabel(cb, 'Cation Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_xwater_by_anion_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_xwater_by_anion_hydrated.fig'));

% Version B: Faceted by Cation, Colored by Anion Radius
figure('Position', [100, 100, 1800, 1200]);

for c_idx = 1:n_cations
    subplot(n_rows, n_cols, c_idx);
hold on; grid on; box on;

    cation = unique_cations{c_idx};
    cation_mask = strcmp({plot_data.cation}, cation);
    cation_data = plot_data(cation_mask);
    
    if ~isempty(cation_data)
        % Get unique salts for this cation
        salts_in_facet = unique({cation_data.salt_name});
        
        % Get anion radii for color mapping
        all_an_radii = [cation_data.an_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_an_radius = 300.0;   % Smallest hydrated anion
        global_max_an_radius = 379.0;  % Largest hydrated anion
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({cation_data.salt_name}, salt);
            salt_data = cation_data(salt_mask);
            
            % Get anion radius for this salt
            an_radius = salt_data(1).an_radius;
            
            % Map radius to color - NO normalization
            t = (an_radius - global_min_an_radius) / (global_max_an_radius - global_min_an_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by x_water for smooth lines
            x_water_vals = [salt_data.x_water];
            ln_aw_vals = [salt_data.ln_aw];
            [x_sorted, sort_idx] = sort(x_water_vals);
            ln_aw_sorted = ln_aw_vals(sort_idx);
            
            plot(x_sorted, ln_aw_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(x_sorted(1), ln_aw_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
    end
end

    xlabel('Mole Fraction Water', 'FontWeight', 'bold');
    ylabel('ln(a_w)', 'FontWeight', 'bold');
    title(sprintf('%s^+', strrep(cation, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(a_w) vs Mole Fraction Water - Faceted by Cation (colored by anion hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
caxis(cb_ax, [300.0, 379.0]);  % Set to actual hydrated anion radius range in pm
ylabel(cb, 'Anion Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_xwater_by_cation_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_aw_vs_xwater_by_cation_hydrated.fig'));

write_output('Generated ln(aw) vs Mole Fraction Water faceted plots\n');

%% PLOT 3: ln(gamma) vs Ionic Strength

% Version A: Faceted by Anion, Colored by Cation Radius
figure('Position', [100, 100, 1800, 1200]);

for a_idx = 1:n_anions
    subplot(n_rows, n_cols, a_idx);
hold on; grid on; box on;

    anion = unique_anions{a_idx};
    anion_mask = strcmp({plot_data.anion}, anion);
    anion_data = plot_data(anion_mask);
    
    if ~isempty(anion_data)
        % Get unique salts for this anion
        salts_in_facet = unique({anion_data.salt_name});
        
        % Get cation radii for color mapping
        all_cat_radii = [anion_data.cat_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_radius = 329.0;   % Smallest hydrated cation
        global_max_radius = 438.0;  % Largest hydrated cation
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({anion_data.salt_name}, salt);
            salt_data = anion_data(salt_mask);
            
            % Get cation radius for this salt
            cat_radius = salt_data(1).cat_radius;
            
            % Map radius to color - NO normalization
            t = (cat_radius - global_min_radius) / (global_max_radius - global_min_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by ionic strength for smooth lines
            I_vals = [salt_data.ionic_strength];
            ln_gamma_vals = [salt_data.ln_gamma];
            [I_sorted, sort_idx] = sort(I_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(I_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at end of line
            text(I_sorted(end), ln_gamma_sorted(end), ['  ' strrep(salt, '_', '\_')], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold');
        end
        
        % Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^-', strrep(anion, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs Ionic Strength - Faceted by Anion (colored by cation hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
caxis(cb_ax, [329.0, 438.0]);  % Set to actual hydrated radius range in pm
ylabel(cb, 'Cation Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_I_by_anion_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_I_by_anion_hydrated.fig'));

% Version B: Faceted by Cation, Colored by Anion Radius
figure('Position', [100, 100, 1800, 1200]);

for c_idx = 1:n_cations
    subplot(n_rows, n_cols, c_idx);
hold on; grid on; box on;

    cation = unique_cations{c_idx};
    cation_mask = strcmp({plot_data.cation}, cation);
    cation_data = plot_data(cation_mask);
    
    if ~isempty(cation_data)
        % Get unique salts for this cation
        salts_in_facet = unique({cation_data.salt_name});
        
        % Get anion radii for color mapping
        all_an_radii = [cation_data.an_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_an_radius = 300.0;   % Smallest hydrated anion
        global_max_an_radius = 379.0;  % Largest hydrated anion
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({cation_data.salt_name}, salt);
            salt_data = cation_data(salt_mask);
            
            % Get anion radius for this salt
            an_radius = salt_data(1).an_radius;
            
            % Map radius to color - NO normalization
            t = (an_radius - global_min_an_radius) / (global_max_an_radius - global_min_an_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by ionic strength for smooth lines
            I_vals = [salt_data.ionic_strength];
            ln_gamma_vals = [salt_data.ln_gamma];
            [I_sorted, sort_idx] = sort(I_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(I_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at end of line
            text(I_sorted(end), ln_gamma_sorted(end), ['  ' strrep(salt, '_', '\_')], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold');
        end
        
        % Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^+', strrep(cation, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs Ionic Strength - Faceted by Cation (colored by anion hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
    caxis(cb_ax, [300.0, 379.0]);  % Set to actual hydrated anion radius range in pm
    ylabel(cb, 'Anion Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_I_by_cation_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_I_by_cation_hydrated.fig'));

write_output('Generated ln(gamma) vs Ionic Strength faceted plots\n');

%% PLOT 4: ln(gamma) vs Mole Fraction Water

% Version A: Faceted by Anion, Colored by Cation Radius
figure('Position', [100, 100, 1800, 1200]);

for a_idx = 1:n_anions
    subplot(n_rows, n_cols, a_idx);
    hold on; grid on; box on;
    
    anion = unique_anions{a_idx};
    anion_mask = strcmp({plot_data.anion}, anion);
    anion_data = plot_data(anion_mask);
    
    if ~isempty(anion_data)
        % Get unique salts for this anion
        salts_in_facet = unique({anion_data.salt_name});
        
        % Get cation radii for color mapping
        all_cat_radii = [anion_data.cat_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_radius = 329.0;   % Smallest hydrated cation
        global_max_radius = 438.0;  % Largest hydrated cation
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({anion_data.salt_name}, salt);
            salt_data = anion_data(salt_mask);
            
            % Get cation radius for this salt
            cat_radius = salt_data(1).cat_radius;
            
            % Map radius to color - NO normalization
            t = (cat_radius - global_min_radius) / (global_max_radius - global_min_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by x_water for smooth lines
            x_water_vals = [salt_data.x_water];
            ln_gamma_vals = [salt_data.ln_gamma];
            [x_sorted, sort_idx] = sort(x_water_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(x_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(x_sorted(1), ln_gamma_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
end

% Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Mole Fraction Water', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^-', strrep(anion, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs Mole Fraction Water - Faceted by Anion (colored by cation hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
    caxis(cb_ax, [329.0, 438.0]);  % Set to actual hydrated radius range in pm
    ylabel(cb, 'Cation Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_xwater_by_anion_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_xwater_by_anion_hydrated.fig'));

% Version B: Faceted by Cation, Colored by Anion Radius
figure('Position', [100, 100, 1800, 1200]);

for c_idx = 1:n_cations
    subplot(n_rows, n_cols, c_idx);
    hold on; grid on; box on;
    
    cation = unique_cations{c_idx};
    cation_mask = strcmp({plot_data.cation}, cation);
    cation_data = plot_data(cation_mask);
    
    if ~isempty(cation_data)
        % Get unique salts for this cation
        salts_in_facet = unique({cation_data.salt_name});
        
        % Get anion radii for color mapping
        all_an_radii = [cation_data.an_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_an_radius = 300.0;   % Smallest hydrated anion
        global_max_an_radius = 379.0;  % Largest hydrated anion
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({cation_data.salt_name}, salt);
            salt_data = cation_data(salt_mask);
            
            % Get anion radius for this salt
            an_radius = salt_data(1).an_radius;
            
            % Map radius to color - NO normalization
            t = (an_radius - global_min_an_radius) / (global_max_an_radius - global_min_an_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by x_water for smooth lines
            x_water_vals = [salt_data.x_water];
            ln_gamma_vals = [salt_data.ln_gamma];
            [x_sorted, sort_idx] = sort(x_water_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(x_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(x_sorted(1), ln_gamma_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
        end
        
        % Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Mole Fraction Water', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^+', strrep(cation, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs Mole Fraction Water - Faceted by Cation (colored by anion hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
    caxis(cb_ax, [300.0, 379.0]);  % Set to actual hydrated anion radius range in pm
    ylabel(cb, 'Anion Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_xwater_by_cation_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_xwater_by_cation_hydrated.fig'));

write_output('Generated ln(gamma) vs Mole Fraction Water faceted plots\n');

%% PLOT 5: ln(gamma) vs RH

% Version A: Faceted by Anion, Colored by Cation Radius
figure('Position', [100, 100, 1800, 1200]);

for a_idx = 1:n_anions
    subplot(n_rows, n_cols, a_idx);
    hold on; grid on; box on;
    
    anion = unique_anions{a_idx};
    anion_mask = strcmp({plot_data.anion}, anion);
    anion_data = plot_data(anion_mask);
    
    if ~isempty(anion_data)
        % Get unique salts for this anion
        salts_in_facet = unique({anion_data.salt_name});
        
        % Get cation radii for color mapping
        all_cat_radii = [anion_data.cat_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_radius = 329.0;   % Smallest hydrated cation
        global_max_radius = 438.0;  % Largest hydrated cation
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({anion_data.salt_name}, salt);
            salt_data = anion_data(salt_mask);
            
            % Get cation radius for this salt
            cat_radius = salt_data(1).cat_radius;
            
            % Map radius to color - NO normalization
            t = (cat_radius - global_min_radius) / (global_max_radius - global_min_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by RH for smooth lines
            RH_vals = [salt_data.RH];
            ln_gamma_vals = [salt_data.ln_gamma];
            [RH_sorted, sort_idx] = sort(RH_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(RH_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(RH_sorted(1), ln_gamma_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
        end
        
        % Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('RH (Water Activity)', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^-', strrep(anion, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs RH - Faceted by Anion (colored by cation hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
    caxis(cb_ax, [329.0, 438.0]);  % Set to actual hydrated radius range in pm
    ylabel(cb, 'Cation Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_RH_by_anion_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_RH_by_anion_hydrated.fig'));

% Version B: Faceted by Cation, Colored by Anion Radius
figure('Position', [100, 100, 1800, 1200]);

for c_idx = 1:n_cations
    subplot(n_rows, n_cols, c_idx);
    hold on; grid on; box on;
    
    cation = unique_cations{c_idx};
    cation_mask = strcmp({plot_data.cation}, cation);
    cation_data = plot_data(cation_mask);
    
    if ~isempty(cation_data)
        % Get unique salts for this cation
        salts_in_facet = unique({cation_data.salt_name});
        
        % Get anion radii for color mapping
        all_an_radii = [cation_data.an_radius];
        
        % Define global radius range for consistent coloring (in pm)
        global_min_an_radius = 300.0;   % Smallest hydrated anion
        global_max_an_radius = 379.0;  % Largest hydrated anion
        
        % Plot each salt as a line
        for s = 1:length(salts_in_facet)
            salt = salts_in_facet{s};
            salt_mask = strcmp({cation_data.salt_name}, salt);
            salt_data = cation_data(salt_mask);
            
            % Get anion radius for this salt
            an_radius = salt_data(1).an_radius;
            
            % Map radius to color - NO normalization
            t = (an_radius - global_min_an_radius) / (global_max_an_radius - global_min_an_radius);
            t = max(0, min(1, t));  % Clamp to [0, 1]
            line_color = get_viridis_color(t);
            
            % Sort by RH for smooth lines
            RH_vals = [salt_data.RH];
            ln_gamma_vals = [salt_data.ln_gamma];
            [RH_sorted, sort_idx] = sort(RH_vals);
            ln_gamma_sorted = ln_gamma_vals(sort_idx);
            
            plot(RH_sorted, ln_gamma_sorted, '-', 'LineWidth', 1.5, 'Color', line_color);
            
            % Add text label at start of line
            text(RH_sorted(1), ln_gamma_sorted(1), [strrep(salt, '_', '\_') '  '], ...
                 'FontSize', 7, 'Color', line_color, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'right');
        end
        
        % Add ideal line
        plot(xlim, [0, 0], 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('RH (Water Activity)', 'FontWeight', 'bold');
    ylabel('ln(\gamma_w)', 'FontWeight', 'bold');
    title(sprintf('%s^+', strrep(cation, '_', '\_')), 'FontWeight', 'bold');
    set(gca, 'FontSize', 9);
end

sgtitle('ln(\gamma_w) vs RH - Faceted by Cation (colored by anion hydrated radius)', ...
        'FontSize', 14, 'FontWeight', 'bold');

% Add colorbar for radius
cb_ax = axes('Position', [0.92, 0.15, 0.02, 0.7], 'Visible', 'off');
colormap(cb_ax, viridis_colormap());
cb = colorbar(cb_ax, 'Location', 'east');
    caxis(cb_ax, [300.0, 379.0]);  % Set to actual hydrated anion radius range in pm
    ylabel(cb, 'Anion Hydrated Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);

saveas(gcf, fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_RH_by_cation_hydrated.png'));
savefig(fullfile(filepath, '../../..', 'figures', 'faceted_ln_gamma_vs_RH_by_cation_hydrated.fig'));

write_output('Generated ln(gamma) vs RH faceted plots\n');
write_output('All faceted plots completed! Generated 10 figure sets.\n\n');

%% Helper function for distinguishable colors (defined at end)
function colors = get_colors_func(n, base_colors)
    if n <= size(base_colors, 1)
        colors = base_colors(1:n, :);
    else
        % Use base colors and interpolate additional ones
        colors = base_colors;
        % Add interpolated colors for larger sets
        for i = size(base_colors, 1)+1:n
            idx = mod(i-1, size(base_colors, 1)) + 1;
            next_idx = mod(i, size(base_colors, 1)) + 1;
            t = 0.5;  % Midpoint interpolation
            colors(i, :) = colors(idx, :) * (1-t) + colors(next_idx, :) * t;
        end
    end
end

%% Helper function for viridis-like colormap
function cmap = viridis_colormap()
    % Viridis-like colormap (blue -> green -> yellow)
    % Define key colors
    blue = [0.267, 0.004, 0.329];   % Dark blue/purple
    teal = [0.127, 0.566, 0.550];   % Teal/green
    yellow = [0.993, 0.906, 0.144]; % Bright yellow
    
    % Create 256-color map interpolating through these points
    n = 256;
    half = round(n/2);
    
    % First half: blue to teal
    cmap1 = zeros(half, 3);
    for i = 1:half
        t = (i-1) / (half-1);
        cmap1(i, :) = blue * (1-t) + teal * t;
    end
    
    % Second half: teal to yellow
    cmap2 = zeros(n-half, 3);
    for i = 1:n-half
        t = (i-1) / (n-half-1);
        cmap2(i, :) = teal * (1-t) + yellow * t;
    end
    
    cmap = [cmap1; cmap2];
end

%% Helper function to get viridis color for a value between 0 and 1
function color = get_viridis_color_func(t)
    % Viridis-like color interpolation (blue -> green -> yellow)
    % Define key colors
    blue = [0.267, 0.004, 0.329];   % Dark blue/purple
    teal = [0.127, 0.566, 0.550];   % Teal/green
    yellow = [0.993, 0.906, 0.144]; % Bright yellow
    
    % Clamp t to [0, 1]
    t = max(0, min(1, t));
    
    % Interpolate through blue -> teal -> yellow
    if t < 0.5
        % First half: blue to teal
        t_scaled = t * 2;
        color = blue * (1-t_scaled) + teal * t_scaled;
    else
        % Second half: teal to yellow
        t_scaled = (t - 0.5) * 2;
        color = teal * (1-t_scaled) + yellow * t_scaled;
    end
end
