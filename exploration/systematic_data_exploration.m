close all
clear
clc

% =========================================================================
% SYSTEMATIC EXPLORATION OF WATER ACTIVITY DATA
% =========================================================================
% This script performs comprehensive analysis of molality, RH, and activity
% coefficient data considering:
% - Salt properties (MW, dissolution enthalpy)
% - Cation and anion identities
% - Ionic radii and charge densities
% - Chemical families (halides, sulfates, nitrates, etc.)
% - Thermodynamic patterns
% =========================================================================

%% Setup
fprintf('\n========================================\n');
fprintf('SYSTEMATIC WATER ACTIVITY DATA EXPLORATION\n');
fprintf('========================================\n\n');

% Add paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'data'));
addpath(fullfile(filepath, '..', 'util'));

% Open output file for writing analysis results
output_filename = fullfile(filepath, 'systematic_exploration_results.txt');
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
write_output('SYSTEMATIC WATER ACTIVITY DATA EXPLORATION\n');
write_output('========================================\n');
write_output('Date: %s\n', datestr(now));
write_output('========================================\n\n');

% Constants
MWw = 18.015; % g/mol
R = 8.314; % J/(mol·K)
T_kelvin = 298.15; % 25°C

%% Load data
write_output('Loading data...\n');
data_file = fullfile(filepath, '..', 'data', 'water_activity_all_salts_combined.csv');
data_table = readtable(data_file);

% Extract unique salts
unique_salts = unique(data_table.Salt);
n_salts = length(unique_salts);
write_output('Loaded data for %d salts\n', n_salts);
write_output('Total data points: %d\n\n', height(data_table));

%% Define comprehensive salt properties database
write_output('Building salt properties database...\n');

% Ionic radii (pm) - Shannon (1976), 6-coordinate unless noted
% Cations
ionic_radii.Li = 76;    % Li+
ionic_radii.Na = 102;   % Na+
ionic_radii.K = 138;    % K+
ionic_radii.Rb = 152;   % Rb+
ionic_radii.Cs = 167;   % Cs+
ionic_radii.Mg = 72;    % Mg2+
ionic_radii.Ca = 100;   % Ca2+
ionic_radii.Sr = 118;   % Sr2+
ionic_radii.Ba = 135;   % Ba2+
ionic_radii.Zn = 74;    % Zn2+
ionic_radii.Cu = 73;    % Cu2+ (high spin)
ionic_radii.Ni = 69;    % Ni2+
ionic_radii.Mn = 83;    % Mn2+ (high spin)
ionic_radii.Ag = 115;   % Ag+
ionic_radii.NH4 = 148;  % NH4+ (effective)
ionic_radii.H = 140;    % H3O+ (effective)

% Anions
ionic_radii.F = 133;    % F-
ionic_radii.Cl = 181;   % Cl-
ionic_radii.Br = 196;   % Br-
ionic_radii.I = 220;    % I-
ionic_radii.OH = 137;   % OH-
ionic_radii.NO3 = 200;  % NO3- (effective for polyatomic)
ionic_radii.SO4 = 244;  % SO42- (effective)
ionic_radii.ClO3 = 195; % ClO3- (effective)
ionic_radii.ClO4 = 240; % ClO4- (effective)

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
        if isfield(ionic_radii, salt_analysis(i).cation)
            salt_analysis(i).cat_radius = ionic_radii.(salt_analysis(i).cation);
        else
            salt_analysis(i).cat_radius = NaN;
        end
        
        if isfield(ionic_radii, salt_analysis(i).anion)
            salt_analysis(i).an_radius = ionic_radii.(salt_analysis(i).anion);
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

%% ANALYSIS 1: Activity Coefficient vs RH by Chemical Family
write_output('=== Analysis 1: Activity Coefficient vs RH by Chemical Family ===\n');

categories = unique({salt_analysis.category});
% Filter out 'unknown' category if present
categories = categories(~strcmp(categories, 'unknown'));
n_categories = length(categories);

% Calculate grid dimensions to fit all categories
n_cols = 3;
n_rows = ceil(n_categories / n_cols);

figure('Position', [100, 100, 1400, n_rows * 300]);

for cat_idx = 1:n_categories
    subplot(n_rows, n_cols, cat_idx);
    hold on; grid on; box on;
    
    cat = categories{cat_idx};
    cat_mask = strcmp({salt_analysis.category}, cat);
    cat_salts = salt_analysis(cat_mask);
    
    colors = get_distinguishable_colors(length(cat_salts));
    
    for s = 1:length(cat_salts)
        % Only plot if we have valid data
        if ~isempty(cat_salts(s).RH) && ~isempty(cat_salts(s).gamma)
            plot(cat_salts(s).RH, cat_salts(s).gamma, '-', ...
                 'LineWidth', 2, 'Color', colors(s, :), ...
                 'DisplayName', strrep(cat_salts(s).name, '_', '\_'));
        else
            % Debug: report if salt has no data
            if strcmp(cat_salts(s).name, 'NH4NO3')
                write_output('Warning: NH4NO3 found but has empty RH or gamma data\n');
            end
        end
    end
    
    % Add ideal line
    plot([0, 1], [1, 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal');
    
    xlabel('RH (Water Activity)', 'FontWeight', 'bold');
    ylabel('\gamma_w (Activity Coefficient)', 'FontWeight', 'bold');
    title(strrep(cat, '_', ' '), 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8);
    xlim([0, 1]);
    ylim([0.4, 1.1]);
    set(gca, 'FontSize', 10);
end

sgtitle('Water Activity Coefficient vs RH by Chemical Family', ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_gamma_vs_RH_by_family.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_gamma_vs_RH_by_family.fig'));

%% ANALYSIS 2: Activity Coefficient vs RH by Cation
write_output('=== Analysis 2: Activity Coefficient vs RH by Cation ===\n');

cations = unique({salt_analysis.cation});
cations = cations(~strcmp(cations, 'Unknown'));

figure('Position', [100, 100, 1600, 1000]);

for cat_idx = 1:min(9, length(cations))
    subplot(3, 3, cat_idx);
    hold on; grid on; box on;
    
    cat = cations{cat_idx};
    cat_mask = strcmp({salt_analysis.cation}, cat);
    cat_salts = salt_analysis(cat_mask);
    
    colors = get_distinguishable_colors(length(cat_salts));
    
    for s = 1:length(cat_salts)
        % Only plot if we have valid data
        if ~isempty(cat_salts(s).RH) && ~isempty(cat_salts(s).gamma)
            plot(cat_salts(s).RH, cat_salts(s).gamma, 'o-', ...
                 'LineWidth', 1.5, 'MarkerSize', 3, 'Color', colors(s, :), ...
                 'DisplayName', strrep(cat_salts(s).name, '_', '\_'));
        end
    end
    
    % Add ideal line
    plot([0, 1], [1, 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Ideal');
    
    xlabel('RH (Water Activity)', 'FontWeight', 'bold');
    ylabel('\gamma_w', 'FontWeight', 'bold');
    title(['Cation: ' strrep(cat, '_', '\_')], 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 7);
    xlim([0, 1]);
    ylim([0.4, 1.1]);
    set(gca, 'FontSize', 9);
end

sgtitle('Water Activity Coefficient vs RH by Cation', ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_gamma_vs_RH_by_cation.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_gamma_vs_RH_by_cation.fig'));

%% ANALYSIS 3: Charge Density Effects
write_output('=== Analysis 3: Charge Density Effects ===\n');

% Filter salts with valid charge density data
valid_cd = ~isnan([salt_analysis.mean_charge_density]);
cd_salts = salt_analysis(valid_cd);

figure('Position', [100, 100, 1400, 500]);

% Plot 1: Mean gamma vs mean charge density
subplot(1, 3, 1);
hold on; grid on; box on;

mean_cd = [cd_salts.mean_charge_density];
mean_gamma = [cd_salts.mean_gamma];
salt_names_cd = {cd_salts.name};

scatter(mean_cd, mean_gamma, 100, 'filled', 'MarkerEdgeColor', 'k');

% Add labels for each point
for i = 1:length(cd_salts)
    text(mean_cd(i), mean_gamma(i), ['  ' strrep(salt_names_cd{i}, '_', '\_')], ...
         'FontSize', 7, 'Interpreter', 'tex');
end

xlabel('Mean Charge Density (e/pm^3)', 'FontWeight', 'bold');
ylabel('Mean \gamma_w', 'FontWeight', 'bold');
title('Activity Coefficient vs Charge Density', 'FontWeight', 'bold');

% Plot 2: Color by enthalpy of solution
subplot(1, 3, 2);
hold on; grid on; box on;

delta_H = [cd_salts.delta_H_soln];
valid_H = ~isnan(delta_H);

scatter(mean_cd(valid_H), mean_gamma(valid_H), 100, delta_H(valid_H), ...
        'filled', 'MarkerEdgeColor', 'k');
colorbar;
colormap(jet);

xlabel('Mean Charge Density (e/pm^3)', 'FontWeight', 'bold');
ylabel('Mean \gamma_w', 'FontWeight', 'bold');
title('Colored by \DeltaH_{soln} (kJ/mol)', 'FontWeight', 'bold');

% Plot 3: Cation charge density vs anion charge density
subplot(1, 3, 3);
hold on; grid on; box on;

cat_cd = [cd_salts.cat_charge_density];
an_cd = [cd_salts.an_charge_density];

scatter(cat_cd, an_cd, 100, mean_gamma, 'filled', 'MarkerEdgeColor', 'k');
colorbar;
ylabel(colorbar, 'Mean \gamma_w', 'FontWeight', 'bold');

xlabel('Cation Charge Density (e/pm^3)', 'FontWeight', 'bold');
ylabel('Anion Charge Density (e/pm^3)', 'FontWeight', 'bold');
title('Cation vs Anion Charge Densities', 'FontWeight', 'bold');

sgtitle('Charge Density Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_charge_density_analysis.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_charge_density_analysis.fig'));

%% ANALYSIS 4: Ionic Radius Trends
write_output('=== Analysis 4: Ionic Radius Trends ===\n');

% --- Part A: Same anion, different cations (Varying Cation Series) ---
halide_series = {
    {'LiCl', 'NaCl', 'KCl', 'RbCl', 'CsCl', 'NH4Cl'}, 'Chlorides (Cl^-)';
    {'LiBr', 'NaBr', 'KBr', 'CsBr'}, 'Bromides (Br^-)';
    {'LiI', 'NaI', 'KI', 'CsI'}, 'Iodides (I^-)';
    {'LiNO3', 'NaNO3', 'KNO3', 'NH4NO3', 'AgNO3'}, 'Nitrates (NO_3^-)';
    {'MgCl2', 'CaCl2', 'SrCl2', 'BaCl2', 'ZnCl2'}, 'Chlorides 2:1 (Cl^-)';
    {'CaBr2', 'SrBr2', 'BaBr2', 'ZnBr2'}, 'Bromides 2:1 (Br^-)';
};

figure('Position', [100, 100, 1600, 900]);
for series_idx = 1:size(halide_series, 1)
    subplot(2, 3, series_idx);
    hold on; grid on; box on;
    
    series_salts = halide_series{series_idx, 1};
    series_name = halide_series{series_idx, 2};
    
    colors = get_distinguishable_colors(length(series_salts));
    
    for s = 1:length(series_salts)
        salt_name = series_salts{s};
        
        % Find in salt_analysis
        idx = find(strcmp({salt_analysis.name}, salt_name));
        if ~isempty(idx) && ~isempty(salt_analysis(idx).molality) && ~isempty(salt_analysis(idx).gamma)
            if ~isnan(salt_analysis(idx).cat_radius)
                plot(salt_analysis(idx).molality, salt_analysis(idx).gamma, 'o-', ...
                     'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(s, :), ...
                     'DisplayName', sprintf('%s (r_{cat}=%d pm)', ...
                            strrep(salt_name, '_', '\_'), ...
                            round(salt_analysis(idx).cat_radius)));
            else
                plot(salt_analysis(idx).molality, salt_analysis(idx).gamma, 'o-', ...
                     'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(s, :), ...
                     'DisplayName', strrep(salt_name, '_', '\_'));
            end
        end
    end
    
    xlabel('Molality (mol/kg)', 'FontWeight', 'bold');
    ylabel('\gamma_w', 'FontWeight', 'bold');
    title(series_name, 'FontWeight', 'bold', 'Interpreter', 'tex');
    legend('Location', 'best', 'FontSize', 7);
    ylim([0.4, 1.1]);
    grid on;
    set(gca, 'FontSize', 10);
end
sgtitle('Activity Coefficient Trends: Cation Size Effects (Constant Anion)', ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_varying_cation.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_varying_cation.fig'));


% --- Part B: Same cation, different anions (Varying Anion Series) ---
cation_series = {
    {'LiCl', 'LiBr', 'LiI', 'LiNO3', 'LiOH', 'LiClO4'}, 'Lithium Salts (Li^+)';
    {'NaCl', 'NaBr', 'NaI', 'NaNO3', 'NaOH', 'NaClO4'}, 'Sodium Salts (Na^+)';
    {'KCl', 'KBr', 'KI', 'KNO3', 'KClO3'}, 'Potassium Salts (K^+)';
    {'NH4Cl', 'NH4NO3', '(NH4)2SO4'}, 'Ammonium Salts (NH_4^+)';
    {'CaCl2', 'CaBr2', 'CaI2', 'Ca(NO3)2'}, 'Calcium Salts (Ca^{2+})';
    {'MgCl2', 'Mg(NO3)2', 'MgSO4'}, 'Magnesium Salts (Mg^{2+})';
    {'ZnCl2', 'ZnBr2', 'ZnI2', 'ZnSO4'}, 'Zinc Salts (Zn^{2+})';
};

% Calculate grid dimensions for cation series (7 series)
n_series = size(cation_series, 1);
n_cols_b = 3;
n_rows_b = ceil(n_series / n_cols_b);
figure('Position', [100, 100, 1600, n_rows_b * 300]);
for series_idx = 1:n_series
    subplot(n_rows_b, n_cols_b, series_idx);
    hold on; grid on; box on;
    
    series_salts = cation_series{series_idx, 1};
    series_name = cation_series{series_idx, 2};
    
    colors = get_distinguishable_colors(length(series_salts));
    
    for s = 1:length(series_salts)
        salt_name = series_salts{s};
        
        % Find in salt_analysis
        idx = find(strcmp({salt_analysis.name}, salt_name));
        if ~isempty(idx) && ~isempty(salt_analysis(idx).molality) && ~isempty(salt_analysis(idx).gamma)
            if ~isnan(salt_analysis(idx).an_radius)
                plot(salt_analysis(idx).molality, salt_analysis(idx).gamma, 'o-', ...
                     'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(s, :), ...
                     'DisplayName', sprintf('%s (r_{an}=%d pm)', ...
                            strrep(salt_name, '_', '\_'), ...
                            round(salt_analysis(idx).an_radius)));
            else
                plot(salt_analysis(idx).molality, salt_analysis(idx).gamma, 'o-', ...
                     'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(s, :), ...
                     'DisplayName', strrep(salt_name, '_', '\_'));
            end
        end
    end
    
    xlabel('Molality (mol/kg)', 'FontWeight', 'bold');
    ylabel('\gamma_w', 'FontWeight', 'bold');
    title(series_name, 'FontWeight', 'bold', 'Interpreter', 'tex');
    legend('Location', 'best', 'FontSize', 7);
    ylim([0.4, 1.1]);
    grid on;
    set(gca, 'FontSize', 10);
end
sgtitle('Activity Coefficient Trends: Anion Size Effects (Constant Cation)', ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_varying_anion.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_varying_anion.fig'));


% --- Part C: Comparative analysis - Correlations & Slopes ---
figure('Position', [100, 100, 1600, 900]); % Increased height for 2 rows

% -------------------------
% Row 1: Mean Gamma Correlations
% -------------------------

% Plot 1: Mean gamma vs cation radius (for halides)
subplot(2, 3, 1);
hold on; grid on; box on;
cat_radii = [];
mean_gammas = [];
salt_labels_cat = {};
anion_types = {};
for i = 1:length(salt_analysis)
    % Only consider simple halides (1:1)
    if strcmp(salt_analysis(i).category, 'endothermic_halide') || ...
       strcmp(salt_analysis(i).category, 'exothermic_halide')
        if ~isnan(salt_analysis(i).cat_radius) && ...
           isfield(salt_analysis(i), 'an_stoich') && ...
           salt_analysis(i).an_stoich == 1 && ...
           salt_analysis(i).cat_stoich == 1
            cat_radii(end+1) = salt_analysis(i).cat_radius;
            mean_gammas(end+1) = salt_analysis(i).mean_gamma;
            salt_labels_cat{end+1} = salt_analysis(i).name;
            anion_types{end+1} = salt_analysis(i).anion;
        end
    end
end
% Color by anion type
unique_anions = unique(anion_types);
anion_colors = get_distinguishable_colors(length(unique_anions));
for a = 1:length(unique_anions)
    anion_mask = strcmp(anion_types, unique_anions{a});
    scatter(cat_radii(anion_mask), mean_gammas(anion_mask), 120, ...
            anion_colors(a, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', [unique_anions{a} '^-']);
end
xlabel('Cation Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Mean \gamma_w', 'FontWeight', 'bold', 'FontSize', 11);
title('Cation Size Effect (1:1 Halides)', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% Plot 2: Mean gamma vs anion radius
subplot(2, 3, 2);
hold on; grid on; box on;
an_radii = [];
mean_gammas_an = [];
salt_labels_an = {};
cation_types = {};
for i = 1:length(salt_analysis)
    if strcmp(salt_analysis(i).category, 'endothermic_halide') || ...
       strcmp(salt_analysis(i).category, 'exothermic_halide')
        if ~isnan(salt_analysis(i).an_radius) && ...
           isfield(salt_analysis(i), 'an_stoich') && ...
           salt_analysis(i).an_stoich == 1 && ...
           salt_analysis(i).cat_stoich == 1
            an_radii(end+1) = salt_analysis(i).an_radius;
            mean_gammas_an(end+1) = salt_analysis(i).mean_gamma;
            salt_labels_an{end+1} = salt_analysis(i).name;
            cation_types{end+1} = salt_analysis(i).cation;
        end
    end
end
% Color by cation type
unique_cations = unique(cation_types);
cation_colors = get_distinguishable_colors(length(unique_cations));
for c = 1:length(unique_cations)
    cation_mask = strcmp(cation_types, unique_cations{c});
    scatter(an_radii(cation_mask), mean_gammas_an(cation_mask), 120, ...
            cation_colors(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', [unique_cations{c} '^+']);
end
xlabel('Anion Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Mean \gamma_w', 'FontWeight', 'bold', 'FontSize', 11);
title('Anion Size Effect (1:1 Halides)', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% Plot 3: Combined - sum of ionic radii
subplot(2, 3, 3);
hold on; grid on; box on;
for i = 1:length(salt_analysis)
    if strcmp(salt_analysis(i).category, 'endothermic_halide') || ...
       strcmp(salt_analysis(i).category, 'exothermic_halide')
        if ~isnan(salt_analysis(i).cat_radius) && ~isnan(salt_analysis(i).an_radius) && ...
           isfield(salt_analysis(i), 'an_stoich') && ...
           salt_analysis(i).an_stoich == 1 && ...
           salt_analysis(i).cat_stoich == 1
            r_sum = salt_analysis(i).cat_radius + salt_analysis(i).an_radius;
            
            % Color by enthalpy
            if ~isnan(salt_analysis(i).delta_H_soln)
                if salt_analysis(i).delta_H_soln < 0
                    marker_color = [0.8, 0.2, 0.2]; % Red for exothermic
                    marker_shape = 'o';
                else
                    marker_color = [0.2, 0.2, 0.8]; % Blue for endothermic
                    marker_shape = 's';
                end
            else
                marker_color = [0.5, 0.5, 0.5];
                marker_shape = 'd';
            end
            
            scatter(r_sum, salt_analysis(i).mean_gamma, 120, marker_color, ...
                    'filled', marker_shape, 'MarkerEdgeColor', 'k');
            text(r_sum, salt_analysis(i).mean_gamma, ...
                 ['  ' strrep(salt_analysis(i).name, '_', '\_')], ...
                 'FontSize', 7);
        end
    end
end
% Add dummy points for legend
scatter(NaN, NaN, 120, [0.8, 0.2, 0.2], 'filled', 'o', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'Exothermic');
scatter(NaN, NaN, 120, [0.2, 0.2, 0.8], 'filled', 's', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'Endothermic');
xlabel('Sum of Ionic Radii (r_{cat} + r_{an}, pm)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Mean \gamma_w', 'FontWeight', 'bold', 'FontSize', 11);
title('Total Ionic Size Effect', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

% -------------------------
% Row 2: Slope (Sensitivity) Trends
% -------------------------

% Plot 4: Slope vs Cation Radius (grouped by Anion Series)
subplot(2, 3, 4);
hold on; grid on; box on;
halide_colors = get_distinguishable_colors(size(halide_series, 1));

for series_idx = 1:size(halide_series, 1)
    series_salts = halide_series{series_idx, 1};
    series_name_full = halide_series{series_idx, 2};
    % Extract simplified name for legend (e.g., "Chlorides")
    series_name_parts = split(series_name_full, '(');
    series_label = strtrim(series_name_parts{1});
    
    radii_vec = [];
    slopes_vec = [];
    
    for s = 1:length(series_salts)
        salt_name = series_salts{s};
        idx = find(strcmp({salt_analysis.name}, salt_name));
        
        if ~isempty(idx) && ~isempty(salt_analysis(idx).molality) && length(salt_analysis(idx).molality) > 1
             if ~isnan(salt_analysis(idx).cat_radius)
                 % Calculate linear slope (dGamma/dMolality approx)
                 p = polyfit(salt_analysis(idx).molality, salt_analysis(idx).gamma, 1);
                 slope_val = p(1);
                 
                 radii_vec(end+1) = salt_analysis(idx).cat_radius;
                 slopes_vec(end+1) = slope_val;
             end
        end
    end
    
    if ~isempty(radii_vec)
        % Sort by radius for cleaner line plotting
        [radii_vec, sort_i] = sort(radii_vec);
        slopes_vec = slopes_vec(sort_i);
        
        plot(radii_vec, slopes_vec, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', halide_colors(series_idx, :), ...
             'MarkerFaceColor', halide_colors(series_idx, :), ...
             'DisplayName', series_label);
    end
end
xlabel('Cation Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Slope (d\gamma_w / dm)', 'FontWeight', 'bold', 'FontSize', 11);
title('Sensitivity: Slope vs Cation Radius', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Plot 5: Slope vs Anion Radius (grouped by Cation Series)
subplot(2, 3, 5);
hold on; grid on; box on;
cation_series_colors = get_distinguishable_colors(size(cation_series, 1));

for series_idx = 1:size(cation_series, 1)
    series_salts = cation_series{series_idx, 1};
    series_name_full = cation_series{series_idx, 2};
    % Extract simplified name for legend
    series_name_parts = split(series_name_full, '(');
    series_label = strtrim(series_name_parts{1});
    
    radii_vec = [];
    slopes_vec = [];
    
    for s = 1:length(series_salts)
        salt_name = series_salts{s};
        idx = find(strcmp({salt_analysis.name}, salt_name));
        
        if ~isempty(idx) && ~isempty(salt_analysis(idx).molality) && length(salt_analysis(idx).molality) > 1
             if ~isnan(salt_analysis(idx).an_radius)
                 % Calculate linear slope (dGamma/dMolality approx)
                 p = polyfit(salt_analysis(idx).molality, salt_analysis(idx).gamma, 1);
                 slope_val = p(1);
                 
                 radii_vec(end+1) = salt_analysis(idx).an_radius;
                 slopes_vec(end+1) = slope_val;
             end
        end
    end
    
    if ~isempty(radii_vec)
        % Sort by radius for cleaner line plotting
        [radii_vec, sort_i] = sort(radii_vec);
        slopes_vec = slopes_vec(sort_i);
        
        plot(radii_vec, slopes_vec, 's-', 'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', cation_series_colors(series_idx, :), ...
             'MarkerFaceColor', cation_series_colors(series_idx, :), ...
             'DisplayName', series_label);
    end
end
xlabel('Anion Radius (pm)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Slope (d\gamma_w / dm)', 'FontWeight', 'bold', 'FontSize', 11);
title('Sensitivity: Slope vs Anion Radius', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 8);
grid on;

sgtitle('Ionic Radius Correlation Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_correlations.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_ionic_radius_correlations.fig'));

%% ANALYSIS 5: Enthalpy of Solution Correlation
write_output('=== Analysis 5: Enthalpy of Solution Correlation ===\n');

% Filter salts with valid enthalpy data
valid_H = ~isnan([salt_analysis.delta_H_soln]);
H_salts = salt_analysis(valid_H);

figure('Position', [100, 100, 1400, 500]);

% Plot 1: Mean gamma vs delta H
subplot(1, 3, 1);
hold on; grid on; box on;

delta_H_vals = [H_salts.delta_H_soln];
mean_gamma_vals = [H_salts.mean_gamma];

% Color by category
cats = {H_salts.category};
unique_cats = unique(cats);
colors = get_distinguishable_colors(length(unique_cats));

for c = 1:length(unique_cats)
    cat_mask = strcmp(cats, unique_cats{c});
    scatter(delta_H_vals(cat_mask), mean_gamma_vals(cat_mask), 100, ...
            colors(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', strrep(unique_cats{c}, '_', ' '));
end

% Divide into endo/exothermic
plot([0, 0], ylim, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
text(20, 0.75, 'Endothermic', 'FontSize', 12, 'FontWeight', 'bold');
text(-60, 0.75, 'Exothermic', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('\DeltaH_{soln} (kJ/mol)', 'FontWeight', 'bold');
ylabel('Mean \gamma_w', 'FontWeight', 'bold');
title('Activity Coefficient vs Dissolution Enthalpy', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 8);

% Plot 2: RH range width vs delta H
subplot(1, 3, 2);
hold on; grid on; box on;

RH_range_width = arrayfun(@(x) diff(x.RH_range), H_salts);

for c = 1:length(unique_cats)
    cat_mask = strcmp(cats, unique_cats{c});
    scatter(delta_H_vals(cat_mask), RH_range_width(cat_mask), 100, ...
            colors(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', strrep(unique_cats{c}, '_', ' '));
end

plot([0, 0], ylim, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

xlabel('\DeltaH_{soln} (kJ/mol)', 'FontWeight', 'bold');
ylabel('RH Range Width', 'FontWeight', 'bold');
title('RH Operating Range vs Dissolution Enthalpy', 'FontWeight', 'bold');

% Plot 3: Deviation from ideality
subplot(1, 3, 3);
hold on; grid on; box on;

dev_from_ideal = [H_salts.mean_deviation_from_ideal];

for c = 1:length(unique_cats)
    cat_mask = strcmp(cats, unique_cats{c});
    scatter(delta_H_vals(cat_mask), dev_from_ideal(cat_mask), 100, ...
            colors(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', strrep(unique_cats{c}, '_', ' '));
end

plot([0, 0], ylim, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
plot(xlim, [0, 0], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlabel('\DeltaH_{soln} (kJ/mol)', 'FontWeight', 'bold');
ylabel('Mean Deviation from Ideal (\gamma_w - 1)', 'FontWeight', 'bold');
title('Non-ideality vs Dissolution Enthalpy', 'FontWeight', 'bold');

sgtitle('Thermodynamic Correlation Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_enthalpy_correlation.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_enthalpy_correlation.fig'));

%% ANALYSIS 6: Ionic Strength Effects
write_output('=== Analysis 6: Ionic Strength Effects ===\n');

% Calculate ionic strength for all salts at multiple molality values
figure('Position', [100, 100, 1400, 500]);

% Plot 1: 1:1 vs 2:1 salts
subplot(1, 3, 1);
hold on; grid on; box on;

% Separate by ionic strength factor (charge type)
IS_factors = unique([salt_analysis.ionic_strength_factor]);
IS_factors = IS_factors(~isnan(IS_factors));

colors = get_distinguishable_colors(length(IS_factors));

for is_idx = 1:length(IS_factors)
    is_val = IS_factors(is_idx);
    is_mask = [salt_analysis.ionic_strength_factor] == is_val;
    is_salts = salt_analysis(is_mask);
    
    % Plot a few representative salts
    n_plot = min(3, length(is_salts));
    for s = 1:n_plot
        ionic_strength = is_salts(s).molality * is_val;
        plot(ionic_strength, is_salts(s).gamma, 'o', ...
             'Color', colors(is_idx, :), 'MarkerSize', 3, ...
             'DisplayName', sprintf('IS factor=%.1f: %s', is_val, ...
                    strrep(is_salts(s).name, '_', '\_')));
    end
end

xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
ylabel('\gamma_w', 'FontWeight', 'bold');
title('Activity Coefficient vs Ionic Strength', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 7);

% Plot 2: Comparison at fixed molality
subplot(1, 3, 2);
hold on; grid on; box on;

target_molality = 1.0; % mol/kg
tol = 0.3;

mean_gamma_at_m = [];
is_factor_vals = [];
salt_labels = {};

for i = 1:length(salt_analysis)
    if ~isnan(salt_analysis(i).ionic_strength_factor)
        % Find data near target molality
        m_mask = abs(salt_analysis(i).molality - target_molality) < tol;
        if sum(m_mask) > 0
            mean_gamma_at_m(end+1) = mean(salt_analysis(i).gamma(m_mask));
            is_factor_vals(end+1) = salt_analysis(i).ionic_strength_factor;
            salt_labels{end+1} = salt_analysis(i).name;
        end
    end
end

scatter(is_factor_vals, mean_gamma_at_m, 100, 'filled', 'MarkerEdgeColor', 'k');

% Add labels
for i = 1:length(salt_labels)
    text(is_factor_vals(i), mean_gamma_at_m(i), ...
         ['  ' strrep(salt_labels{i}, '_', '\_')], ...
         'FontSize', 7);
end

xlabel('Ionic Strength Factor', 'FontWeight', 'bold');
ylabel('\gamma_w (at m \approx 1)', 'FontWeight', 'bold');
title(sprintf('Activity Coefficient at m = %.1f mol/kg', target_molality), ...
      'FontWeight', 'bold');

% Plot 3: Debye-Hückel comparison
subplot(1, 3, 3);
hold on; grid on; box on;

% Select a few 1:1 salts for Debye-Hückel comparison
dh_salts = {'NaCl', 'KCl', 'LiCl'};
colors_dh = get_distinguishable_colors(length(dh_salts));

for s = 1:length(dh_salts)
    idx = find(strcmp({salt_analysis.name}, dh_salts{s}));
    if ~isempty(idx)
        I = salt_analysis(idx).molality * salt_analysis(idx).ionic_strength_factor;
        
        % Extended Debye-Hückel for water activity coefficient (approximate)
        % This is a simplified version; gamma_w decreases with increasing I
        A_DH = 0.509; % Debye-Hückel constant at 25°C (kg^0.5/mol^0.5)
        
        % Plot experimental - only if we have valid data
        if ~isempty(I) && ~isempty(salt_analysis(idx).gamma)
            plot(I, salt_analysis(idx).gamma, 'o-', ...
                 'Color', colors_dh(s, :), 'LineWidth', 2, 'MarkerSize', 4, ...
                 'DisplayName', [strrep(dh_salts{s}, '_', '\_') ' (exp)']);
        end
    end
end

xlabel('Ionic Strength (mol/kg)', 'FontWeight', 'bold');
ylabel('\gamma_w', 'FontWeight', 'bold');
title('Experimental Data', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);

sgtitle('Ionic Strength Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_ionic_strength.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_ionic_strength.fig'));

%% ANALYSIS 7: Molecular Weight Effects
write_output('=== Analysis 7: Molecular Weight Effects ===\n');

figure('Position', [100, 100, 1400, 500]);

% Plot 1: Mean gamma vs MW
subplot(1, 3, 1);
hold on; grid on; box on;

MW_vals = [salt_analysis.MW];
mean_gamma_all = [salt_analysis.mean_gamma];

scatter(MW_vals, mean_gamma_all, 100, 'filled', 'MarkerEdgeColor', 'k');

xlabel('Molecular Weight (g/mol)', 'FontWeight', 'bold');
ylabel('Mean \gamma_w', 'FontWeight', 'bold');
title('Activity Coefficient vs Molecular Weight', 'FontWeight', 'bold');

% Plot 2: Molality range vs MW
subplot(1, 3, 2);
hold on; grid on; box on;

max_molality = arrayfun(@(x) max(x.molality), salt_analysis);

scatter(MW_vals, max_molality, 100, 'filled', 'MarkerEdgeColor', 'k');

xlabel('Molecular Weight (g/mol)', 'FontWeight', 'bold');
ylabel('Maximum Molality (mol/kg)', 'FontWeight', 'bold');
title('Solubility Range vs Molecular Weight', 'FontWeight', 'bold');

% Plot 3: Comparison by chemical family
subplot(1, 3, 3);
hold on; grid on; box on;

cats_all = {salt_analysis.category};
unique_cats_all = unique(cats_all);
colors_mw = get_distinguishable_colors(length(unique_cats_all));

for c = 1:length(unique_cats_all)
    cat_mask = strcmp(cats_all, unique_cats_all{c});
    scatter(MW_vals(cat_mask), mean_gamma_all(cat_mask), 100, ...
            colors_mw(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', strrep(unique_cats_all{c}, '_', ' '));
end

xlabel('Molecular Weight (g/mol)', 'FontWeight', 'bold');
ylabel('Mean \gamma_w', 'FontWeight', 'bold');
title('MW Effect by Chemical Family', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 8);

sgtitle('Molecular Weight Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_molecular_weight.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_molecular_weight.fig'));

%% ANALYSIS 8: Electrolyte Type and Ion Stoichiometry
write_output('=== Analysis 8: Electrolyte Type and Ion Stoichiometry ===\n');

% Calculate total number of ions for each salt
for i = 1:length(salt_analysis)
    % Check if stoichiometry fields exist and are valid
    has_cat_stoich = isfield(salt_analysis(i), 'cat_stoich') && ...
                     ~isempty(salt_analysis(i).cat_stoich) && ...
                     isscalar(salt_analysis(i).cat_stoich) && ...
                     ~isnan(salt_analysis(i).cat_stoich);
    
    has_an_stoich = isfield(salt_analysis(i), 'an_stoich') && ...
                    ~isempty(salt_analysis(i).an_stoich) && ...
                    isscalar(salt_analysis(i).an_stoich) && ...
                    ~isnan(salt_analysis(i).an_stoich);
    
    if has_cat_stoich && has_an_stoich
        salt_analysis(i).n_cations = salt_analysis(i).cat_stoich;
        salt_analysis(i).n_anions = salt_analysis(i).an_stoich;
        salt_analysis(i).total_ions = salt_analysis(i).cat_stoich + salt_analysis(i).an_stoich;
        
        % Classify electrolyte type
        cat_str = sprintf('%d', salt_analysis(i).cat_stoich);
        an_str = sprintf('%d', salt_analysis(i).an_stoich);
        cat_charge_str = sprintf('%d', salt_analysis(i).cat_charge);
        an_charge_str = sprintf('%d', salt_analysis(i).an_charge);
        
        salt_analysis(i).electrolyte_type = sprintf('%s:%s (z_{+}=%s, z_{-}=%s)', ...
            cat_str, an_str, cat_charge_str, an_charge_str);
        
        % Simple classification
        if salt_analysis(i).cat_charge == 1 && salt_analysis(i).an_charge == 1
            salt_analysis(i).electrolyte_class = '1:1';
        elseif salt_analysis(i).cat_charge == 2 && salt_analysis(i).an_charge == 1
            salt_analysis(i).electrolyte_class = '2:1';
        elseif salt_analysis(i).cat_charge == 1 && salt_analysis(i).an_charge == 2
            salt_analysis(i).electrolyte_class = '1:2';
        elseif salt_analysis(i).cat_charge == 2 && salt_analysis(i).an_charge == 2
            salt_analysis(i).electrolyte_class = '2:2';
        else
            salt_analysis(i).electrolyte_class = 'other';
        end
    else
        salt_analysis(i).n_cations = NaN;
        salt_analysis(i).n_anions = NaN;
        salt_analysis(i).total_ions = NaN;
        salt_analysis(i).electrolyte_type = 'unknown';
        salt_analysis(i).electrolyte_class = 'unknown';
    end
end

% Filter salts with valid ion data
valid_ion_mask = ~isnan([salt_analysis.total_ions]);
ion_salts = salt_analysis(valid_ion_mask);

figure('Position', [100, 100, 1600, 1000]);

% Plot 1: Activity coefficient vs total ions at fixed molality
subplot(2, 3, 1);
hold on; grid on; box on;

target_m = 2.0;
tol_m = 0.5;

total_ions_vec = [];
gamma_vec = [];
salt_labels_ion = {};
colors_by_type = [];

electrolyte_classes = unique({ion_salts.electrolyte_class});
class_colors = get_distinguishable_colors(length(electrolyte_classes));

for i = 1:length(ion_salts)
    % Find data near target molality
    m_mask = abs(ion_salts(i).molality - target_m) < tol_m;
    if sum(m_mask) > 5  % Need reasonable number of points
        total_ions_vec(end+1) = ion_salts(i).total_ions;
        gamma_vec(end+1) = mean(ion_salts(i).gamma(m_mask));
        salt_labels_ion{end+1} = ion_salts(i).name;
        
        % Assign color by electrolyte class
        class_idx = find(strcmp(electrolyte_classes, ion_salts(i).electrolyte_class));
        colors_by_type(end+1, :) = class_colors(class_idx, :);
    end
end

for i = 1:length(total_ions_vec)
    scatter(total_ions_vec(i), gamma_vec(i), 150, colors_by_type(i, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    text(total_ions_vec(i), gamma_vec(i), ['  ' strrep(salt_labels_ion{i}, '_', '\_')], ...
         'FontSize', 8);
end

% Add legend for electrolyte classes
for c = 1:length(electrolyte_classes)
    scatter(NaN, NaN, 150, class_colors(c, :), 'filled', 'MarkerEdgeColor', 'k', ...
            'DisplayName', electrolyte_classes{c});
end
legend('Location', 'best', 'FontSize', 10);

xlabel('Total Number of Ions (ν_{+} + ν_{-})', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('\gamma_w (at m ≈ 2 mol/kg)', 'FontWeight', 'bold', 'FontSize', 12);
title('Activity Coefficient vs Total Ion Number', 'FontWeight', 'bold');

% Plot 2: 3D scatter - n_cations vs n_anions vs gamma
subplot(2, 3, 2);
hold on; grid on; box on;

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    n_cat = [class_salts.n_cations];
    n_an = [class_salts.n_anions];
    mean_gamma_class = [class_salts.mean_gamma];
    
    scatter3(n_cat, n_an, mean_gamma_class, 150, class_colors(c, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', electrolyte_classes{c});
end

xlabel('Number of Cations (ν_{+})', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('Number of Anions (ν_{-})', 'FontWeight', 'bold', 'FontSize', 11);
zlabel('Mean \gamma_w', 'FontWeight', 'bold', 'FontSize', 11);
title('Ion Stoichiometry vs Activity Coefficient', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
view(45, 30);
grid on;

% Plot 3: Activity coefficient by electrolyte class (box plot style)
subplot(2, 3, 3);
hold on; grid on; box on;

class_data = {};
class_labels = {};
for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    % Collect all gamma values for this class
    all_gamma = [];
    for s = 1:length(class_salts)
        all_gamma = [all_gamma; class_salts(s).gamma];
    end
    class_data{c} = all_gamma;
    class_labels{c} = electrolyte_classes{c};
end

% Create violin plot style visualization
for c = 1:length(class_data)
    data = class_data{c};
    
    % Calculate statistics
    q25 = prctile(data, 25);
    q50 = prctile(data, 50);
    q75 = prctile(data, 75);
    mu = mean(data);
    
    % Plot box
    x_pos = c;
    box_width = 0.3;
    
    % Draw quartile box
    rectangle('Position', [x_pos-box_width/2, q25, box_width, q75-q25], ...
              'FaceColor', class_colors(c, :), 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Draw median line
    plot([x_pos-box_width/2, x_pos+box_width/2], [q50, q50], 'k-', 'LineWidth', 2);
    
    % Draw mean marker
    scatter(x_pos, mu, 80, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    
    % Draw whiskers (1.5*IQR)
    iqr = q75 - q25;
    lower_whisker = max(min(data), q25 - 1.5*iqr);
    upper_whisker = min(max(data), q75 + 1.5*iqr);
    plot([x_pos, x_pos], [q25, lower_whisker], 'k-', 'LineWidth', 1.5);
    plot([x_pos, x_pos], [q75, upper_whisker], 'k-', 'LineWidth', 1.5);
    plot([x_pos-box_width/4, x_pos+box_width/4], [lower_whisker, lower_whisker], 'k-', 'LineWidth', 1.5);
    plot([x_pos-box_width/4, x_pos+box_width/4], [upper_whisker, upper_whisker], 'k-', 'LineWidth', 1.5);
end

xlim([0.5, length(class_data)+0.5]);
xticks(1:length(class_data));
xticklabels(class_labels);
ylabel('\gamma_w Distribution', 'FontWeight', 'bold', 'FontSize', 12);
title('Activity Coefficient by Electrolyte Type', 'FontWeight', 'bold');
grid on;

% Plot 4: Molality range comparison by electrolyte type
subplot(2, 3, 4);
hold on; grid on; box on;

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    for s = 1:length(class_salts)
        % Only plot if we have valid data
        if ~isempty(class_salts(s).molality) && ~isempty(class_salts(s).gamma)
            plot(class_salts(s).molality, class_salts(s).gamma, '-', ...
                 'Color', [class_colors(c, :), 0.6], 'LineWidth', 1.5);
        end
    end
    
    % Add dummy for legend
    plot(NaN, NaN, '-', 'Color', class_colors(c, :), 'LineWidth', 3, ...
         'DisplayName', electrolyte_classes{c});
end

xlabel('Molality (mol/kg)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('\gamma_w', 'FontWeight', 'bold', 'FontSize', 12);
title('Activity Coefficient Curves by Electrolyte Type', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;

% Plot 5: Ion number effect on RH range
subplot(2, 3, 5);
hold on; grid on; box on;

total_ions_all = [ion_salts.total_ions];
RH_range_widths = arrayfun(@(x) diff(x.RH_range), ion_salts);

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    
    scatter(total_ions_all(class_mask), RH_range_widths(class_mask), 150, ...
            class_colors(c, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', electrolyte_classes{c});
end

xlabel('Total Number of Ions', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('RH Range Width', 'FontWeight', 'bold', 'FontSize', 12);
title('Operating Range vs Ion Stoichiometry', 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;

% Plot 6: Heatmap of mean gamma by cation/anion numbers
subplot(2, 3, 6);
hold on; box on;

max_cat = max([ion_salts.n_cations]);
max_an = max([ion_salts.n_anions]);

% Create grid - use sum and count approach instead of nanmean
gamma_sum_grid = zeros(max_cat, max_an);
count_grid = zeros(max_cat, max_an);

for i = 1:length(ion_salts)
    nc = ion_salts(i).n_cations;
    na = ion_salts(i).n_anions;
    gamma_sum_grid(nc, na) = gamma_sum_grid(nc, na) + ion_salts(i).mean_gamma;
    count_grid(nc, na) = count_grid(nc, na) + 1;
end

% Calculate mean where we have data
gamma_grid = nan(max_cat, max_an);
for nc = 1:max_cat
    for na = 1:max_an
        if count_grid(nc, na) > 0
            gamma_grid(nc, na) = gamma_sum_grid(nc, na) / count_grid(nc, na);
        end
    end
end

% Plot heatmap
imagesc(gamma_grid');
colorbar;
colormap(jet);

% Add text annotations
for nc = 1:max_cat
    for na = 1:max_an
        if count_grid(nc, na) > 0
            text(nc, na, sprintf('%.3f\n(n=%d)', gamma_grid(nc, na), count_grid(nc, na)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end
end

xlabel('Number of Cations (ν_{+})', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Number of Anions (ν_{-})', 'FontWeight', 'bold', 'FontSize', 12);
title('Mean \gamma_w Heatmap by Ion Numbers', 'FontWeight', 'bold');
xticks(1:max_cat);
yticks(1:max_an);
axis equal tight;

sgtitle('Electrolyte Type and Ion Stoichiometry Analysis', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_electrolyte_type.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_electrolyte_type.fig'));

% Additional figure: Activity Coefficient vs RH by Electrolyte Type
figure('Position', [100, 100, 1600, 600]);

% Plot 1: All data overlaid, color by electrolyte type
subplot(1, 2, 1);
hold on; grid on; box on;

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    % Plot all salts in this class
    for s = 1:length(class_salts)
        % Only plot if we have valid data
        if ~isempty(class_salts(s).RH) && ~isempty(class_salts(s).gamma)
            plot(class_salts(s).RH, class_salts(s).gamma, '-', ...
                 'Color', [class_colors(c, :), 0.5], 'LineWidth', 1.5);
        end
    end
    
    % Add thick dummy line for legend
    plot(NaN, NaN, '-', 'Color', class_colors(c, :), 'LineWidth', 4, ...
         'DisplayName', sprintf('%s (n=%d)', electrolyte_classes{c}, sum(class_mask)));
end

% Add ideal line
plot([0, 1], [1, 1], 'k--', 'LineWidth', 2.5, 'DisplayName', 'Ideal (\gamma_w = 1)');

xlabel('RH (Water Activity)', 'FontWeight', 'bold', 'FontSize', 13);
ylabel('\gamma_w (Activity Coefficient)', 'FontWeight', 'bold', 'FontSize', 13);
title('Activity Coefficient vs RH by Electrolyte Type', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
xlim([0, 1]);
ylim([0.4, 1.1]);
set(gca, 'FontSize', 11);

% Plot 2: Mean trends by electrolyte type
subplot(1, 2, 2);
hold on; grid on; box on;

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    % Create binned average for smoother visualization
    RH_bins = linspace(0, 1, 50);
    gamma_mean = zeros(size(RH_bins));
    gamma_std = zeros(size(RH_bins));
    gamma_count = zeros(size(RH_bins));
    
    for s = 1:length(class_salts)
        for i = 1:length(RH_bins)-1
            mask = class_salts(s).RH >= RH_bins(i) & class_salts(s).RH < RH_bins(i+1);
            if sum(mask) > 0
                gamma_mean(i) = gamma_mean(i) + sum(class_salts(s).gamma(mask));
                gamma_count(i) = gamma_count(i) + sum(mask);
            end
        end
    end
    
    % Calculate average where we have data
    valid_bins = gamma_count > 0;
    gamma_mean(valid_bins) = gamma_mean(valid_bins) ./ gamma_count(valid_bins);
    gamma_mean(~valid_bins) = NaN;
    
    % Plot mean trend with thicker line
    plot(RH_bins, gamma_mean, '-', 'Color', class_colors(c, :), ...
         'LineWidth', 3.5, 'DisplayName', electrolyte_classes{c});
end

% Add ideal line
plot([0, 1], [1, 1], 'k--', 'LineWidth', 2.5, 'DisplayName', 'Ideal (\gamma_w = 1)');

xlabel('RH (Water Activity)', 'FontWeight', 'bold', 'FontSize', 13);
ylabel('Mean \gamma_w', 'FontWeight', 'bold', 'FontSize', 13);
title('Mean Trends by Electrolyte Type', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
xlim([0, 1]);
ylim([0.4, 1.1]);
set(gca, 'FontSize', 11);

sgtitle('Activity Coefficient vs RH Analysis by Electrolyte Type', ...
        'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_electrolyte_gamma_vs_RH.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_electrolyte_gamma_vs_RH.fig'));

% Print statistics by electrolyte type
write_output('\nELECTROLYTE TYPE STATISTICS:\n');
write_output('%-10s %-15s %-15s %-15s %-15s\n', 'Type', 'Count', 'Mean γ_w', 'Std γ_w', 'Mean RH Range');
write_output('%s\n', repmat('-', 1, 70));

for c = 1:length(electrolyte_classes)
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    all_gamma = [];
    all_rh_ranges = [];
    for s = 1:length(class_salts)
        all_gamma = [all_gamma; class_salts(s).gamma];
        all_rh_ranges(end+1) = diff(class_salts(s).RH_range);
    end
    
    write_output('%-10s %-15d %-15.4f %-15.4f %-15.4f\n', ...
            electrolyte_classes{c}, length(class_salts), ...
            mean(all_gamma), std(all_gamma), mean(all_rh_ranges));
end
write_output('\n');

% Print detailed breakdown by electrolyte type
write_output('DETAILED SALT BREAKDOWN BY ELECTROLYTE TYPE:\n');
for c = 1:length(electrolyte_classes)
    write_output('\n%s Electrolytes:\n', upper(electrolyte_classes{c}));
    write_output('%-20s %-8s %-8s %-8s %-12s %-12s\n', ...
            'Salt', 'ν_+', 'ν_-', 'ν_total', 'Mean γ_w', 'RH Range');
    write_output('%s\n', repmat('-', 1, 75));
    
    class_mask = strcmp({ion_salts.electrolyte_class}, electrolyte_classes{c});
    class_salts = ion_salts(class_mask);
    
    for s = 1:length(class_salts)
        write_output('%-20s %-8d %-8d %-8d %-12.4f [%.3f-%.3f]\n', ...
                class_salts(s).name, ...
                class_salts(s).n_cations, ...
                class_salts(s).n_anions, ...
                class_salts(s).total_ions, ...
                class_salts(s).mean_gamma, ...
                class_salts(s).RH_range(1), ...
                class_salts(s).RH_range(2));
    end
end
write_output('\n');

%% ANALYSIS 9: Multi-dimensional Summary
write_output('=== Analysis 9: Multi-dimensional Summary ===\n');

figure('Position', [100, 100, 1600, 500]);

% Filter to salts with complete data
complete_data_mask = ~isnan([salt_analysis.mean_charge_density]) & ...
                     ~isnan([salt_analysis.delta_H_soln]);
comp_salts = salt_analysis(complete_data_mask);

if length(comp_salts) > 0
    % Plot 1: 3D scatter - MW vs Charge Density vs Gamma (colored by enthalpy)
    subplot(1, 3, 1);
    hold on; grid on; box on;
    
    MW_comp = [comp_salts.MW];
    cd_comp = [comp_salts.mean_charge_density];
    gamma_comp = [comp_salts.mean_gamma];
    H_comp = [comp_salts.delta_H_soln];
    
    scatter3(MW_comp, cd_comp, gamma_comp, 100, H_comp, 'filled', ...
             'MarkerEdgeColor', 'k');
    colorbar;
    colormap(jet);
    
    xlabel('MW (g/mol)', 'FontWeight', 'bold');
    ylabel('Charge Density (e/pm^3)', 'FontWeight', 'bold');
    zlabel('Mean \gamma_w', 'FontWeight', 'bold');
    title('3D: MW vs CD vs \gamma_w', 'FontWeight', 'bold');
    view(45, 30);
    
    % Plot 2: Charge density vs enthalpy (colored by gamma)
    subplot(1, 3, 2);
    hold on; grid on; box on;
    
    scatter(cd_comp, H_comp, 100, gamma_comp, 'filled', 'MarkerEdgeColor', 'k');
    cb = colorbar;
    ylabel(cb, 'Mean \gamma_w', 'FontWeight', 'bold');
    
    xlabel('Mean Charge Density (e/pm^3)', 'FontWeight', 'bold');
    ylabel('\DeltaH_{soln} (kJ/mol)', 'FontWeight', 'bold');
    title('Charge Density vs Enthalpy', 'FontWeight', 'bold');
    
    % Plot 3: Parallel coordinates-style
    subplot(1, 3, 3);
    hold on; box on;
    
    % Normalize features for parallel coordinates
    MW_norm = (MW_comp - min(MW_comp)) / (max(MW_comp) - min(MW_comp));
    cd_norm = (cd_comp - min(cd_comp)) / (max(cd_comp) - min(cd_comp));
    H_norm = (H_comp - min(H_comp)) / (max(H_comp) - min(H_comp));
    gamma_norm = (gamma_comp - min(gamma_comp)) / (max(gamma_comp) - min(gamma_comp));
    
    for i = 1:length(comp_salts)
        plot([1, 2, 3, 4], [MW_norm(i), cd_norm(i), H_norm(i), gamma_norm(i)], ...
             'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
    end
    
    xlim([0.5, 4.5]);
    ylim([0, 1]);
    xticks([1, 2, 3, 4]);
    xticklabels({'MW', 'Charge Density', '\DeltaH_{soln}', '\gamma_w'});
    ylabel('Normalized Value', 'FontWeight', 'bold');
    title('Normalized Feature Comparison', 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
end

sgtitle('Multi-dimensional Feature Analysis', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_multidimensional.png'));
savefig(fullfile(filepath, '..', 'figures', 'exploration_multidimensional.fig'));


%% ANALYSIS 10: Kosmotrope vs Chaotrope Effect at RH 75%
write_output('=== Analysis 10: Kosmotrope vs Chaotrope Effect (RH ~ 0.75) ===\n');
kc_groups = {};
kc_gammas = [];

figure('Position', [100, 100, 1000, 600]);
hold on; grid on; box on;

% 1. Extract Data
for i = 1:length(salt_analysis)
    % Check if mean_gamma exists and is valid (not NaN)
    if ~isnan(salt_analysis(i).mean_gamma)
        kc_groups{end+1} = salt_analysis(i).kc_pair_type;
        kc_gammas(end+1) = salt_analysis(i).mean_gamma;
    end
end

% 2. Manual Plotting (Jittered Scatter + Mean Bar)
if ~isempty(kc_gammas)
    unique_types = unique(kc_groups);
    x_positions = 1:length(unique_types);
    
    % Loop through each group
    for j = 1:length(unique_types)
        type_mask = strcmp(kc_groups, unique_types{j});
        vals = kc_gammas(type_mask);
        
        % A. Scatter points with random jitter
        x_vals = repmat(x_positions(j), size(vals)) + (rand(size(vals))-0.5)*0.2;
        scatter(x_vals, vals, 50, 'filled', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerFaceAlpha', 0.6);
        
        % B. Plot Mean Line (Red horizontal bar)
        mean_val = mean(vals);
        line([x_positions(j)-0.25, x_positions(j)+0.25], [mean_val, mean_val], ...
             'Color', 'r', 'LineWidth', 3);
         
        % Add text for N count
        text(x_positions(j), max(vals) + 0.05, sprintf('n=%d', length(vals)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    % 3. Formatting
    set(gca, 'XTick', x_positions, 'XTickLabel', unique_types);
    ylabel('Mean \gamma_w (All valid RH data)', 'FontWeight', 'bold');
    xlabel('Ion Pair Type (Cation-Anion)', 'FontWeight', 'bold');
    
    title({'Mean Activity Coefficient by Kosmotrope/Chaotrope Pair', ...
           'k=Kosmotrope, c=Chaotrope, u=Unknown (Red Bar = Mean)'}, ...
          'FontWeight', 'bold');
    
    saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_kosmotrope_chaotrope.png'));
end
write_output('Generated Kosmotrope/Chaotrope analysis.\n\n');

%% ANALYSIS 11: Correlation Matrix
write_output('=== Analysis 11: Correlation Matrix ===\n');

% Select target RH for activity coefficient comparison
target_RH = 0.75;
RH_tolerance = 0.05;

% Initialize arrays for correlation analysis
correlation_data = [];
correlation_labels = {};
salt_names_corr = {};

% Features to include:
% 1. Cation radius
% 2. Anion radius
% 3. Mean ionic radius (cat + an)/2
% 4. Cation charge
% 5. Anion charge
% 6. Cation charge density
% 7. Anion charge density
% 8. Mean charge density
% 9. Enthalpy of solution
% 10. Molecular weight
% 11. Ionic strength factor
% 12. Total number of ions
% 13. Activity coefficient at target RH

feature_names = {
    'Cation Radius (pm)', ...
    'Anion Radius (pm)', ...
    'Mean Radius (pm)', ...
    'Radius Difference (pm)', ...
    'Cation Charge', ...
    'Anion Charge', ...
    'Cation Charge Density (e/pm^3)', ...
    'Anion Charge Density (e/pm^3)', ...
    'Mean Charge Density (e/pm^3)', ...
    'ΔH_{soln} (kJ/mol)', ...
    'Molecular Weight (g/mol)', ...
    'Ionic Strength Factor', ...
    'Total Ion Number', ...
    'γ_w at RH≈0.75'
};

write_output('Extracting features for correlation analysis...\n');
write_output('Target RH: %.2f ± %.2f\n', target_RH, RH_tolerance);

% Extract features for each salt
for i = 1:length(salt_analysis)
    % Only include salts with sufficient data
    if ~strcmp(salt_analysis(i).category, 'unknown') && ...
       ~isnan(salt_analysis(i).cat_radius) && ...
       ~isnan(salt_analysis(i).an_radius) && ...
       ~isnan(salt_analysis(i).delta_H_soln) && ...
       ~isnan(salt_analysis(i).ionic_strength_factor) && ...
       ~isempty(salt_analysis(i).RH)
        
        % Find gamma at target RH
        RH_mask = abs(salt_analysis(i).RH - target_RH) < RH_tolerance;
        
        if sum(RH_mask) > 0
            gamma_at_target_RH = mean(salt_analysis(i).gamma(RH_mask));
            
            % Calculate mean radius AND radius difference
            mean_radius = (salt_analysis(i).cat_radius + salt_analysis(i).an_radius) / 2;
            radius_diff = abs(salt_analysis(i).cat_radius - salt_analysis(i).an_radius); % <--- ADD THIS LINE
            
            % Build feature vector
            feature_vec = [
                salt_analysis(i).cat_radius;           % 1
                salt_analysis(i).an_radius;            % 2
                mean_radius;                           % 3
                radius_diff;                           % 4
                salt_analysis(i).cat_charge;           % 5
                salt_analysis(i).an_charge;            % 6
                salt_analysis(i).cat_charge_density;   % 7
                salt_analysis(i).an_charge_density;    % 8
                salt_analysis(i).mean_charge_density;  % 9
                salt_analysis(i).delta_H_soln;         % 10
                salt_analysis(i).MW;                   % 11
                salt_analysis(i).ionic_strength_factor;% 12
                salt_analysis(i).total_ions;           % 13
                gamma_at_target_RH                     % 14
            ];
            
            % Check for any NaN values
            if ~any(isnan(feature_vec))
                correlation_data = [correlation_data; feature_vec'];
                salt_names_corr{end+1} = salt_analysis(i).name;
            end
        end
    end
end

write_output('Successfully extracted features for %d salts\n\n', length(salt_names_corr));

% Calculate correlation matrix
if ~isempty(correlation_data)
    corr_matrix = corrcoef(correlation_data);
    
    % Print correlation matrix to output file
    write_output('CORRELATION MATRIX:\n');
    write_output('(Pearson correlation coefficients)\n\n');
    
    % Print header
    write_output('%-35s', 'Feature');
    for j = 1:length(feature_names)
        write_output(' %6d', j);
    end
    write_output('\n');
    write_output('%s\n', repmat('-', 1, 35 + 7*length(feature_names)));
    
    % Print rows
    for i = 1:length(feature_names)
        write_output('%-2d. %-30s', i, feature_names{i});
        for j = 1:length(feature_names)
            write_output(' %6.3f', corr_matrix(i, j));
        end
        write_output('\n');
    end
    write_output('\n');
    
    % Find and report strongest correlations (excluding diagonal)
    write_output('STRONGEST CORRELATIONS (|r| > 0.7):\n');
    write_output('%-35s %-35s %10s\n', 'Feature 1', 'Feature 2', 'r');
    write_output('%s\n', repmat('-', 1, 82));
    
    strong_corr_count = 0;
    for i = 1:length(feature_names)
        for j = i+1:length(feature_names)
            if abs(corr_matrix(i, j)) > 0.7
                write_output('%-35s %-35s %10.4f\n', ...
                    feature_names{i}, feature_names{j}, corr_matrix(i, j));
                strong_corr_count = strong_corr_count + 1;
            end
        end
    end
    
    if strong_corr_count == 0
        write_output('No correlations with |r| > 0.7 found\n');
    end
    write_output('\n');
    
    % Find correlations with activity coefficient
    write_output('CORRELATIONS WITH ACTIVITY COEFFICIENT (γ_w at RH≈%.2f):\n', target_RH);
    write_output('%-40s %10s\n', 'Feature', 'r');
    write_output('%s\n', repmat('-', 1, 52));
    
    gamma_idx = length(feature_names); % Last feature is gamma
    [~, sort_idx] = sort(abs(corr_matrix(gamma_idx, 1:end-1)), 'descend');
    
    for i = 1:length(sort_idx)
        idx = sort_idx(i);
        write_output('%-40s %10.4f\n', feature_names{idx}, corr_matrix(gamma_idx, idx));
    end
    write_output('\n');
    
    % Create correlation heatmap
    figure('Position', [100, 100, 1400, 1200]);
    
    % Main correlation heatmap
    subplot(2, 2, [1, 2]);
    imagesc(corr_matrix);
    colormap(jet);
    caxis([-1, 1]);
    cb = colorbar;
    ylabel(cb, 'Correlation Coefficient (r)', 'FontWeight', 'bold', 'FontSize', 11);
    
    % Add text annotations
    for i = 1:size(corr_matrix, 1)
        for j = 1:size(corr_matrix, 2)
            % Color text based on background
            if abs(corr_matrix(i, j)) > 0.6
                text_color = 'w';
            else
                text_color = 'k';
            end
            
            text(j, i, sprintf('%.2f', corr_matrix(i, j)), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontSize', 8, ...
                 'FontWeight', 'bold', ...
                 'Color', text_color);
        end
    end
    
    % Set axis labels
    xticks(1:length(feature_names));
    yticks(1:length(feature_names));
    xticklabels(feature_names);
    yticklabels(feature_names);
    xtickangle(45);
    
    title('Correlation Matrix: Salt Properties and Activity Coefficients', ...
          'FontWeight', 'bold', 'FontSize', 14);
    axis equal tight;
    set(gca, 'FontSize', 9);
    
    % Subplot: Bar chart of correlations with gamma
    subplot(2, 2, 3);
    hold on; grid on; box on;
    
    gamma_corr = corr_matrix(gamma_idx, 1:end-1);
    [sorted_corr, sort_idx] = sort(gamma_corr);
    sorted_names = feature_names(sort_idx);
    
    colors_bar = zeros(length(sorted_corr), 3);
    for i = 1:length(sorted_corr)
        if sorted_corr(i) > 0
            colors_bar(i, :) = [0.2, 0.6, 0.2]; % Green for positive
        else
            colors_bar(i, :) = [0.8, 0.2, 0.2]; % Red for negative
        end
    end
    
    for i = 1:length(sorted_corr)
        barh(i, sorted_corr(i), 'FaceColor', colors_bar(i, :), 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    yticks(1:length(sorted_names));
    yticklabels(sorted_names);
    xlabel('Correlation with γ_w', 'FontWeight', 'bold', 'FontSize', 12);
    title(sprintf('Correlations with γ_w (RH≈%.2f)', target_RH), ...
          'FontWeight', 'bold', 'FontSize', 12);
    xlim([-1, 1]);
    plot([0, 0], ylim, 'k--', 'LineWidth', 2);
    set(gca, 'FontSize', 9);
    
    % Subplot: Scatter plot of top correlation with gamma
    subplot(2, 2, 4);
    hold on; grid on; box on;
    
    [max_corr, max_idx] = max(abs(gamma_corr));
    
    % Get electrolyte class for coloring
    electrolyte_colors_corr = [];
    electrolyte_class_corr = {};
    unique_classes_corr = {};
    
    for i = 1:length(salt_names_corr)
        idx = find(strcmp({salt_analysis.name}, salt_names_corr{i}));
        if ~isempty(idx) && isfield(salt_analysis(idx), 'electrolyte_class')
            electrolyte_class_corr{i} = salt_analysis(idx).electrolyte_class;
            
            if ~any(strcmp(unique_classes_corr, electrolyte_class_corr{i}))
                unique_classes_corr{end+1} = electrolyte_class_corr{i};
            end
        else
            electrolyte_class_corr{i} = 'unknown';
        end
    end
    
    class_colors_scatter = get_distinguishable_colors(length(unique_classes_corr));
    
    % Plot each electrolyte class
    for c = 1:length(unique_classes_corr)
        class_mask = strcmp(electrolyte_class_corr, unique_classes_corr{c});
        
        scatter(correlation_data(class_mask, max_idx), ...
                correlation_data(class_mask, gamma_idx), ...
                120, class_colors_scatter(c, :), 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
                'DisplayName', unique_classes_corr{c});
    end
    
    % Add regression line
    p = polyfit(correlation_data(:, max_idx), correlation_data(:, gamma_idx), 1);
    x_fit = linspace(min(correlation_data(:, max_idx)), max(correlation_data(:, max_idx)), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'k--', 'LineWidth', 2.5, 'DisplayName', ...
         sprintf('Linear fit (r=%.3f)', gamma_corr(max_idx)));
    
    xlabel(feature_names{max_idx}, 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('γ_w', 'FontWeight', 'bold', 'FontSize', 12);
    title(sprintf('Strongest Predictor of γ_w (|r|=%.3f)', max_corr), ...
          'FontWeight', 'bold', 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 9);
    set(gca, 'FontSize', 10);
    
    sgtitle('Correlation Analysis of Salt Properties', ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_correlation_matrix.png'));
    savefig(fullfile(filepath, '..', 'figures', 'exploration_correlation_matrix.fig'));
    
    % Additional figure: Correlation matrix for different RH values
    figure('Position', [100, 100, 1600, 1000]);
    
    RH_targets = [0.5, 0.6, 0.7, 0.8, 0.9];
    RH_tol = 0.05;
    
    for rh_idx = 1:length(RH_targets)
        target_rh_val = RH_targets(rh_idx);
        
        % Extract data for this RH
        corr_data_rh = [];
        for i = 1:length(salt_analysis)
            if ~strcmp(salt_analysis(i).category, 'unknown') && ...
               ~isnan(salt_analysis(i).cat_radius) && ...
               ~isnan(salt_analysis(i).an_radius) && ...
               ~isnan(salt_analysis(i).delta_H_soln) && ...
               ~isnan(salt_analysis(i).ionic_strength_factor) && ...
               ~isempty(salt_analysis(i).RH)
                
                RH_mask = abs(salt_analysis(i).RH - target_rh_val) < RH_tol;
                
                if sum(RH_mask) > 0
                    gamma_rh = mean(salt_analysis(i).gamma(RH_mask));
                    mean_radius = (salt_analysis(i).cat_radius + salt_analysis(i).an_radius) / 2;
                    radius_diff = abs(salt_analysis(i).cat_radius - salt_analysis(i).an_radius);
                    
                    feature_vec = [
                        salt_analysis(i).cat_radius;
                        salt_analysis(i).an_radius;
                        mean_radius;
                        radius_diff;
                        salt_analysis(i).cat_charge;
                        salt_analysis(i).an_charge;
                        salt_analysis(i).cat_charge_density;
                        salt_analysis(i).an_charge_density;
                        salt_analysis(i).mean_charge_density;
                        salt_analysis(i).delta_H_soln;
                        salt_analysis(i).MW;
                        salt_analysis(i).ionic_strength_factor;
                        salt_analysis(i).total_ions;
                        gamma_rh
                    ];
                    
                    if ~any(isnan(feature_vec))
                        corr_data_rh = [corr_data_rh; feature_vec'];
                    end
                end
            end
        end
        
        if size(corr_data_rh, 1) >= 3  % Need at least 3 samples
            subplot(2, 3, rh_idx);
            
            corr_mat_rh = corrcoef(corr_data_rh);
            gamma_corr_rh = corr_mat_rh(end, 1:end-1);
            
            [sorted_corr_rh, sort_idx_rh] = sort(gamma_corr_rh);
            sorted_names_rh = feature_names(sort_idx_rh);
            
            colors_bar_rh = zeros(length(sorted_corr_rh), 3);
            for i = 1:length(sorted_corr_rh)
                if sorted_corr_rh(i) > 0
                    colors_bar_rh(i, :) = [0.2, 0.6, 0.2];
                else
                    colors_bar_rh(i, :) = [0.8, 0.2, 0.2];
                end
            end
            
            hold on; grid on; box on;
            for i = 1:length(sorted_corr_rh)
                barh(i, sorted_corr_rh(i), 'FaceColor', colors_bar_rh(i, :), ...
                     'EdgeColor', 'k', 'LineWidth', 1);
            end
            
            yticks(1:length(sorted_names_rh));
            yticklabels(sorted_names_rh);
            xlabel('Correlation', 'FontWeight', 'bold');
            title(sprintf('RH ≈ %.2f (n=%d)', target_rh_val, size(corr_data_rh, 1)), ...
                  'FontWeight', 'bold');
            xlim([-1, 1]);
            plot([0, 0], ylim, 'k--', 'LineWidth', 1.5);
            set(gca, 'FontSize', 8);
        else
            subplot(2, 3, rh_idx);
            text(0.5, 0.5, sprintf('Insufficient data\nfor RH ≈ %.2f', target_rh_val), ...
                 'HorizontalAlignment', 'center', 'FontSize', 12);
            axis off;
        end
    end
    
    sgtitle('Feature Correlations with γ_w at Different RH Values', ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    saveas(gcf, fullfile(filepath, '..', 'figures', 'exploration_correlation_vs_RH.png'));
    savefig(fullfile(filepath, '..', 'figures', 'exploration_correlation_vs_RH.fig'));
    
    write_output('Correlation analysis complete. Generated 2 figures.\n\n');
else
    write_output('WARNING: Insufficient data for correlation analysis.\n\n');
end

%% Generate Summary Statistics Report
write_output('\n=== GENERATING SUMMARY STATISTICS ===\n\n');

write_output('OVERALL DATASET STATISTICS:\n');
write_output('Total number of salts: %d\n', n_salts);
write_output('Total data points: %d\n', height(data_table));
write_output('Average points per salt: %.1f\n', height(data_table) / n_salts);
write_output('\n');

write_output('PROPERTY RANGES:\n');
write_output('Molality: %.2f - %.2f mol/kg\n', min(data_table.Molality_mol_per_kg), ...
        max(data_table.Molality_mol_per_kg));
write_output('RH: %.4f - %.4f\n', min(data_table.RH_Water_Activity), ...
        max(data_table.RH_Water_Activity));
write_output('Activity Coefficient: %.4f - %.4f\n', ...
        min(data_table.Activity_Coefficient_Water), ...
        max(data_table.Activity_Coefficient_Water));
write_output('Molecular Weight: %.2f - %.2f g/mol\n', min([salt_analysis.MW]), ...
        max([salt_analysis.MW]));
write_output('\n');

write_output('CHARGE DENSITY STATISTICS (for salts with data):\n');
valid_cd_mask = ~isnan([salt_analysis.mean_charge_density]);
write_output('Number of salts with charge density data: %d\n', sum(valid_cd_mask));
write_output('Mean charge density range: %.2e - %.2e e/pm^3\n', ...
        min([salt_analysis(valid_cd_mask).mean_charge_density]), ...
        max([salt_analysis(valid_cd_mask).mean_charge_density]));
write_output('\n');

write_output('ENTHALPY OF SOLUTION STATISTICS (for salts with data):\n');
valid_H_mask = ~isnan([salt_analysis.delta_H_soln]);
H_vals = [salt_analysis(valid_H_mask).delta_H_soln];
write_output('Number of salts with enthalpy data: %d\n', sum(valid_H_mask));
write_output('Enthalpy range: %.1f - %.1f kJ/mol\n', min(H_vals), max(H_vals));
write_output('Exothermic salts (ΔH < 0): %d\n', sum(H_vals < 0));
write_output('Endothermic salts (ΔH > 0): %d\n', sum(H_vals > 0));
write_output('Mean enthalpy: %.1f kJ/mol\n', mean(H_vals));
write_output('\n');

write_output('CHEMICAL FAMILY DISTRIBUTION:\n');
cats_all = {salt_analysis.category};
unique_cats_all = unique(cats_all);
for c = 1:length(unique_cats_all)
    count = sum(strcmp(cats_all, unique_cats_all{c}));
    write_output('%s: %d salts\n', unique_cats_all{c}, count);
end
write_output('\n');

write_output('TOP 5 SALTS BY DEVIATION FROM IDEALITY (|γ_w - 1|):\n');
[~, sort_idx] = sort(abs([salt_analysis.mean_deviation_from_ideal]), 'descend');
for i = 1:min(5, length(salt_analysis))
    idx = sort_idx(i);
    write_output('%d. %s: γ_w = %.4f (deviation = %.4f)\n', i, ...
            salt_analysis(idx).name, salt_analysis(idx).mean_gamma, ...
            salt_analysis(idx).mean_deviation_from_ideal);
end
write_output('\n');

%% Save analysis results to file
write_output('=== SAVING ANALYSIS RESULTS ===\n');

% Create table for export
export_data = struct2table(salt_analysis);

% Save as CSV
writetable(export_data, fullfile(filepath, 'salt_properties_analysis.csv'));
write_output('Saved: salt_properties_analysis.csv\n');

% Save workspace
save(fullfile(filepath, 'systematic_exploration_workspace.mat'));
write_output('Saved: systematic_exploration_workspace.mat\n');

write_output('\n========================================\n');
write_output('EXPLORATION COMPLETE!\n');
write_output('========================================\n');
write_output('Generated %d figures\n', 14);
write_output('All figures saved to ../figures/\n');
write_output('Analysis data saved to exploration/\n');
write_output('Text output saved to: systematic_exploration_results.txt\n\n');

% Close output file
if fid_out ~= 1
    fclose(fid_out);
    fprintf('Analysis results saved to: %s\n', output_filename);
else
    fprintf('Analysis completed (output file could not be created).\n');
end

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
