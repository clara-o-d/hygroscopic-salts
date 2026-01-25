%% Functional Clustering: Cluster Entire Water Activity Curves
% This script implements functional data analysis to cluster water activity
% curves based on their complete functional form, not just pointwise values.
%
% Clustering approaches:
%   1. Functional k-means: Distance metric = integrated squared difference
%   2. Hierarchical clustering: Based on spline/PCHIP coefficients
%
% Expected cluster types:
%   - "Strong binders": Steep early drop in ln(aw)
%   - "Late binders": Flat initially, then sudden drop
%   - "Inefficient salts": Weak depression even at high molality
%
% Physical validation:
%   - Test enrichment of high charge density ions in clusters
%   - Correlation with hydration enthalpy
%   - Analysis by ion type and chemical family
%
% Author: Generated for AWH ML Discussion
% Date: 2026-01-25

clear; clc; close all;

% Load required packages for Octave
pkg load statistics

% Set up graphics for headless operation (no display)
try
    % Try available graphics toolkits
    available_toolkits = available_graphics_toolkits();
    if ~isempty(available_toolkits)
        graphics_toolkit(available_toolkits{1});
    end
    set(0, 'DefaultFigureVisible', 'off');
catch
    warning('Could not set up graphics toolkit. Figures may not be saved.');
end

%% Configuration
N_CLUSTERS = 3;           % Number of clusters to identify
M_GRID_MIN = 0;           % Minimum molality for standardized grid
M_GRID_MAX = 6;           % Maximum molality for standardized grid
N_GRID_POINTS = 200;      % Number of points in standardized grid
MIN_MOLALITY_COVERAGE = 3.0;  % Minimum max molality for inclusion
MIN_DATA_POINTS = 5;      % Minimum data points required
OUTPUT_DIR = '../../../figures/clustering/';
DATA_FILE = '../../../data/water_activity_all_salts_combined.csv';
PROPERTIES_FILE = '../../../data/baseline_with_ion_properties_legacy.csv';

% Create output directory
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

fprintf('\n');
fprintf('================================================================\n');
fprintf('   FUNCTIONAL CLUSTERING OF WATER ACTIVITY CURVES\n');
fprintf('================================================================\n\n');

%% Load Data
fprintf('Loading data...\n');

% Read CSV file manually (Octave-compatible)
fid = fopen(DATA_FILE, 'r');
if fid == -1
    error('Cannot open data file: %s', DATA_FILE);
end

% Read header
header_line = fgetl(fid);
header_cols = strsplit(header_line, ',');

% Find column indices
salt_col = find(strcmp(header_cols, 'Salt'));
rh_col = find(strcmp(header_cols, 'RH_Water_Activity'));
molality_col = find(strcmp(header_cols, 'Molality_mol_per_kg'));

% Read data
data_raw = textscan(fid, repmat('%s', 1, length(header_cols)), 'Delimiter', ',');
fclose(fid);

% Extract relevant columns
salts_all = data_raw{salt_col};
rh_all = cellfun(@str2double, data_raw{rh_col});
molality_all = cellfun(@str2double, data_raw{molality_col});

% Get unique salts
salts = unique(salts_all);
n_salts = length(salts);
fprintf('  - %d unique salts found\n', n_salts);
fprintf('  - %d total data points\n\n', length(salts_all));

% Load ion properties for physical validation
fprintf('Loading ion properties...\n');
fid_props = fopen(PROPERTIES_FILE, 'r');
if fid_props == -1
    warning('Cannot open properties file: %s', PROPERTIES_FILE);
    props_available = false;
    fprintf('  - Properties not loaded (file not accessible)\n\n');
else
    % Read properties header
    props_header = fgetl(fid_props);
    props_cols = strsplit(props_header, ',');
    
    % Find key column indices
    electrolyte_col = find(strcmp(props_cols, 'electrolyte'));
    cation_col = find(strcmp(props_cols, 'cation'));
    anion_col = find(strcmp(props_cols, 'anion'));
    r_M_col = find(strcmp(props_cols, 'r_M_angstrom'));
    r_X_col = find(strcmp(props_cols, 'r_X_angstrom'));
    dG_cation_col = find(strcmp(props_cols, 'cation_1_delta_G_hydration'));
    dG_anion_col = find(strcmp(props_cols, 'anion_1_delta_G_hydration'));
    elec_type_col = find(strcmp(props_cols, 'electrolyte_type'));
    
    % Read properties data
    props_raw = textscan(fid_props, repmat('%s', 1, length(props_cols)), 'Delimiter', ',');
    fclose(fid_props);
    
    % Create properties struct array
    n_props = length(props_raw{1});
    props = struct();
    
    for col_idx = 1:length(props_cols)
        col_name = props_cols{col_idx};
        col_data = props_raw{col_idx};
        
        % Try to convert to numeric if possible
        is_numeric = true;
        for row_idx = 1:length(col_data)
            val = str2double(col_data{row_idx});
            if isnan(val) && ~strcmp(col_data{row_idx}, 'NaN') && ~isempty(col_data{row_idx})
                is_numeric = false;
                break;
            end
        end
        
        if is_numeric
            props.(col_name) = cellfun(@str2double, col_data);
        else
            props.(col_name) = col_data;
        end
    end
    
    props_available = true;
    fprintf('  - Properties loaded for %d salts\n\n', n_props);
end

%% Create Standardized Functional Grid
fprintf('Creating standardized functional representations...\n');
fprintf('  - Grid range: [%.1f, %.1f] mol/kg\n', M_GRID_MIN, M_GRID_MAX);
fprintf('  - Grid points: %d\n', N_GRID_POINTS);
fprintf('  - Minimum coverage required: %.1f mol/kg\n\n', MIN_MOLALITY_COVERAGE);

% Standardized molality grid
m_grid = linspace(M_GRID_MIN, M_GRID_MAX, N_GRID_POINTS)';

% Storage for functional data
ln_aw_curves = NaN(N_GRID_POINTS, n_salts);
salt_names_included = cell(n_salts, 1);
max_molalities = zeros(n_salts, 1);
n_points_per_salt = zeros(n_salts, 1);
spline_coeffs = cell(n_salts, 1);  % Store spline coefficients for method 2

fprintf('Processing salts:\n');
fprintf('--------------------------------------------------\n');

valid_salt_idx = false(n_salts, 1);

for i = 1:n_salts
    salt_name = salts{i};
    
    % Extract data for this salt
    idx = strcmp(salts_all, salt_name);
    m_salt = molality_all(idx);
    aw_salt = rh_all(idx);
    
    % Sort by molality
    [m_sorted, sort_idx] = sort(m_salt);
    aw_sorted = aw_salt(sort_idx);
    ln_aw_sorted = log(aw_sorted);
    
    max_m = max(m_sorted);
    n_pts = length(m_sorted);
    
    % Store basic info
    salt_names_included{i} = salt_name;
    max_molalities(i) = max_m;
    n_points_per_salt(i) = n_pts;
    
    % Check quality criteria
    if n_pts < MIN_DATA_POINTS
        fprintf('[%3d/%3d] %-15s SKIP (only %d points)\n', i, n_salts, salt_name, n_pts);
        continue;
    end
    
    if max_m < MIN_MOLALITY_COVERAGE
        fprintf('[%3d/%3d] %-15s SKIP (max m=%.2f < %.2f)\n', ...
                i, n_salts, salt_name, max_m, MIN_MOLALITY_COVERAGE);
        continue;
    end
    
    % Interpolate onto standardized grid using shape-preserving PCHIP
    try
        % Only interpolate within the data range
        m_grid_valid = (m_grid >= min(m_sorted)) & (m_grid <= max(m_sorted));
        ln_aw_interp = NaN(N_GRID_POINTS, 1);
        
        % Use PCHIP (shape-preserving cubic interpolation)
        ln_aw_interp(m_grid_valid) = pchip(m_sorted, ln_aw_sorted, m_grid(m_grid_valid));
        
        % Store the curve
        ln_aw_curves(:, i) = ln_aw_interp;
        
        % Extract spline coefficients for hierarchical clustering
        pp = pchip(m_sorted, ln_aw_sorted);
        spline_coeffs{i} = pp.coefs(:);  % Flatten coefficient matrix
        
        valid_salt_idx(i) = true;
        fprintf('[%3d/%3d] %-15s OK (m: %.2f-%.2f, %d pts)\n', ...
                i, n_salts, salt_name, min(m_sorted), max_m, n_pts);
        
    catch ME
        fprintf('[%3d/%3d] %-15s ERROR: %s\n', i, n_salts, salt_name, ME.message);
    end
end

% Filter to valid salts only
ln_aw_curves = ln_aw_curves(:, valid_salt_idx);
salt_names_included = salt_names_included(valid_salt_idx);
max_molalities = max_molalities(valid_salt_idx);
n_points_per_salt = n_points_per_salt(valid_salt_idx);
spline_coeffs = spline_coeffs(valid_salt_idx);
n_valid = sum(valid_salt_idx);

fprintf('--------------------------------------------------\n');
fprintf('Valid salts for clustering: %d / %d\n\n', n_valid, n_salts);

if n_valid < N_CLUSTERS
    error('Not enough valid salts (%d) for %d clusters', n_valid, N_CLUSTERS);
end

%% Method 1: Functional K-Means with Integrated Squared Difference
fprintf('================================================================\n');
fprintf('METHOD 1: FUNCTIONAL K-MEANS CLUSTERING\n');
fprintf('================================================================\n');
fprintf('Distance metric: Integrated squared difference of ln(aw) curves\n\n');

% Compute pairwise functional distances
fprintf('Computing pairwise functional distances...\n');
dist_matrix = zeros(n_valid, n_valid);

for i = 1:n_valid
    for j = (i+1):n_valid
        % Find overlapping region where both curves are defined
        valid_both = ~isnan(ln_aw_curves(:, i)) & ~isnan(ln_aw_curves(:, j));
        
        if sum(valid_both) > 0
            % Integrated squared difference over the overlap
            diff_squared = (ln_aw_curves(valid_both, i) - ln_aw_curves(valid_both, j)).^2;
            m_overlap = m_grid(valid_both);
            dm = m_overlap(2) - m_overlap(1);
            
            % Trapezoidal integration
            integrated_dist = trapz(m_overlap, diff_squared);
            
            % Normalize by overlap length to handle different ranges
            % This prevents bias toward salts with larger molality ranges
            overlap_length = max(m_overlap) - min(m_overlap);
            if overlap_length > 0
                integrated_dist = integrated_dist / overlap_length;
            end
            
            dist_matrix(i, j) = sqrt(integrated_dist);
            dist_matrix(j, i) = dist_matrix(i, j);
        else
            % No overlap - assign large distance
            dist_matrix(i, j) = 1e6;
            dist_matrix(j, i) = 1e6;
        end
    end
    
    if mod(i, 10) == 0
        fprintf('  Progress: %d/%d salts\n', i, n_valid);
    end
end

fprintf('Distance matrix computed (%dx%d)\n\n', n_valid, n_valid);

% Perform k-means clustering using the distance matrix
fprintf('Performing k-means clustering (k=%d)...\n', N_CLUSTERS);

% Use k-medoids (PAM) which works directly with distance matrix
% Custom implementation for Octave compatibility
rng(42);  % For reproducibility

% Simple k-medoids implementation
fprintf('  Using custom k-medoids (PAM) algorithm...\n');
best_cost = inf;
best_labels = [];
best_medoids = [];

% Try multiple random initializations
n_replicates = 20;
for rep = 1:n_replicates
    % Random initialization
    medoid_idx = randperm(n_valid, N_CLUSTERS);
    
    max_iter = 100;
    converged = false;
    
    for iter = 1:max_iter
        % Assignment step: assign each point to nearest medoid
        labels = zeros(n_valid, 1);
        for i = 1:n_valid
            [~, labels(i)] = min(dist_matrix(i, medoid_idx));
        end
        
        % Update step: find best medoid for each cluster
        new_medoid_idx = medoid_idx;
        for k = 1:N_CLUSTERS
            cluster_members = find(labels == k);
            if isempty(cluster_members)
                continue;
            end
            
            % Find point in cluster that minimizes total distance to others
            min_total_dist = inf;
            best_medoid = medoid_idx(k);
            for candidate = cluster_members'
                total_dist = sum(dist_matrix(candidate, cluster_members));
                if total_dist < min_total_dist
                    min_total_dist = total_dist;
                    best_medoid = candidate;
                end
            end
            new_medoid_idx(k) = best_medoid;
        end
        
        % Check convergence
        if isequal(sort(medoid_idx), sort(new_medoid_idx))
            converged = true;
            break;
        end
        medoid_idx = new_medoid_idx;
    end
    
    % Compute total cost
    total_cost = 0;
    for i = 1:n_valid
        total_cost = total_cost + dist_matrix(i, medoid_idx(labels(i)));
    end
    
    % Keep best solution
    if total_cost < best_cost
        best_cost = total_cost;
        best_labels = labels;
        best_medoids = medoid_idx;
    end
    
    if mod(rep, 5) == 0
        fprintf('  Replicate %d/%d complete (best cost: %.4f)\n', rep, n_replicates, best_cost);
    end
end

cluster_labels_fkm = best_labels;
centroids_idx = best_medoids;

fprintf('Clustering complete! (final cost: %.4f)\n\n', best_cost);

% Display cluster assignments
fprintf('Cluster sizes:\n');
for k = 1:N_CLUSTERS
    n_in_cluster = sum(cluster_labels_fkm == k);
    fprintf('  Cluster %d: %d salts (%.1f%%)\n', k, n_in_cluster, ...
            100*n_in_cluster/n_valid);
end
fprintf('\n');

%% Method 2: Hierarchical Clustering on Spline Coefficients
fprintf('================================================================\n');
fprintf('METHOD 2: HIERARCHICAL CLUSTERING ON SPLINE COEFFICIENTS\n');
fprintf('================================================================\n\n');

% Pad spline coefficient vectors to same length (use longest)
coeff_lengths = cellfun(@length, spline_coeffs);
max_coeff_length = max(coeff_lengths);

fprintf('Preparing spline coefficient matrix...\n');
fprintf('  - Max coefficient vector length: %d\n', max_coeff_length);

spline_matrix = zeros(n_valid, max_coeff_length);
for i = 1:n_valid
    n_coeffs = length(spline_coeffs{i});
    spline_matrix(i, 1:n_coeffs) = spline_coeffs{i}';
    % Remaining entries stay as 0 (padding)
end

% Standardize coefficients (important since different positions may have different scales)
spline_matrix_std = (spline_matrix - mean(spline_matrix, 1)) ./ (std(spline_matrix, 0, 1) + 1e-10);

fprintf('\nPerforming hierarchical clustering...\n');
fprintf('  - Linkage method: Ward (minimizes within-cluster variance)\n');
fprintf('  - Distance metric: Euclidean\n\n');

% Compute linkage
Z = linkage(spline_matrix_std, 'ward');

% Cut dendrogram to get N_CLUSTERS
cluster_labels_hier = cluster(Z, 'maxclust', N_CLUSTERS);

fprintf('Hierarchical clustering complete!\n\n');

% Display cluster sizes
fprintf('Cluster sizes:\n');
for k = 1:N_CLUSTERS
    n_in_cluster = sum(cluster_labels_hier == k);
    fprintf('  Cluster %d: %d salts (%.1f%%)\n', k, n_in_cluster, ...
            100*n_in_cluster/n_valid);
end
fprintf('\n');

%% Visualize Clusters
fprintf('================================================================\n');
fprintf('VISUALIZATION\n');
fprintf('================================================================\n\n');

% Use functional k-means results for main analysis
cluster_labels = cluster_labels_fkm;

% Create figure showing all curves colored by cluster
try
    figure('Position', [100, 100, 1400, 800]);

% Plot 1: Functional k-means clusters
subplot(2, 3, 1:2);
hold on;
colors = lines(N_CLUSTERS);

for k = 1:N_CLUSTERS
    cluster_idx = find(cluster_labels == k);
    
    for i = 1:length(cluster_idx)
        salt_idx = cluster_idx(i);
        valid_pts = ~isnan(ln_aw_curves(:, salt_idx));
        
        plot(m_grid(valid_pts), ln_aw_curves(valid_pts, salt_idx), ...
            'Color', [colors(k,:), 0.3], 'LineWidth', 1);
    end
end

% Plot cluster centroids (medoids)
for k = 1:N_CLUSTERS
    medoid_idx = centroids_idx(k);
    valid_pts = ~isnan(ln_aw_curves(:, medoid_idx));
    plot(m_grid(valid_pts), ln_aw_curves(valid_pts, medoid_idx), ...
        'Color', colors(k,:), 'LineWidth', 3, 'DisplayName', ...
        sprintf('Cluster %d (n=%d)', k, sum(cluster_labels == k)));
end

xlabel('Molality (mol/kg)', 'FontSize', 12);
ylabel('ln(a_w)', 'FontSize', 12);
title('Functional K-Means Clustering of ln(a_w) Curves', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 10);
grid on;
box on;
hold off;

% Plot 2: Dendrogram for hierarchical clustering
subplot(2, 3, 3);
dendrogram(Z, 0, 'ColorThreshold', median(Z(end-N_CLUSTERS+2:end-N_CLUSTERS+1, 3)));
title('Hierarchical Clustering Dendrogram', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Salt Index', 'FontSize', 10);
ylabel('Distance (Ward)', 'FontSize', 10);
set(gca, 'XTickLabel', []);

% Plot 3-5: Individual cluster profiles with salt names
for k = 1:N_CLUSTERS
    subplot(2, 3, 3 + k);
    hold on;
    
    cluster_idx = find(cluster_labels == k);
    
    % Plot all curves in this cluster
    for i = 1:length(cluster_idx)
        salt_idx = cluster_idx(i);
        valid_pts = ~isnan(ln_aw_curves(:, salt_idx));
        
        plot(m_grid(valid_pts), ln_aw_curves(valid_pts, salt_idx), ...
            'Color', [colors(k,:), 0.5], 'LineWidth', 1.5);
    end
    
    % Highlight medoid
    medoid_idx = centroids_idx(k);
    valid_pts = ~isnan(ln_aw_curves(:, medoid_idx));
    plot(m_grid(valid_pts), ln_aw_curves(valid_pts, medoid_idx), ...
        'k-', 'LineWidth', 3);
    
    xlabel('Molality (mol/kg)', 'FontSize', 10);
    ylabel('ln(a_w)', 'FontSize', 10);
    title(sprintf('Cluster %d (n=%d)\nMedoid: %s', k, length(cluster_idx), ...
                  salt_names_included{medoid_idx}), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    box on;
    hold off;
end

    sgtitle('Functional Clustering of Water Activity Curves', 'FontSize', 16, 'FontWeight', 'bold');

    % Save figure
    saveas(gcf, fullfile(OUTPUT_DIR, 'functional_clustering_overview.png'));
    saveas(gcf, fullfile(OUTPUT_DIR, 'functional_clustering_overview.fig'));
    fprintf('Saved: functional_clustering_overview.png/fig\n\n');
catch ME
    fprintf('Warning: Could not create overview figure: %s\n\n', ME.message);
end

%% Characterize Cluster Functional Profiles
fprintf('================================================================\n');
fprintf('CLUSTER CHARACTERIZATION\n');
fprintf('================================================================\n\n');

cluster_types = cell(N_CLUSTERS, 1);

for k = 1:N_CLUSTERS
    fprintf('--- Cluster %d ---\n', k);
    cluster_idx = find(cluster_labels == k);
    
    % Compute average curve for this cluster
    cluster_curves = ln_aw_curves(:, cluster_idx);
    avg_curve = mean(cluster_curves, 2, 'omitnan');
    
    % Analyze functional shape
    valid_avg = ~isnan(avg_curve);
    m_valid = m_grid(valid_avg);
    ln_aw_valid = avg_curve(valid_avg);
    
    if length(m_valid) >= 10
        % Initial slope (first 20% of data)
        n_initial = max(3, floor(0.2 * length(m_valid)));
        p_init = polyfit(m_valid(1:n_initial), ln_aw_valid(1:n_initial), 1);
        initial_slope = p_init(1);
        
        % Late slope (last 20% of data)
        n_late = max(3, floor(0.2 * length(m_valid)));
        idx_late = (length(m_valid) - n_late + 1):length(m_valid);
        p_late = polyfit(m_valid(idx_late), ln_aw_valid(idx_late), 1);
        late_slope = p_late(1);
        
        % Total depression at m=6 (if available)
        if max(m_valid) >= 5.5
            ln_aw_at_6 = interp1(m_valid, ln_aw_valid, 6, 'linear', 'extrap');
            total_depression = -ln_aw_at_6;  % More positive = more depression
        else
            ln_aw_at_max = ln_aw_valid(end);
            total_depression = -ln_aw_at_max;
        end
        
        % Classify cluster type
        if initial_slope < -0.15 && total_depression > 0.4
            cluster_type = 'STRONG BINDERS';
        elseif abs(initial_slope) < 0.1 && abs(late_slope) > 0.15
            cluster_type = 'LATE BINDERS';
        elseif total_depression < 0.3
            cluster_type = 'INEFFICIENT SALTS';
        else
            cluster_type = 'INTERMEDIATE';
        end
        
        cluster_types{k} = cluster_type;
        
        fprintf('  Type: %s\n', cluster_type);
        fprintf('  Initial slope (m→0): %.4f\n', initial_slope);
        fprintf('  Late slope (high m): %.4f\n', late_slope);
        fprintf('  Total depression: %.4f\n', total_depression);
    else
        cluster_types{k} = 'UNDEFINED';
        fprintf('  Type: UNDEFINED (insufficient data)\n');
    end
    
    % List salts in cluster
    fprintf('  Salts (%d):\n', length(cluster_idx));
    for i = 1:min(10, length(cluster_idx))  % Show first 10
        fprintf('    - %s\n', salt_names_included{cluster_idx(i)});
    end
    if length(cluster_idx) > 10
        fprintf('    ... and %d more\n', length(cluster_idx) - 10);
    end
    fprintf('\n');
end

%% Physical Validation: Ion Properties
fprintf('================================================================\n');
fprintf('PHYSICAL VALIDATION\n');
fprintf('================================================================\n\n');

% Match salts with property database
fprintf('Matching salts with ion property database...\n');

% Create lookup table for cluster assignments
cluster_results = struct();
for i = 1:n_valid
    cluster_results(i).Salt = salt_names_included{i};
    cluster_results(i).Cluster_FKM = cluster_labels(i);
end

% Manual merge with properties (Octave-compatible)
props_matched = struct();
n_matched = 0;

for i = 1:n_valid
    salt_name = salt_names_included{i};
    
    % Find matching row in props
    match_idx = [];
    if props_available && isfield(props, 'electrolyte')
        for j = 1:length(props.electrolyte)
            if strcmp(props.electrolyte{j}, salt_name)
                match_idx = j;
                break;
            end
        end
    end
    
    % Copy cluster info
    props_matched(i).Salt = salt_name;
    props_matched(i).Cluster_FKM = cluster_labels(i);
    
    % Copy properties if match found
    if ~isempty(match_idx)
        n_matched = n_matched + 1;
        
        % Copy all property fields
        prop_fields = fieldnames(props);
        for f = 1:length(prop_fields)
            field_name = prop_fields{f};
            if iscell(props.(field_name))
                props_matched(i).(field_name) = props.(field_name){match_idx};
            else
                props_matched(i).(field_name) = props.(field_name)(match_idx);
            end
        end
    else
        % Fill with NaN/empty for missing matches
        props_matched(i).cation = '';
        props_matched(i).anion = '';
        props_matched(i).r_M_angstrom = NaN;
        props_matched(i).r_X_angstrom = NaN;
        props_matched(i).cation_1_delta_G_hydration = NaN;
        props_matched(i).anion_1_delta_G_hydration = NaN;
        props_matched(i).electrolyte_type = '';
    end
end

fprintf('  - Matched %d / %d salts with property data\n\n', n_matched, n_valid);

if n_matched < n_valid * 0.5
    warning('Less than 50%% of salts matched with property database');
end

% Analyze charge density enrichment
fprintf('Testing charge density enrichment in clusters...\n\n');

try
    figure('Position', [100, 100, 1400, 500]);

for k = 1:N_CLUSTERS
    % Extract cluster members
    cluster_idx_for_props = find([props_matched.Cluster_FKM] == k);
    cluster_props = props_matched(cluster_idx_for_props);
    
    if length(cluster_props) < 2
        fprintf('Cluster %d: Insufficient data for analysis\n\n', k);
        continue;
    end
    
    % Calculate charge densities (charge / radius)
    n_cluster = length(cluster_props);
    cation_radius = zeros(n_cluster, 1);
    anion_radius = zeros(n_cluster, 1);
    cation_charge = ones(n_cluster, 1);  % Default
    anion_charge = ones(n_cluster, 1);
    
    for i = 1:n_cluster
        % Extract radii
        if isfield(cluster_props, 'r_M_angstrom') && ~isempty(cluster_props(i).r_M_angstrom)
            cation_radius(i) = cluster_props(i).r_M_angstrom;
        else
            cation_radius(i) = NaN;
        end
        
        if isfield(cluster_props, 'r_X_angstrom') && ~isempty(cluster_props(i).r_X_angstrom)
            anion_radius(i) = cluster_props(i).r_X_angstrom;
        else
            anion_radius(i) = NaN;
        end
        
        % Parse charges from electrolyte type (e.g., "1-1", "2-1")
        if isfield(cluster_props, 'electrolyte_type') && ~isempty(cluster_props(i).electrolyte_type)
            type_str = cluster_props(i).electrolyte_type;
            if ischar(type_str) && contains(type_str, '-')
                parts = strsplit(type_str, '-');
                if length(parts) == 2
                    cation_charge(i) = str2double(parts{1});
                    anion_charge(i) = str2double(parts{2});
                end
            end
        end
    end
    
    % Calculate charge densities
    cation_charge_density = cation_charge ./ cation_radius;
    anion_charge_density = anion_charge ./ anion_radius;
    avg_charge_density = (cation_charge_density + anion_charge_density) / 2;
    
    % Plot charge density distribution
    subplot(1, N_CLUSTERS, k);
    valid_cd = avg_charge_density(~isnan(avg_charge_density));
    if ~isempty(valid_cd)
        histogram(valid_cd, 15, 'FaceColor', colors(k,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    xlabel('Charge Density (e/Å)', 'FontSize', 11);
    ylabel('Count', 'FontSize', 11);
    title(sprintf('Cluster %d: %s\n(n=%d)', k, cluster_types{k}, length(cluster_props)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    box on;
    
    % Statistics
    if ~isempty(valid_cd)
        mean_cd = mean(valid_cd);
        median_cd = median(valid_cd);
        
        fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
        fprintf('  Mean charge density: %.4f e/Å\n', mean_cd);
        fprintf('  Median charge density: %.4f e/Å\n', median_cd);
        fprintf('  Range: [%.4f, %.4f] e/Å\n', min(valid_cd), max(valid_cd));
    else
        fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
        fprintf('  No valid charge density data\n');
    end
    fprintf('\n');
end

    sgtitle('Charge Density Distribution by Cluster', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(OUTPUT_DIR, 'cluster_charge_density.png'));
    saveas(gcf, fullfile(OUTPUT_DIR, 'cluster_charge_density.fig'));
    fprintf('Saved: cluster_charge_density.png/fig\n\n');
catch ME
    fprintf('Warning: Could not create charge density figure: %s\n\n', ME.message);
end

%% Hydration Enthalpy Analysis
fprintf('Testing correlation with hydration enthalpy...\n\n');

try
    figure('Position', [100, 100, 1400, 500]);

for k = 1:N_CLUSTERS
    % Extract cluster members
    cluster_idx_for_props = find([props_matched.Cluster_FKM] == k);
    cluster_props = props_matched(cluster_idx_for_props);
    
    if length(cluster_props) < 2
        continue;
    end
    
    % Extract hydration energies
    n_cluster = length(cluster_props);
    cation_dG_hyd = zeros(n_cluster, 1);
    anion_dG_hyd = zeros(n_cluster, 1);
    
    for i = 1:n_cluster
        if isfield(cluster_props, 'cation_1_delta_G_hydration') && ~isempty(cluster_props(i).cation_1_delta_G_hydration)
            cation_dG_hyd(i) = cluster_props(i).cation_1_delta_G_hydration;
        else
            cation_dG_hyd(i) = NaN;
        end
        
        if isfield(cluster_props, 'anion_1_delta_G_hydration') && ~isempty(cluster_props(i).anion_1_delta_G_hydration)
            anion_dG_hyd(i) = cluster_props(i).anion_1_delta_G_hydration;
        else
            anion_dG_hyd(i) = NaN;
        end
    end
    
    % Average hydration free energy (more negative = stronger hydration)
    avg_dG_hyd = (cation_dG_hyd + anion_dG_hyd) / 2;
    
    subplot(1, N_CLUSTERS, k);
    valid_dG = avg_dG_hyd(~isnan(avg_dG_hyd));
    if ~isempty(valid_dG)
        histogram(valid_dG, 15, 'FaceColor', colors(k,:), 'EdgeColor', 'k', 'FaceAlpha', 0.7);
    end
    xlabel('ΔG_{hydration} (kJ/mol)', 'FontSize', 11);
    ylabel('Count', 'FontSize', 11);
    title(sprintf('Cluster %d: %s\n(n=%d)', k, cluster_types{k}, length(cluster_props)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    box on;
    
    % Statistics
    if ~isempty(valid_dG)
        mean_dG = mean(valid_dG);
        median_dG = median(valid_dG);
        
        fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
        fprintf('  Mean ΔG_hydration: %.2f kJ/mol\n', mean_dG);
        fprintf('  Median ΔG_hydration: %.2f kJ/mol\n', median_dG);
        fprintf('  Range: [%.2f, %.2f] kJ/mol\n', min(valid_dG), max(valid_dG));
    else
        fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
        fprintf('  No valid hydration energy data\n');
    end
    fprintf('\n');
end

    sgtitle('Hydration Free Energy Distribution by Cluster', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(OUTPUT_DIR, 'cluster_hydration_energy.png'));
    saveas(gcf, fullfile(OUTPUT_DIR, 'cluster_hydration_energy.fig'));
    fprintf('Saved: cluster_hydration_energy.png/fig\n\n');
catch ME
    fprintf('Warning: Could not create hydration energy figure: %s\n\n', ME.message);
end

%% Ion Type Analysis
fprintf('Analyzing ion type composition by cluster...\n\n');

for k = 1:N_CLUSTERS
    % Extract cluster members
    cluster_idx_for_props = find([props_matched.Cluster_FKM] == k);
    cluster_props = props_matched(cluster_idx_for_props);
    
    if length(cluster_props) < 1
        continue;
    end
    
    fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
    fprintf('  Anion composition:\n');
    
    % Count anions
    anion_list = {};
    anion_count = [];
    for i = 1:length(cluster_props)
        if isfield(cluster_props, 'anion') && ~isempty(cluster_props(i).anion)
            anion = cluster_props(i).anion;
            
            % Find or add anion
            idx = find(strcmp(anion_list, anion));
            if isempty(idx)
                anion_list{end+1} = anion;
                anion_count(end+1) = 1;
            else
                anion_count(idx) = anion_count(idx) + 1;
            end
        end
    end
    
    % Sort and display
    [anion_count_sorted, sort_idx] = sort(anion_count, 'descend');
    anion_list_sorted = anion_list(sort_idx);
    for i = 1:length(anion_list_sorted)
        pct = 100 * anion_count_sorted(i) / length(cluster_props);
        fprintf('    %s: %d (%.1f%%)\n', anion_list_sorted{i}, anion_count_sorted(i), pct);
    end
    fprintf('\n');
end

% Cation analysis
for k = 1:N_CLUSTERS
    % Extract cluster members
    cluster_idx_for_props = find([props_matched.Cluster_FKM] == k);
    cluster_props = props_matched(cluster_idx_for_props);
    
    if length(cluster_props) < 1
        continue;
    end
    
    fprintf('--- Cluster %d: %s ---\n', k, cluster_types{k});
    fprintf('  Cation composition:\n');
    
    % Count cations
    cation_list = {};
    cation_count = [];
    for i = 1:length(cluster_props)
        if isfield(cluster_props, 'cation') && ~isempty(cluster_props(i).cation)
            cation = cluster_props(i).cation;
            
            % Find or add cation
            idx = find(strcmp(cation_list, cation));
            if isempty(idx)
                cation_list{end+1} = cation;
                cation_count(end+1) = 1;
            else
                cation_count(idx) = cation_count(idx) + 1;
            end
        end
    end
    
    % Sort and display
    [cation_count_sorted, sort_idx] = sort(cation_count, 'descend');
    cation_list_sorted = cation_list(sort_idx);
    for i = 1:length(cation_list_sorted)
        pct = 100 * cation_count_sorted(i) / length(cluster_props);
        fprintf('    %s: %d (%.1f%%)\n', cation_list_sorted{i}, cation_count_sorted(i), pct);
    end
    fprintf('\n');
end

%% Export Results
fprintf('================================================================\n');
fprintf('EXPORTING RESULTS\n');
fprintf('================================================================\n\n');

% Save results to CSV (Octave-compatible)
csv_filename = fullfile(OUTPUT_DIR, 'clustering_results.csv');
fid = fopen(csv_filename, 'w');
fprintf(fid, 'Salt,Cluster_FunctionalKMeans,Cluster_Hierarchical,Max_Molality,N_DataPoints,Cluster_Type\n');
for i = 1:n_valid
    k = cluster_labels(i);
    fprintf(fid, '%s,%d,%d,%.4f,%d,%s\n', ...
        salt_names_included{i}, cluster_labels(i), cluster_labels_hier(i), ...
        max_molalities(i), n_points_per_salt(i), cluster_types{k});
end
fclose(fid);
fprintf('Saved: clustering_results.csv\n');

% Save detailed cluster assignments with medoid info
csv_detailed = fullfile(OUTPUT_DIR, 'clustering_detailed.csv');
fid = fopen(csv_detailed, 'w');
fprintf(fid, 'Salt,Cluster,Type,Medoid,DistanceToMedoid\n');
for i = 1:n_valid
    k = cluster_labels(i);
    dist_to_medoid = dist_matrix(i, centroids_idx(k));
    fprintf(fid, '%s,%d,%s,%s,%.6f\n', ...
        salt_names_included{i}, k, cluster_types{k}, ...
        salt_names_included{centroids_idx(k)}, dist_to_medoid);
end
fclose(fid);
fprintf('Saved: clustering_detailed.csv\n\n');

%% Summary Statistics
fprintf('================================================================\n');
fprintf('SUMMARY\n');
fprintf('================================================================\n\n');

fprintf('Clustering Methods Compared:\n');
fprintf('  1. Functional k-means (integrated squared difference)\n');
fprintf('  2. Hierarchical clustering (Ward linkage on spline coefficients)\n\n');

fprintf('Cluster Type Identification:\n');
for k = 1:N_CLUSTERS
    fprintf('  Cluster %d: %s (%d salts)\n', k, cluster_types{k}, sum(cluster_labels == k));
end
fprintf('\n');

% Compare the two clustering methods
agreement = sum(cluster_labels == cluster_labels_hier) / n_valid * 100;
fprintf('Agreement between methods: %.1f%%\n\n', agreement);

fprintf('Physical validation:\n');
fprintf('  - Charge density distributions computed\n');
fprintf('  - Hydration free energy correlations analyzed\n');
fprintf('  - Ion type enrichment examined\n\n');

fprintf('All results saved to: %s\n', OUTPUT_DIR);
fprintf('\n================================================================\n');
fprintf('FUNCTIONAL CLUSTERING COMPLETE\n');
fprintf('================================================================\n\n');
