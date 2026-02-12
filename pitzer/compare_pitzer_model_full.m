close all 
clear
clc 

%% 1. CONFIGURATION & RUN LIST
% -------------------------------------------------------------------------
% Path setup
addpath('../calculate_mf');
addpath('../data');
csv_path = '../data/baseline_with_ion_properties_legacy.csv';

% Output Setup
output_dir = '../figures/pitzer_comparisons';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

T = 25;       % Temperature in Celsius
MWw = 18.015; % Molecular weight of water

% -------------------------------------------------------------------------
% DATA INPUT (from canonical load_salt_data)
% -------------------------------------------------------------------------
raw_data = load_salt_data();

%% 2. PARSE DATA & LOAD PITZER CSV
% -------------------------------------------------------------------------
if ~exist(csv_path, 'file'), error('CSV file not found: %s', csv_path); end

fprintf('Loading Pitzer database...\n');
opts = detectImportOptions(csv_path);
opts.VariableNamingRule = 'preserve';
pitzer_table = readtable(csv_path, opts);
db_names = pitzer_table.electrolyte; 

missing_salts = {}; 
salts_data = struct([]); 

fprintf('Processing %d salts...\n', length(raw_data));

for i = 1:length(raw_data)
    row = raw_data{i};
    s = struct(); 
    s.Name     = row{1};
    s.RH_Range = [row{3}, row{4}];
    s.Func     = row{5};
    s.Color    = rand(1, 3) * 0.7; 
    
    % --- Determine Family ---
    name = s.Name;
    if contains(name, 'SO4'),       s.Family = 'Sulfates';
    elseif contains(name, 'ClO'),   s.Family = 'Chlorates';
    elseif contains(name, 'NO3'),   s.Family = 'Nitrates';
    elseif contains(name, 'OH'),    s.Family = 'Hydroxides';
    elseif contains(name, 'Cl'),    s.Family = 'Chlorides';
    elseif contains(name, 'Br'),    s.Family = 'Bromides';
    elseif contains(name, 'I'),     s.Family = 'Iodides';
    else,                           s.Family = 'Other';
    end

    % --- Fetch Parameters ---
    row_idx = find(strcmp(db_names, s.Name), 1);
    
    if isempty(row_idx) && strcmp(s.Name, 'NH42SO4')
        row_idx = find(strcmp(db_names, '(NH4)2SO4'), 1);
    end
    
    if isempty(row_idx)
        missing_salts{end+1} = s.Name; 
        continue; 
    end
    
    db_row = pitzer_table(row_idx, :);
    s.MW = db_row.molecule_molecular_weight;
    
    % --- Pitzer Parameters (Force 0 if NaN) ---
    s.Beta0 = db_row.B_MX_0_original;
    s.Beta1 = db_row.B_MX_1_original;
    s.Beta2 = db_row.B_MX_2_original;
    s.Cphi  = db_row.C_MX_phi_original;
    
    if isnan(s.Beta0), s.Beta0 = 0; end
    if isnan(s.Beta1), s.Beta1 = 0; end
    if isnan(s.Beta2), s.Beta2 = 0; end
    if isnan(s.Cphi),  s.Cphi  = 0; end
    
    % --- DEBUG: Verify NaCl params against literature ---
    if strcmp(s.Name, 'NaCl')
        fprintf('  [DEBUG] NaCl Params: b0=%.4f, b1=%.4f, Cphi=%.4f\n', ...
            s.Beta0, s.Beta1, s.Cphi);
        if s.Beta0 > 1 || s.Beta0 < 0.01
             warning('NaCl Beta0 seems off. Check if CSV columns are correct.');
        end
    end
    
    % --- Calc Nu & Alpha ---
    e_type = char(db_row.electrolyte_type); 
    parts = sscanf(e_type, '%d-%d');
    if isempty(parts), parts = [1; 1]; end 
    
    % Stoichiometric sum (v = vM + vX)
    % Example: 1-1 -> 2 ions. 2-1 -> 3 ions.
    val_lcm = lcm(parts(1), parts(2));
    nu_cat = val_lcm / parts(1);
    nu_an  = val_lcm / parts(2);
    s.nu = nu_cat + nu_an;
    
    % Alpha assignment
    if parts(1) == 2 && parts(2) == 2
        s.Alpha1 = 1.4; s.Alpha2 = 12.0;
    else
        s.Alpha1 = 2.0;
        % Alpha2 only matters if Beta2 is non-zero
        s.Alpha2 = (abs(s.Beta2) > 1e-9) * 12.0;
    end

    % --- Generate Data ---
    s.RH_data = linspace(s.RH_Range(1), s.RH_Range(2), 100);
    
    if exist(s.Func, 'file')
        try
            if nargin(s.Func) == 2
                s.mf_salt = feval(s.Func, s.RH_data, T);
            else
                s.mf_salt = feval(s.Func, s.RH_data);
            end
        catch
             s.mf_salt = zeros(size(s.RH_data));
             for k=1:length(s.RH_data)
                 if nargin(s.Func) == 2
                    s.mf_salt(k) = feval(s.Func, s.RH_data(k), T);
                 else
                    s.mf_salt(k) = feval(s.Func, s.RH_data(k));
                 end
             end
        end
    else
        continue;
    end
    
    s.mf_water = 1 - s.mf_salt;
    s.molality = (s.mf_salt / s.MW) ./ (s.mf_water / 1000);
    
    % --- Pitzer Calc ---
    s.aw_pitzer = zeros(size(s.molality));
    for k = 1:length(s.molality)
        s.aw_pitzer(k) = pitzer_water_activity(...
            s.molality(k), s.nu, s.Beta0, s.Beta1, s.Beta2, s.Cphi, ...
            s.Alpha1, s.Alpha2, T);
    end
    
    s.rel_error = 100 * (s.RH_data - s.aw_pitzer) ./ s.RH_data;
    
    if isempty(salts_data), salts_data = s; else, salts_data(end+1) = s; end
end

%% 3. PLOTTING (With Aggregation)
% -------------------------------------------------------------------------
if isempty(salts_data), error('No data.'); end

all_families = {salts_data.Family};
unique_families = unique(all_families);

% Count members per family
fam_counts = zeros(size(unique_families));
for i = 1:length(unique_families)
    fam_counts(i) = sum(strcmp(all_families, unique_families{i}));
end

% Sort families by size
[sorted_counts, sort_idx] = sort(fam_counts, 'ascend');
sorted_families = unique_families(sort_idx);

% Identify "Small" families (The bottom 4 smallest)
num_small = min(4, length(sorted_families));
small_families = sorted_families(1:num_small);
large_families = sorted_families(num_small+1:end);

% Create Combined Group for plotting
combined_name = strjoin(small_families, ', ');
fprintf('\nCombining small families into one figure: %s\n', combined_name);

% --- Plot Large Families Individually ---
for f = 1:length(large_families)
    plot_family(salts_data, large_families{f}, large_families{f}, output_dir);
end

% --- Plot Combined Small Families ---
plot_family(salts_data, small_families, combined_name, output_dir);

%% 4. SUMMARY
fprintf('\n=== STATISTICAL SUMMARY ===\n');
fprintf('%-10s | %-12s | %-10s\n', 'Salt', 'Mean Err (%)', 'Max Err(%)');
fprintf('------------------------------------------\n');
for i = 1:length(salts_data)
    s = salts_data(i);
    errs = abs(s.rel_error(~isnan(s.rel_error)));
    fprintf('%-10s | %8.3f     | %7.3f\n', s.Name, mean(errs), max(errs));
end

%% HELPER FUNCTIONS
% -------------------------------------------------------------------------
function plot_family(all_data, fam_filter, fig_title, out_dir)
    % fam_filter can be a string (one family) or cell array (multiple)
    if ischar(fam_filter), fam_filter = {fam_filter}; end
    
    % Filter data
    mask = ismember({all_data.Family}, fam_filter);
    fam_salts = all_data(mask);
    n = length(fam_salts);
    if n == 0, return; end
    
    % Setup Figure
    nCols = 3; 
    nRows = ceil(n/nCols);
    fig_h = min(1200, 300 * nRows);
    
    % Safe filename (replace spaces/commas)
    fname_safe = regexprep(fig_title, '[, ]+', '_');
    if length(fname_safe) > 50, fname_safe = 'Combined_Small_Families'; end

    % Plot 1: Water Activity
    f1 = figure('Name', [fig_title ' - aw'], 'Position', [50, 50, 1400, fig_h]);
    for k = 1:n
        s = fam_salts(k);
        subplot(nRows, nCols, k); hold on; grid on; box on;
        plot(s.molality, s.RH_data, '.-', 'LineWidth', 1, 'Color', s.Color);
        plot(s.molality, s.aw_pitzer, 'k--', 'LineWidth', 1.2);
        title(s.Name, 'Interpreter', 'none', 'FontSize', 10);
        if k==1, legend('Data', 'Pitzer','Location','best'); end
    end
    sgtitle(fig_title, 'Interpreter', 'none', 'FontSize', 12);
    save_fig(f1, out_dir, ['Aw_' fname_safe]);
    
    % Plot 2: Residuals
    f2 = figure('Name', [fig_title ' - Resid'], 'Position', [100, 100, 1400, fig_h]);
    for k = 1:n
        s = fam_salts(k);
        subplot(nRows, nCols, k); hold on; grid on; box on;
        plot(s.molality, s.rel_error, '.-', 'LineWidth', 1, 'Color', s.Color);
        yline(0, 'k-');
        title(s.Name, 'Interpreter', 'none', 'FontSize', 10);
        if k==1, ylabel('Error (%)'); end
    end
    sgtitle([fig_title ' (Error)'], 'Interpreter', 'none', 'FontSize', 12);
    save_fig(f2, out_dir, ['Resid_' fname_safe]);
end


%% 4. STATISTICAL SUMMARY
% -------------------------------------------------------------------------
fprintf('\n===========================================================\n');
fprintf('                STATISTICAL SUMMARY                        \n');
fprintf('===========================================================\n');
fprintf('%-10s | %-12s | %-10s | %-10s\n', 'Salt', 'Mean Err (%)', 'RMSE (%)', 'Max Err(%)');
fprintf('-----------|--------------|------------|------------\n');

for i = 1:length(salts_data)
    s = salts_data(i);
    errs = s.rel_error(~isnan(s.rel_error) & ~isinf(s.rel_error));
    fprintf('%-10s | %8.3f     | %7.3f    | %7.3f\n', ...
        s.Name, mean(abs(errs)), rms(errs), max(abs(errs)));
end

fprintf('\n===========================================================\n');
fprintf('                MISSING DATA REPORT                        \n');
fprintf('===========================================================\n');
if isempty(missing_salts)
    fprintf('All salts were found in the database.\n');
else
    fprintf('The following %d salts were NOT found in the CSV:\n', length(missing_salts));
    for i = 1:length(missing_salts)
        fprintf('  - %s\n', missing_salts{i});
    end
end
fprintf('===========================================================\n');

function save_fig(h, dir, name)
    set(h, 'Color', 'w');
    fname = fullfile(dir, name);
    print(h, fname, '-dpng', '-r150');
    close(h);
end