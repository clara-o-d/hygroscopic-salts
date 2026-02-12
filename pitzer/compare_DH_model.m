close all 
clear
clc 

%% 1. CONFIGURATION & RUN LIST
% -------------------------------------------------------------------------
addpath('../calculate_mf');
addpath('../data');
csv_path = '../data/baseline_with_ion_properties_legacy.csv';
output_dir = '../figures/dh_comparison';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Constants for Water at 25 deg C
T = 25;           
MWw = 18.015;     
A_DH = 0.51;      % A parameter (mol^-1/2 dm^3/2) 
B_DH = 3.29;      % B parameter (nm^-1 mol^-1/2 dm^3/2)
m_water_moles = 1000 / MWw; % Moles of water in 1kg (~55.51)

% Default "Zoom" target (in molal)
target_zoom_m = 1.5; 

% -------------------------------------------------------------------------
% DATA INPUT (from canonical load_salt_data)
% -------------------------------------------------------------------------
raw_data = load_salt_data();

%% 2. PROCESSING
% -------------------------------------------------------------------------
if ~exist(csv_path, 'file'), error('CSV file not found: %s', csv_path); end
fprintf('Loading Pitzer database for Ion Properties...\n');
opts = detectImportOptions(csv_path);
opts.VariableNamingRule = 'preserve';
pitzer_table = readtable(csv_path, opts);
db_names = pitzer_table.electrolyte; 
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
    
    % --- DB Lookup ---
    row_idx = find(strcmp(db_names, s.Name), 1);
    if isempty(row_idx) && strcmp(s.Name, 'NH42SO4')
        row_idx = find(strcmp(db_names, '(NH4)2SO4'), 1);
    end
    if isempty(row_idx), continue; end
    
    db_row = pitzer_table(row_idx, :);
    s.MW = db_row.molecule_molecular_weight;
    
    % --- Charge & Stoichiometry ---
    e_type = char(db_row.electrolyte_type); 
    parts = sscanf(e_type, '%d-%d');
    if isempty(parts), parts = [1; 1]; end 
    s.z_cat = parts(1);
    s.z_an  = parts(2);
    val_lcm = lcm(s.z_cat, s.z_an);
    s.nu_cat = val_lcm / s.z_cat;
    s.nu_an  = val_lcm / s.z_an;
    s.nu = s.nu_cat + s.nu_an;
    
    % --- Generate Experimental Data ---
    s.RH_data = linspace(s.RH_Range(1), s.RH_Range(2), 100);
    try
        if nargin(s.Func) == 2
            s.mf_salt = feval(s.Func, s.RH_data, T);
        else
            s.mf_salt = feval(s.Func, s.RH_data);
        end
    catch
         s.mf_salt = zeros(size(s.RH_data));
         for k=1:length(s.RH_data)
             if nargin(s.Func) == 2, s.mf_salt(k) = feval(s.Func, s.RH_data(k), T);
             else, s.mf_salt(k) = feval(s.Func, s.RH_data(k)); end
         end
    end
    
    s.mf_water = 1 - s.mf_salt;
    s.molality = (s.mf_salt / s.MW) ./ (s.mf_water / 1000);
    s.aw_exp = s.RH_data; 
    
    % --- NAIVE FITTING OF a0 (Ion Size Parameter) ---
    fit_idx = min(2, length(s.molality)); 
    m_fit   = s.molality(fit_idx);
    aw_fit  = s.aw_exp(fit_idx);
    
    phi_target = -log(aw_fit) / ((MWw/1000) * s.nu * m_fit);
    
    A_ln = A_DH * 2.303; 
    prod_z = abs(s.z_cat * s.z_an);
    I_fit = 0.5 * (m_fit*s.nu_cat*s.z_cat^2 + m_fit*s.nu_an*s.z_an^2);
    
    obj_fun = @(a0) calc_phi_extended(I_fit, A_ln, B_DH, a0, prod_z) - phi_target;
    
    try
        s.a0_fitted = fzero(obj_fun, [0.01, 3.0]); 
    catch
        s.a0_fitted = 0.4; 
    end
    
    % --- FULL CURVE CALCULATION ---
    s.aw_dh = zeros(size(s.molality));
    
    for k = 1:length(s.molality)
        m = s.molality(k);
        if m < 1e-6
            s.aw_dh(k) = 1.0; 
            continue; 
        end
        
        I = 0.5 * ( (m * s.nu_cat * s.z_cat^2) + (m * s.nu_an * s.z_an^2) );
        phi_calc = calc_phi_extended(I, A_ln, B_DH, s.a0_fitted, prod_z);
        
        ln_aw = - (MWw / 1000) * s.nu * m * phi_calc;
        s.aw_dh(k) = exp(ln_aw);
    end
    
    % --- CALCULATE ACTIVITY COEFFICIENTS (Water) ---
    % xw = n_w / (n_w + n_ions)
    s.xw = m_water_moles ./ (m_water_moles + s.nu * s.molality);
    
    s.gamma_w_exp = s.aw_exp ./ s.xw;
    s.gamma_w_dh  = s.aw_dh  ./ s.xw;

    s.rel_error = 100 * (s.aw_exp - s.aw_dh) ./ s.aw_exp;
    
    if isempty(salts_data), salts_data = s; else, salts_data(end+1) = s; end
end

%% 3. PLOTTING
% -------------------------------------------------------------------------
if isempty(salts_data), error('No data.'); end

all_families = {salts_data.Family};
unique_families = unique(all_families);

% Count and Sort
fam_counts = zeros(size(unique_families));
for i = 1:length(unique_families)
    fam_counts(i) = sum(strcmp(all_families, unique_families{i}));
end
[~, sort_idx] = sort(fam_counts, 'ascend');
sorted_families = unique_families(sort_idx);

num_small = min(4, length(sorted_families));
small_families = sorted_families(1:num_small);
large_families = sorted_families(num_small+1:end);
combined_name = strjoin(small_families, ', ');

% --- PLOT LOOP ---
% 1. Large Families
for f = 1:length(large_families)
    % Standard aw/resid plots
    plot_family(salts_data, large_families{f}, large_families{f}, output_dir, target_zoom_m);
    % New Gamma plot
    plot_gamma(salts_data, large_families{f}, large_families{f}, output_dir, target_zoom_m);
end

% 2. Combined Small Families
plot_family(salts_data, small_families, combined_name, output_dir, target_zoom_m);
plot_gamma(salts_data, small_families, combined_name, output_dir, target_zoom_m);

%% 4. SUMMARY
fprintf('\n=== EXTENDED DEBYE-HUCKEL SUMMARY (Fitted to 2nd Point) ===\n');
fprintf('%-10s | %-10s | %-12s | %-10s\n', 'Salt', 'a0 (nm)', 'Mean Err (%)', 'Max Err(%)');
fprintf('------------------------------------------------------\n');
for i = 1:length(salts_data)
    s = salts_data(i);
    errs = abs(s.rel_error(~isnan(s.rel_error)));
    fprintf('%-10s | %6.3f     | %8.3f     | %7.3f\n', ...
        s.Name, s.a0_fitted, mean(errs), max(errs));
end

%% HELPER FUNCTIONS
% -------------------------------------------------------------------------
function phi = calc_phi_extended(I, A, B, a, z_prod)
    x = B * a * sqrt(I);
    if x < 1e-4, sigma = 1 - (3*x)/4; 
    else, sigma = (3 / x^3) * ( (1+x) - 2*log(1+x) - (1/(1+x)) ); end
    phi = 1 - (A * z_prod * sqrt(I) / 3) * sigma;
end

function [final_xlim] = get_smart_xlim(all_data_subset, target_zoom)
    % Determine common xlim for the whole figure to keep it consistent? 
    % Or per subplot? The previous code did per subplot. We stick to per subplot.
    % This function is just a logic holder if needed, but we do it inline.
    final_xlim = target_zoom; 
end

function plot_family(all_data, fam_filter, fig_title, out_dir, target_zoom)
    if ischar(fam_filter), fam_filter = {fam_filter}; end
    mask = ismember({all_data.Family}, fam_filter);
    fam_salts = all_data(mask);
    n = length(fam_salts);
    if n == 0, return; end
    
    nCols = 3; nRows = ceil(n/nCols);
    fig_h = min(1200, 300 * nRows);
    fname_safe = regexprep(fig_title, '[, ]+', '_');
    if length(fname_safe)>50, fname_safe='Combined'; end
    
    f1 = figure('Name', [fig_title ' - aw'], 'Position', [50, 50, 1400, fig_h]);
    
    for k = 1:n
        s = fam_salts(k);
        subplot(nRows, nCols, k); hold on; grid on; box on;
        plot(s.molality, s.aw_exp, '.-', 'LineWidth', 1, 'Color', s.Color);
        plot(s.molality, s.aw_dh, 'k--', 'LineWidth', 1.2);
        title(sprintf('%s (a_0=%.2fnm)', s.Name, s.a0_fitted), 'Interpreter', 'none', 'FontSize', 10);
        
        % Smart Zoom
        zoom_boundary = max(target_zoom, max(s.molality) * 0.5);
        final_xlim = min(max(s.molality), zoom_boundary);
        xlim([0, max(final_xlim, 0.1)]);
        
        if k==1, legend('Data', 'DH Extended','Location','best'); end
    end
    sgtitle(fig_title, 'Interpreter', 'none', 'FontSize', 12);
    set(f1, 'Color', 'w');
    print(f1, fullfile(out_dir, ['Aw_' fname_safe]), '-dpng', '-r150');
    close(f1);
end

function plot_gamma(all_data, fam_filter, fig_title, out_dir, target_zoom)
    if ischar(fam_filter), fam_filter = {fam_filter}; end
    mask = ismember({all_data.Family}, fam_filter);
    fam_salts = all_data(mask);
    n = length(fam_salts);
    if n == 0, return; end
    
    nCols = 3; nRows = ceil(n/nCols);
    fig_h = min(1200, 300 * nRows);
    fname_safe = regexprep(fig_title, '[, ]+', '_');
    if length(fname_safe)>50, fname_safe='Combined'; end
    
    f1 = figure('Name', [fig_title ' - Gamma_w'], 'Position', [70, 70, 1400, fig_h]);
    
    for k = 1:n
        s = fam_salts(k);
        subplot(nRows, nCols, k); hold on; grid on; box on;
        
        plot(s.molality, s.gamma_w_exp, '.-', 'LineWidth', 1, 'Color', s.Color);
        plot(s.molality, s.gamma_w_dh, 'k--', 'LineWidth', 1.2);
        
        title(sprintf('%s \\gamma_w', s.Name), 'Interpreter', 'tex', 'FontSize', 10);
        ylabel('\gamma_w = a_w / x_w');
        
        % Smart Zoom
        zoom_boundary = max(target_zoom, max(s.molality) * 0.5);
        final_xlim = min(max(s.molality), zoom_boundary);
        xlim([0, max(final_xlim, 0.1)]);
        
        if k==1, legend('Data', 'DH Extended','Location','best'); end
    end
    sgtitle([fig_title ' - Activity Coeff'], 'Interpreter', 'none', 'FontSize', 12);
    set(f1, 'Color', 'w');
    print(f1, fullfile(out_dir, ['Gamma_' fname_safe]), '-dpng', '-r150');
    close(f1);
end