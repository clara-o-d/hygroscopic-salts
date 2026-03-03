close all
clear
clc

%% Setup paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'calculate_activity_mixtures'));

% Output directory
fig_out_dir = fullfile(filepath, '..', 'figures', 'mixture_facets');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% Constants
T       = 25;           % Temperature (°C)
MWw     = 18.015;       % Molecular weight of water (g/mol)
n_water = 1000 / MWw;  % Moles of water per kg water (~55.51 mol/kg)

y_mol  = 0.5;   % Molar fraction of salt 2 for "constant molality ratio" path
r_mass = 0.5;   % Mass fraction of salt 2 for "constant mass fraction ratio" path
n_pts  = 250;   % Number of points along each mixture curve

%% -----------------------------------------------------------------------
%  Load Pitzer mixing parameters (theta, psi) from ion_activity_data.xlsx
% -----------------------------------------------------------------------
ion_raw = readcell(fullfile(filepath, '..', 'data', 'ion_activity_data.xlsx'));
% Row 1 = header; data cols (1-indexed): 1=System, 2=Exptl, 3=MaxI, 5=theta, 6=psi
ion_map = containers.Map('KeyType','char','ValueType','any');
for ii = 2:size(ion_raw, 1)
    sys = ion_raw{ii, 1};
    if ~(ischar(sys) || isstring(sys)), continue; end
    th  = ion_raw{ii, 5};
    ps  = ion_raw{ii, 6};
    ex  = ion_raw{ii, 2};
    mx  = ion_raw{ii, 3};
    if ~isnumeric(th) || ~isnumeric(mx), continue; end
    if ~isnumeric(ps), ps = NaN; end
    ion_map(char(sys)) = struct('theta', th, 'psi', ps, ...
                                'exptl', char(ex), 'maxI', mx);
end

%% -----------------------------------------------------------------------
%  Define all 20 mixtures
%  Fields per struct:
%    name, display
%    mix_func  - calculate_activity function name (takes mf1, mf2)
%    s1_name, s1_MW, s1_nu, s1_func, s1_needsT, s1_RH  (pure salt 1)
%    s2_name, s2_MW, s2_nu, s2_func, s2_needsT, s2_RH  (pure salt 2)
%    mf1_max, mf2_max  - calibrated upper bounds of mixture function
% -----------------------------------------------------------------------
mixtures = {};

% 1. NaCl + LiCl
mixtures{end+1} = struct( ...
    'name','NaCl_LiCl', 'display','NaCl + LiCl', 'ion_sys','LiCl-NaCl', ...
    'mix_func','calculate_activity_NaCl_LiCl', ...
    's1_name','NaCl',   's1_MW',58.443,  's1_nu',2, 's1_func','calculate_mf_NaCl',  's1_needsT',false, 's1_RH',[0.765 0.99], ...
    's2_name','LiCl',   's2_MW',42.394,  's2_nu',2, 's2_func','calculate_mf_LiCl',  's2_needsT',true,  's2_RH',[0.12  0.97], ...
    'mf1_max',0.2366, 'mf2_max',0.1549);

% 2. NaNO3 + LiNO3
mixtures{end+1} = struct( ...
    'name','NaNO3_LiNO3', 'display','NaNO_3 + LiNO_3', 'ion_sys','LiNO3-NaNO3', ...
    'mix_func','calculate_activity_NaNO3_LiNO3', ...
    's1_name','NaNO3',  's1_MW',85.00,   's1_nu',2, 's1_func','calculate_mf_NaNO3', 's1_needsT',false, 's1_RH',[0.971 0.995], ...
    's2_name','LiNO3',  's2_MW',68.95,   's2_nu',2, 's2_func','calculate_mf_LiNO3', 's2_needsT',false, 's2_RH',[0.736 0.99], ...
    'mf1_max',0.3599, 'mf2_max',0.2471);

% 3. NaC2H3O2 + LiC2H3O2  (no pure-salt calculate_mf functions available)
mixtures{end+1} = struct( ...
    'name','NaC2H3O2_LiC2H3O2', 'display','NaC_2H_3O_2 + LiC_2H_3O_2', 'ion_sys','LiOAc-NaOAc', ...
    'mix_func','calculate_activity_NaC2H3O2_LiC2H3O2', ...
    's1_name','NaAc',   's1_MW',82.034,  's1_nu',2, 's1_func','', 's1_needsT',false, 's1_RH',[NaN NaN], ...
    's2_name','LiAc',   's2_MW',65.98,   's2_nu',2, 's2_func','', 's2_needsT',false, 's2_RH',[NaN NaN], ...
    'mf1_max',0.2437, 'mf2_max',0.2244);

% 4. LiCl + BaCl2
mixtures{end+1} = struct( ...
    'name','LiCl_BaCl2', 'display','LiCl + BaCl_2', 'ion_sys','LiCl-BaCl2', ...
    'mix_func','calculate_activity_LiCl_BaCl2', ...
    's1_name','LiCl',   's1_MW',42.394,  's1_nu',2, 's1_func','calculate_mf_LiCl',  's1_needsT',true,  's1_RH',[0.12  0.97], ...
    's2_name','BaCl2',  's2_MW',208.23,  's2_nu',3, 's2_func','calculate_mf_BaCl2', 's2_needsT',false, 's2_RH',[0.908 0.96], ...
    'mf1_max',0.1620, 'mf2_max',0.2104);

% 5. LiCl + MgCl2
mixtures{end+1} = struct( ...
    'name','LiCl_MgCl2', 'display','LiCl + MgCl_2', 'ion_sys','', ...
    'mix_func','calculate_activity_LiCl_MgCl2', ...
    's1_name','LiCl',   's1_MW',42.394,  's1_nu',2, 's1_func','calculate_mf_LiCl',  's1_needsT',true,  's1_RH',[0.12  0.97], ...
    's2_name','MgCl2',  's2_MW',95.211,  's2_nu',3, 's2_func','calculate_mf_MgCl',  's2_needsT',false, 's2_RH',[0.33  0.97], ...
    'mf1_max',0.2762, 'mf2_max',0.1860);

% 6. LiCl + CaCl2
mixtures{end+1} = struct( ...
    'name','LiCl_CaCl2', 'display','LiCl + CaCl_2', 'ion_sys','', ...
    'mix_func','calculate_activity_LiCl_CaCl2', ...
    's1_name','LiCl',   's1_MW',42.394,  's1_nu',2, 's1_func','calculate_mf_LiCl',  's1_needsT',true,  's1_RH',[0.12  0.97], ...
    's2_name','CaCl2',  's2_MW',110.984, 's2_nu',3, 's2_func','calculate_mf_CaCl',  's2_needsT',true,  's2_RH',[0.31  0.97], ...
    'mf1_max',0.3851, 'mf2_max',0.3879);

% 7. HClO4 + LiClO4
mixtures{end+1} = struct( ...
    'name','HClO4_LiClO4', 'display','HClO_4 + LiClO_4', 'ion_sys','HClO4-LiClO4', ...
    'mix_func','calculate_activity_HClO4_LiClO4', ...
    's1_name','HClO4',  's1_MW',100.46,  's1_nu',2, 's1_func','calculate_mf_HClO4', 's1_needsT',false, 's1_RH',[0.634 0.997], ...
    's2_name','LiClO4', 's2_MW',106.39,  's2_nu',2, 's2_func','calculate_mf_LiClO4','s2_needsT',false, 's2_RH',[0.779 0.987], ...
    'mf1_max',0.2497, 'mf2_max',0.2614);

% 8. HClO4 + NaClO4
mixtures{end+1} = struct( ...
    'name','HClO4_NaClO4', 'display','HClO_4 + NaClO_4', 'ion_sys','HClO4-NaClO4', ...
    'mix_func','calculate_activity_HClO4_NaClO4', ...
    's1_name','HClO4',  's1_MW',100.46,  's1_nu',2, 's1_func','calculate_mf_HClO4', 's1_needsT',false, 's1_RH',[0.634 0.997], ...
    's2_name','NaClO4', 's2_MW',122.44,  's2_nu',2, 's2_func','calculate_mf_NaClO4','s2_needsT',false, 's2_RH',[0.778 0.99], ...
    'mf1_max',0.4191, 'mf2_max',0.5598);

% 9. LiClO4 + NaClO4
mixtures{end+1} = struct( ...
    'name','LiClO4_NaClO4', 'display','LiClO_4 + NaClO_4', 'ion_sys','LiClO4-NaClO4', ...
    'mix_func','calculate_activity_LiClO4_NaClO4', ...
    's1_name','LiClO4', 's1_MW',106.39,  's1_nu',2, 's1_func','calculate_mf_LiClO4','s1_needsT',false, 's1_RH',[0.779 0.987], ...
    's2_name','NaClO4', 's2_MW',122.44,  's2_nu',2, 's2_func','calculate_mf_NaClO4','s2_needsT',false, 's2_RH',[0.778 0.99], ...
    'mf1_max',0.2723, 'mf2_max',0.3502);

% 10. HCl + NaCl
mixtures{end+1} = struct( ...
    'name','HCl_NaCl', 'display','HCl + NaCl', 'ion_sys','HCl-NaCl', ...
    'mix_func','calculate_activity_HCl_NaCl', ...
    's1_name','HCl',    's1_MW',36.461,  's1_nu',2, 's1_func','calculate_mf_HCl',   's1_needsT',false, 's1_RH',[0.17  0.97], ...
    's2_name','NaCl',   's2_MW',58.443,  's2_nu',2, 's2_func','calculate_mf_NaCl',  's2_needsT',false, 's2_RH',[0.765 0.99], ...
    'mf1_max',0.0649, 'mf2_max',0.1347);

% 11. HCl + BaCl2
mixtures{end+1} = struct( ...
    'name','HCl_BaCl2', 'display','HCl + BaCl_2', 'ion_sys','HCl-BaCl2', ...
    'mix_func','calculate_activity_HCl_BaCl2', ...
    's1_name','HCl',    's1_MW',36.461,  's1_nu',2, 's1_func','calculate_mf_HCl',   's1_needsT',false, 's1_RH',[0.17  0.97], ...
    's2_name','BaCl2',  's2_MW',208.23,  's2_nu',3, 's2_func','calculate_mf_BaCl2', 's2_needsT',false, 's2_RH',[0.908 0.96], ...
    'mf1_max',0.0650, 'mf2_max',0.2071);

% 12. HCl + NaClO4
mixtures{end+1} = struct( ...
    'name','HCl_NaClO4', 'display','HCl + NaClO_4', 'ion_sys','', ...
    'mix_func','calculate_activity_HCl_NaClO4', ...
    's1_name','HCl',    's1_MW',36.461,  's1_nu',2, 's1_func','calculate_mf_HCl',   's1_needsT',false, 's1_RH',[0.17  0.97], ...
    's2_name','NaClO4', 's2_MW',122.44,  's2_nu',2, 's2_func','calculate_mf_NaClO4','s2_needsT',false, 's2_RH',[0.778 0.99], ...
    'mf1_max',0.0182, 'mf2_max',0.6272);

% 13. HCl + Ba(ClO4)2
mixtures{end+1} = struct( ...
    'name','HCl_BaClO42', 'display','HCl + Ba(ClO_4)_2', 'ion_sys','', ...
    'mix_func','calculate_activity_HCl_BaClO42', ...
    's1_name','HCl',    's1_MW',36.461,  's1_nu',2, 's1_func','calculate_mf_HCl',   's1_needsT',false, 's1_RH',[0.17  0.97], ...
    's2_name','BaClO42','s2_MW',336.23,  's2_nu',3, 's2_func','calculate_mf_BaClO42','s2_needsT',false,'s2_RH',[0.561 0.995], ...
    'mf1_max',0.0506, 'mf2_max',0.5049);

% 14. HCl + Na2SO4
mixtures{end+1} = struct( ...
    'name','HCl_Na2SO4', 'display','HCl + Na_2SO_4', 'ion_sys','', ...
    'mix_func','calculate_activity_HCl_Na2SO4', ...
    's1_name','HCl',    's1_MW',36.461,  's1_nu',2, 's1_func','calculate_mf_HCl',   's1_needsT',false, 's1_RH',[0.17  0.97], ...
    's2_name','Na2SO4', 's2_MW',142.04,  's2_nu',3, 's2_func','calculate_mf_Na2SO4','s2_needsT',false, 's2_RH',[0.90  0.995], ...
    'mf1_max',0.0344, 'mf2_max',0.2073);

% 15. LiCl + KCl
mixtures{end+1} = struct( ...
    'name','LiCl_KCl', 'display','LiCl + KCl', 'ion_sys','LiCl-KCl', ...
    'mix_func','calculate_activity_LiCl_KCl', ...
    's1_name','LiCl',   's1_MW',42.394,  's1_nu',2, 's1_func','calculate_mf_LiCl',  's1_needsT',true,  's1_RH',[0.12  0.97], ...
    's2_name','KCl',    's2_MW',74.551,  's2_nu',2, 's2_func','calculate_mf_KCl',   's2_needsT',false, 's2_RH',[0.855 0.99], ...
    'mf1_max',0.1450, 'mf2_max',0.2297);

% 16. LiCl + CsCl
mixtures{end+1} = struct( ...
    'name','LiCl_CsCl', 'display','LiCl + CsCl', 'ion_sys','LiCl-CsCl', ...
    'mix_func','calculate_activity_LiCl_CsCl', ...
    's1_name','LiCl',   's1_MW',42.394,  's1_nu',2, 's1_func','calculate_mf_LiCl',  's1_needsT',true,  's1_RH',[0.12  0.97], ...
    's2_name','CsCl',   's2_MW',168.363, 's2_nu',2, 's2_func','calculate_mf_CsCl',  's2_needsT',false, 's2_RH',[0.82  0.99], ...
    'mf1_max',0.2028, 'mf2_max',0.5025);

% 17. NaCl + CsCl
mixtures{end+1} = struct( ...
    'name','NaCl_CsCl', 'display','NaCl + CsCl', 'ion_sys','NaCl-CsCl', ...
    'mix_func','calculate_activity_NaCl_CsCl', ...
    's1_name','NaCl',   's1_MW',58.443,  's1_nu',2, 's1_func','calculate_mf_NaCl',  's1_needsT',false, 's1_RH',[0.765 0.99], ...
    's2_name','CsCl',   's2_MW',168.363, 's2_nu',2, 's2_func','calculate_mf_CsCl',  's2_needsT',false, 's2_RH',[0.82  0.99], ...
    'mf1_max',0.2375, 'mf2_max',0.4730);

% 18. MgCl2 + NaCl
mixtures{end+1} = struct( ...
    'name','MgCl2_NaCl', 'display','MgCl_2 + NaCl', 'ion_sys','NaCl-MgCl2', ...
    'mix_func','calculate_activity_MgCl2_NaCl', ...
    's1_name','MgCl2',  's1_MW',95.211,  's1_nu',3, 's1_func','calculate_mf_MgCl',  's1_needsT',false, 's1_RH',[0.33  0.97], ...
    's2_name','NaCl',   's2_MW',58.443,  's2_nu',2, 's2_func','calculate_mf_NaCl',  's2_needsT',false, 's2_RH',[0.765 0.99], ...
    'mf1_max',0.1033, 'mf2_max',0.1746);

% 19. MgCl2 + CaCl2
mixtures{end+1} = struct( ...
    'name','MgCl2_CaCl2', 'display','MgCl_2 + CaCl_2', 'ion_sys','MgCl2-CaCl2', ...
    'mix_func','calculate_activity_MgCl2_CaCl2', ...
    's1_name','MgCl2',  's1_MW',95.211,  's1_nu',3, 's1_func','calculate_mf_MgCl',  's1_needsT',false, 's1_RH',[0.33  0.97], ...
    's2_name','CaCl2',  's2_MW',110.984, 's2_nu',3, 's2_func','calculate_mf_CaCl',  's2_needsT',true,  's2_RH',[0.31  0.97], ...
    'mf1_max',0.2222, 'mf2_max',0.2498);

% 20. NH4Cl + LiCl
mixtures{end+1} = struct( ...
    'name','NH4Cl_LiCl', 'display','NH_4Cl + LiCl', 'ion_sys','', ...
    'mix_func','calculate_activity_NH4Cl_LiCl', ...
    's1_name','NH4Cl',  's1_MW',53.491,  's1_nu',2, 's1_func','calculate_mf_NH4Cl', 's1_needsT',false, 's1_RH',[0.815 0.99], ...
    's2_name','LiCl',   's2_MW',42.394,  's2_nu',2, 's2_func','calculate_mf_LiCl',  's2_needsT',true,  's2_RH',[0.12  0.97], ...
    'mf1_max',0.1762, 'mf2_max',0.1450);

%% -----------------------------------------------------------------------
%  Plot loop
% -----------------------------------------------------------------------
fprintf('Generating faceted plots for %d mixtures...\n', length(mixtures));

warning('off', 'all');

for mi = 1:length(mixtures)
    mix = mixtures{mi};
    fprintf('  [%2d/%2d]  %s\n', mi, length(mixtures), mix.display);

    % --- pure salt curves ---
    [m_s1, aw_s1, xw_s1, mfw_s1, gw_s1] = ...
        compute_pure_salt(mix.s1_func, mix.s1_MW, mix.s1_nu, ...
                          mix.s1_needsT, T, mix.s1_RH, n_water, n_pts);
    [m_s2, aw_s2, xw_s2, mfw_s2, gw_s2] = ...
        compute_pure_salt(mix.s2_func, mix.s2_MW, mix.s2_nu, ...
                          mix.s2_needsT, T, mix.s2_RH, n_water, n_pts);

    % --- mixture paths ---
    [m_mx, aw_mx, xw_mx, mfw_mx, gw_mx] = ...
        compute_mix_const_ymol(mix.mix_func, ...
            mix.s1_MW, mix.s1_nu, mix.s2_MW, mix.s2_nu, ...
            y_mol, mix.mf1_max, mix.mf2_max, n_water, n_pts);
    [m_mr, aw_mr, xw_mr, mfw_mr, gw_mr] = ...
        compute_mix_const_rmass(mix.mix_func, ...
            mix.s1_MW, mix.s1_nu, mix.s2_MW, mix.s2_nu, ...
            r_mass, mix.mf1_max, mix.mf2_max, n_water, n_pts);

    % -------------------------------------------------------------------
    %  Additive isotherms (rows 2 & 3)
    %
    %  Evaluated on a common a_w grid spanning the overlap of the two
    %  pure-salt ranges, so the additive curve lands BETWEEN the two
    %  pure-salt curves on the plot.
    %
    %  x_w additive:   mole-fraction weighted avg at each a_w
    %    xw_add  = (1-y_mol)*xw1(aw) + y_mol*xw2(aw)
    %
    %  mf_w additive:  mass-conservative W-space weighted avg at each a_w
    %    W_i = mfw_i/(1-mfw_i)  (kg water per kg dry salt i)
    %    W_add = (1-r_mass)*W1(aw) + r_mass*W2(aw)
    %    mfw_add = W_add / (1 + W_add)
    %
    %  Both are strictly between the two pure-salt values at every a_w.
    % -------------------------------------------------------------------
    has_s1 = ~isempty(mix.s1_func) && ~any(isnan(mix.s1_RH));
    has_s2 = ~isempty(mix.s2_func) && ~any(isnan(mix.s2_RH));
    can_add = has_s1 && has_s2;

    aw_add_vec  = [];
    xw_add_vec  = [];
    mfw_add_vec = [];

    if can_add
        % sort pure-salt arrays by a_w ascending for interp1
        v1 = isfinite(aw_s1) & isfinite(xw_s1) & isfinite(mfw_s1);
        v2 = isfinite(aw_s2) & isfinite(xw_s2) & isfinite(mfw_s2);
        [aw1_s, i1s] = sort(aw_s1(v1));
        xw1_s  = xw_s1(v1);   xw1_s  = xw1_s(i1s);
        mfw1_s = mfw_s1(v1);  mfw1_s = mfw1_s(i1s);
        [aw2_s, i2s] = sort(aw_s2(v2));
        xw2_s  = xw_s2(v2);   xw2_s  = xw2_s(i2s);
        mfw2_s = mfw_s2(v2);  mfw2_s = mfw2_s(i2s);

        if numel(aw1_s) > 1 && numel(aw2_s) > 1
            aw_lo = max(min(aw1_s), min(aw2_s));
            aw_hi = min(max(aw1_s), max(aw2_s));

            if aw_hi > aw_lo
                aw_add_vec = linspace(aw_lo, aw_hi, n_pts);

                xw1_at  = interp1(aw1_s, xw1_s,  aw_add_vec, 'linear', NaN);
                xw2_at  = interp1(aw2_s, xw2_s,  aw_add_vec, 'linear', NaN);
                mfw1_at = interp1(aw1_s, mfw1_s, aw_add_vec, 'linear', NaN);
                mfw2_at = interp1(aw2_s, mfw2_s, aw_add_vec, 'linear', NaN);

                % x_w additive: mole-fraction weighted
                xw_add_vec = (1-y_mol)*xw1_at + y_mol*xw2_at;

                % mf_w additive: W-space weighted (mass-conservative)
                W1 = mfw1_at ./ max(1 - mfw1_at, 1e-9);
                W2 = mfw2_at ./ max(1 - mfw2_at, 1e-9);
                W_add = (1-r_mass)*W1 + r_mass*W2;
                mfw_add_vec = W_add ./ (1 + W_add);
            end
        end
    end

    % ---------------------------------------------------------------
    %  Create figure  (3 rows x 2 cols = 6 subplots)
    %  All three curves (salt 1, salt 2, mixture) are overlaid, plus
    %  an "additive" reference line on rows 2 and 3.
    %  Col 1: Constant molality ratio     (y_mol  = 0.5)
    %  Col 2: Constant mass-frac ratio    (r_mass = 0.5)
    %  Row 1: gamma_w  vs  total molality
    %  Row 2: x_w      vs  a_w
    %  Row 3: mf_w     vs  a_w
    % ---------------------------------------------------------------
    fig = figure('Position', [50 50 900 900], 'Color', 'w');

    % colour / line-style palette
    col_s1  = [0.12  0.47  0.71];   % blue   – pure salt 1
    col_s2  = [0.84  0.15  0.16];   % red    – pure salt 2
    col_mix = [0.17  0.63  0.17];   % green  – mixture
    col_add = [0.50  0.50  0.50];   % gray   – additive reference
    lw_pure = 1.8;
    lw_mix  = 2.5;
    lw_add  = 1.8;
    ls_pure = '--';
    ls_mix  = '-';
    ls_add  = '-.';    % dash-dot for additive

    y_full_labels = {'\gamma_w  (water activity coeff.)', ...
                     'x_w  (mole fraction H_2O, ionic)', ...
                     'mf_w  (mass fraction H_2O)'};
    x_label_row1   = 'm_{total}  (mol kg^{-1} H_2O)';
    x_label_rows23 = 'a_w';

    col_titles = { ...
        sprintf('Const. molality ratio  (y_2 = %.1f)', y_mol), ...
        sprintf('Const. mass-frac ratio  (r_2 = %.1f)', r_mass)};

    m_all_vals = [m_s1(:); m_s2(:); m_mx(:); m_mr(:)];
    m_all_vals = m_all_vals(isfinite(m_all_vals));
    m_xmax = max([m_all_vals; 1]) * 1.08;

    % bundle per condition: {s1_data, s2_data, mix_data}
    % additive is the same for both columns (common a_w grid, between the two salts)
    cond_data = { ...
        { {m_s1,aw_s1,xw_s1,mfw_s1,gw_s1}, ...
          {m_s2,aw_s2,xw_s2,mfw_s2,gw_s2}, ...
          {m_mx,aw_mx,xw_mx,mfw_mx,gw_mx} }, ...
        { {m_s1,aw_s1,xw_s1,mfw_s1,gw_s1}, ...
          {m_s2,aw_s2,xw_s2,mfw_s2,gw_s2}, ...
          {m_mr,aw_mr,xw_mr,mfw_mr,gw_mr} }, ...
    };
    curve_colors = {col_s1, col_s2, col_mix};
    curve_lw     = {lw_pure, lw_pure, lw_mix};
    curve_ls     = {ls_pure, ls_pure, ls_mix};
    legend_names = {mix.s1_name, mix.s2_name, 'Mixture', 'Additive (no interact.)'};

    for row = 1:3
        for col = 1:2
            ax = subplot(3, 2, (row-1)*2 + col);
            hold(ax, 'on');

            cd = cond_data{col};

            % --- three main curves (salt 1, salt 2, mixture) ---
            for k = 1:3
                d    = cd{k};
                m_v  = d{1};  aw_v = d{2};
                xw_v = d{3};  mfw_v = d{4};  gw_v = d{5};

                if row == 1
                    plot(ax, m_v, gw_v, curve_ls{k}, ...
                         'Color', curve_colors{k}, 'LineWidth', curve_lw{k}, ...
                         'DisplayName', legend_names{k});
                elseif row == 2
                    plot(ax, aw_v, xw_v, curve_ls{k}, ...
                         'Color', curve_colors{k}, 'LineWidth', curve_lw{k}, ...
                         'DisplayName', legend_names{k});
                else
                    plot(ax, aw_v, mfw_v, curve_ls{k}, ...
                         'Color', curve_colors{k}, 'LineWidth', curve_lw{k}, ...
                         'DisplayName', legend_names{k});
                end
            end

            % --- additive reference line (rows 2 and 3 only) ---
            % Same curve for both columns: blended pure-salt isotherms on
            % a common a_w grid, guaranteed to lie between the two salts.
            if row == 2 && ~isempty(aw_add_vec)
                plot(ax, aw_add_vec, xw_add_vec, ls_add, ...
                     'Color', col_add, 'LineWidth', lw_add, ...
                     'DisplayName', legend_names{4});
            elseif row == 3 && ~isempty(aw_add_vec)
                plot(ax, aw_add_vec, mfw_add_vec, ls_add, ...
                     'Color', col_add, 'LineWidth', lw_add, ...
                     'DisplayName', legend_names{4});
            end

            % axis cosmetics
            if row == 1
                xlabel(ax, x_label_row1,    'FontSize', 9);
                ylabel(ax, y_full_labels{1}, 'FontSize', 9);
                xlim(ax, [0, m_xmax]);
                yline(ax, 1, ':k', 'Alpha', 0.5, 'LineWidth', 1.0, ...
                      'HandleVisibility', 'off');
            elseif row == 2
                xlabel(ax, x_label_rows23,  'FontSize', 9);
                ylabel(ax, y_full_labels{2}, 'FontSize', 9);
                xlim(ax, [0 1]);  ylim(ax, [0 1]);
            else
                xlabel(ax, x_label_rows23,  'FontSize', 9);
                ylabel(ax, y_full_labels{3}, 'FontSize', 9);
                xlim(ax, [0 1]);  ylim(ax, [0 1]);
            end

            grid(ax, 'on');  box(ax, 'on');
            set(ax, 'FontSize', 9);

            if row == 1
                title(ax, col_titles{col}, 'FontSize', 10, 'FontWeight', 'bold');
            end

            % legend: top-right (row 1, col 2) for the 3 salt curves;
            %         middle-right (row 2, col 2) to also show additive
            if row == 1 && col == 2
                legend(ax, 'Location', 'northeast', 'FontSize', 8);
            elseif row == 2 && col == 2 && ~isempty(aw_add_vec)
                legend(ax, 'Location', 'northwest', 'FontSize', 8);
            end
        end
    end

    sgtitle(fig, mix.display, 'FontSize', 13, 'FontWeight', 'bold', ...
            'Interpreter', 'tex');

    % ---- Pitzer mixing parameters (theta, psi) from ion_activity_data ----
    if ~isempty(mix.ion_sys) && isKey(ion_map, mix.ion_sys)
        ip = ion_map(mix.ion_sys);
        if isnan(ip.psi)
            pitzer_str = sprintf('Pitzer (ion\\_activity\\_data):  \\theta = %g   (\\psi not reported)   [%s, I_{max} = %g]', ...
                                 ip.theta, ip.exptl, ip.maxI);
        else
            pitzer_str = sprintf('Pitzer (ion\\_activity\\_data):  \\theta = %g,   \\psi = %g   [%s, I_{max} = %g]', ...
                                 ip.theta, ip.psi, ip.exptl, ip.maxI);
        end
        annotation(fig, 'textbox', [0.05 0.955 0.90 0.025], ...
            'String', pitzer_str, ...
            'FontSize', 9, 'FontWeight', 'normal', ...
            'EdgeColor', [0.7 0.7 0.7], 'BackgroundColor', [0.97 0.97 1.0], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Interpreter', 'tex', 'Margin', 2);
    end

    out_png = fullfile(fig_out_dir, [mix.name '_facet.png']);
    out_fig = fullfile(fig_out_dir, [mix.name '_facet.fig']);
    exportgraphics(fig, out_png, 'Resolution', 150);
    savefig(fig, out_fig);
    close(fig);
end

warning('on', 'all');
fprintf('\nAll facet figures saved to:\n  %s\n', fig_out_dir);


%% ========================================================================
%  Helper functions
%% ========================================================================

function [m_out, aw_out, xw_out, mfw_out, gw_out] = ...
        compute_pure_salt(func_name, MW, nu, needsT, T, RH_range, n_water, n_pts)
% Compute pure-salt thermodynamic curves from a calculate_mf function.

    if isempty(func_name) || any(isnan(RH_range))
        m_out = NaN; aw_out = NaN; xw_out = NaN; mfw_out = NaN; gw_out = NaN;
        return;
    end

    RH_vec = linspace(RH_range(1)*1.001, RH_range(2)*0.999, n_pts);
    m_out   = NaN(1, n_pts);
    aw_out  = NaN(1, n_pts);
    xw_out  = NaN(1, n_pts);
    mfw_out = NaN(1, n_pts);
    gw_out  = NaN(1, n_pts);

    for i = 1:n_pts
        RH = RH_vec(i);
        try
            if needsT
                mf_salt = feval(func_name, RH, T);
            else
                mf_salt = feval(func_name, RH);
            end
            if isnan(mf_salt) || mf_salt <= 0 || mf_salt >= 1
                continue;
            end
            mf_w    = 1 - mf_salt;
            m_salt  = mf_salt * 1000 / (MW * mf_w);
            x_w     = n_water / (n_water + nu * m_salt);
            gamma_w = RH / x_w;
            m_out(i)   = m_salt;
            aw_out(i)  = RH;
            xw_out(i)  = x_w;
            mfw_out(i) = mf_w;
            gw_out(i)  = gamma_w;
        catch
        end
    end
end


function [m_out, aw_out, xw_out, mfw_out, gw_out] = ...
        compute_mix_const_ymol(mix_func, MW1, nu1, MW2, nu2, ...
                               y_mol, mf1_max, mf2_max, n_water, n_pts)
% Mixture curves at constant molality ratio  y_mol = m2/(m1+m2).

    m_test = linspace(0.05, 80, 5000);
    y2 = y_mol;  y1 = 1 - y_mol;
    mf1_t = (y1.*m_test.*MW1) ./ (1000 + y1.*m_test.*MW1 + y2.*m_test.*MW2);
    mf2_t = (y2.*m_test.*MW2) ./ (1000 + y1.*m_test.*MW1 + y2.*m_test.*MW2);

    idx_ok = find(mf1_t <= mf1_max & mf2_t <= mf2_max, 1, 'last');
    if isempty(idx_ok) || idx_ok < 2
        m_out = NaN; aw_out = NaN; xw_out = NaN; mfw_out = NaN; gw_out = NaN;
        return;
    end
    m_max = m_test(idx_ok);
    m_min = max(0.01, m_test(max(find(mf1_t > 0 & mf2_t > 0, 1, 'first'), 1)));
    m_vec = linspace(m_min, m_max, n_pts);

    m_out   = NaN(1,n_pts); aw_out  = NaN(1,n_pts);
    xw_out  = NaN(1,n_pts); mfw_out = NaN(1,n_pts); gw_out  = NaN(1,n_pts);

    for i = 1:n_pts
        m_total = m_vec(i);
        m1 = y1*m_total;  m2 = y2*m_total;
        denom = 1000 + m1*MW1 + m2*MW2;
        mf1 = m1*MW1/denom;  mf2 = m2*MW2/denom;
        mfw = 1 - mf1 - mf2;
        try
            aw = feval(mix_func, mf1, mf2);
            if ~isfinite(aw) || aw <= 0 || aw > 1.02, continue; end
            x_w = n_water / (n_water + nu1*m1 + nu2*m2);
            gw  = aw / x_w;
            m_out(i) = m_total; aw_out(i) = aw;
            xw_out(i) = x_w;   mfw_out(i) = mfw; gw_out(i) = gw;
        catch
        end
    end
end


function [m_out, aw_out, xw_out, mfw_out, gw_out] = ...
        compute_mix_const_rmass(mix_func, MW1, nu1, MW2, nu2, ...
                                r_mass, mf1_max, mf2_max, n_water, n_pts)
% Mixture curves at constant mass-fraction ratio  r_mass = mf2/(mf1+mf2).

    r2 = r_mass;  r1 = 1 - r_mass;
    mft_max = min(mf1_max/r1, mf2_max/r2);
    mft_min = 0.005;
    if mft_max <= mft_min
        m_out = NaN; aw_out = NaN; xw_out = NaN; mfw_out = NaN; gw_out = NaN;
        return;
    end

    mft_vec = linspace(mft_min, mft_max*0.999, n_pts);
    m_out   = NaN(1,n_pts); aw_out  = NaN(1,n_pts);
    xw_out  = NaN(1,n_pts); mfw_out = NaN(1,n_pts); gw_out  = NaN(1,n_pts);

    for i = 1:n_pts
        mft = mft_vec(i);
        mf1 = r1*mft;  mf2 = r2*mft;  mfw = 1 - mft;
        if mfw <= 0, continue; end
        m1 = mf1*1000/(MW1*mfw);  m2 = mf2*1000/(MW2*mfw);
        m_total = m1 + m2;
        try
            aw = feval(mix_func, mf1, mf2);
            if ~isfinite(aw) || aw <= 0 || aw > 1.02, continue; end
            x_w = n_water / (n_water + nu1*m1 + nu2*m2);
            gw  = aw / x_w;
            m_out(i) = m_total; aw_out(i) = aw;
            xw_out(i) = x_w;   mfw_out(i) = mfw; gw_out(i) = gw;
        catch
        end
    end
end
