% ISOTHERMS_SCREENING_ALL_SALTS
% Consolidated script for water uptake isotherms across all salts.
% Uses load_salt_data() as single source of truth. Replaces the former
% Isotherms_screening_endothermic.m and Isotherms_screening_exothermic.m.
%
% For solubility-humidity validation, run validate_solubility_vs_rh.m separately.

close all
clear
clc

% --- Setup Paths ---
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

output_dir = fullfile(filepath, '..', 'figures', 'uptake');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% --- Constants ---
MWw = 18.015;   % Molecular weight of water (g/mol)
R   = 8.314;    % Gas constant (J / mol K)
T_K = 298.15;   % Temperature (K)

% --- Data Input (from canonical load_salt_data) ---
raw_salt_data = load_salt_data();
salt_data = cell(size(raw_salt_data));
for i = 1:length(raw_salt_data)
    r = raw_salt_data{i};
    salt_data{i} = {r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, 0, r{12}, r{13}, r{9}};
    % {Name, MW, RH_min, RH_max, FuncName, func_args, should_plot, n_cat, n_an, T_salt}
end

%% LOAD REFERENCE DATA (DRH from baseline for plotting)
baselineFile = fullfile(filepath, '..', 'data', 'baseline_numeric_only.csv');
ref_map_drh = containers.Map();

if isfile(baselineFile)
    opts = detectImportOptions(baselineFile);
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(baselineFile, opts);
    varNames = tbl.Properties.VariableNames;
    idx_elec = find(strcmpi(varNames, 'electrolyte_stripped'), 1);
    idx_drh  = find(ismember(lower(varNames), {'deliquescence_humidity','drh','deliquescencehumidity'}), 1);

    if ~isempty(idx_elec) && ~isempty(idx_drh)
        raw_names = string(tbl.(varNames{idx_elec}));
        clean_names = erase(raw_names, ["(", ")"]);
        clean_names = strtrim(clean_names);
        drh_vals = tbl.(varNames{idx_drh});
        for k = 1:length(clean_names)
            ref_map_drh(upper(clean_names{k})) = drh_vals(k);
        end
    end
end

%% PROCESSING LOOP
results = struct();
for i = 1:length(salt_data)
    row = salt_data{i};
    name       = row{1};
    MW_salt    = row{2};
    rh_min_calc = row{3};
    rh_max     = row{4};
    funcStr    = row{5};
    func_args  = row{6};
    should_plot = row{7};
    nu         = row{8} + row{9};
    T_salt     = row{10};

    % --- Determine Best DRH ---
    best_drh = rh_min_calc;
    is_drh_from_csv = false;

    if isKey(ref_map_drh, upper(name))
        csv_val = ref_map_drh(upper(name));
        if ~isnan(csv_val) && csv_val > 0
            if csv_val > 1.0
                csv_val = csv_val / 100.0;
            end
            best_drh = csv_val;
            is_drh_from_csv = true;
        end
    end

    % --- Calculations (with temperature support for func_args == 1) ---
    lineStyle = '-';
    if should_plot == 1, lineStyle = '--'; end

    calc_func = str2func(funcStr);
    rh_vec = linspace(rh_min_calc, rh_max, 100);
    mf_vec = zeros(size(rh_vec));

    for j = 1:length(rh_vec)
        try
            if func_args == 1
                mf_vec(j) = calc_func(rh_vec(j), T_salt);
            else
                mf_vec(j) = calc_func(rh_vec(j));
            end
        catch
            mf_vec(j) = NaN;
        end
    end

    u_gg = (1 ./ mf_vec) - 1;
    u_mm = u_gg * (MW_salt / (nu * MWw));
    cleanName = regexprep(name, '(\d+)', '_$1');

    results(i).Name = name;
    results(i).LegendName = cleanName;
    results(i).MW = MW_salt;
    results(i).Nu = nu;
    results(i).RH = rh_vec;
    results(i).U_gg = u_gg;
    results(i).U_mm = u_mm;
    results(i).LineStyle = lineStyle;
    results(i).should_plot = should_plot;
    results(i).BestDRH = best_drh;
    results(i).IsDRHFromCSV = is_drh_from_csv;
    results(i).CalcFunc = calc_func;
    results(i).FuncArgs = func_args;
    results(i).T_salt = T_salt;

    RH_target = 0.90;
    try
        if func_args == 1
            mf_90 = calc_func(RH_target, T_salt);
        else
            mf_90 = calc_func(RH_target);
        end
        u_gg_90 = (1 / mf_90) - 1;
        u_mm_90 = u_gg_90 * (MW_salt / (nu * MWw));
        fprintf('%-20s | U_gg(90%%) = %8.4f g/g | U_mm(90%%) = %8.4f mol/mol\n', ...
            name, u_gg_90, u_mm_90);
    catch
        fprintf('%-20s | RH=90%% evaluation failed\n', name);
    end
end

%% PLOTTING 1 & 2 (Standard Uptake Curves)
colors = lines(length(results));

% --- Figure 1: g/g Uptake ---
fig1 = figure('Position', [100, 100, 900, 700]);
hold on
for i = 1:length(results)
    plot(results(i).RH, results(i).U_gg, 'LineWidth', 1.5, ...
         'LineStyle', results(i).LineStyle, 'Color', colors(i,:));
    text(results(i).RH(1), results(i).U_gg(1), ['  ' results(i).LegendName], ...
        'Color', colors(i,:), 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'tex');
end
xlabel('Relative Humidity (RH)'); ylabel('Uptake (g/g)');
title('Salt Water Uptake (Mass Basis)'); grid on; set(gcf,'color','w'); xlim([0, 1]); ylim([0, 15]);
save_fig(fig1, output_dir, 'uptake_rh_vs_gg');

% --- Figure 2: mol/mol Uptake ---
fig2 = figure('Position', [150, 150, 900, 700]);
hold on
for i = 1:length(results)
    plot(results(i).RH * 100, results(i).U_mm, 'LineWidth', 1.5, ...
         'LineStyle', results(i).LineStyle, 'Color', colors(i,:));
    text(results(i).RH(1)*100, results(i).U_mm(1), ['  ' results(i).LegendName], ...
        'Color', colors(i,:), 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'tex');
end
xlabel('Relative Humidity (%)'); ylabel('Uptake (mol/mol)');
title('Salt Water Uptake (Molar Basis)'); grid on; set(gcf,'color','w'); xlim([0, 100]); ylim([0, 15]);
save_fig(fig2, output_dir, 'uptake_rh_vs_molmol');

%% ENTHALPY VS LOG(DRH) PLOT
enthalpyMap = containers.Map();
enthalpyMap('NaCl') = 3.88; enthalpyMap('KCl') = 17.22; enthalpyMap('NH4Cl') = 14.77;
enthalpyMap('CsCl') = 17.78; enthalpyMap('NaNO3') = 20.50; enthalpyMap('AgNO3') = 22.59;
enthalpyMap('KI') = 20.33; enthalpyMap('KNO3') = 34.89; enthalpyMap('NaClO4') = 13.88;
enthalpyMap('KClO3') = 41.38; enthalpyMap('KBr') = 19.87; enthalpyMap('RbCl') = 17.28;
enthalpyMap('K2SO4') = 23.8; enthalpyMap('NH42SO4') = 6.6;

group_csv.x = []; group_csv.y = []; group_csv.lbl = {};
group_calc.x = []; group_calc.y = []; group_calc.lbl = {};

for i = 1:length(results)
    sName = results(i).Name;
    if isKey(enthalpyMap, sName)
        dH = enthalpyMap(sName);
        drh_val = results(i).BestDRH;
        x_plot = log(drh_val);
        y_plot = dH;
        if results(i).IsDRHFromCSV
            group_csv.x(end+1) = x_plot; group_csv.y(end+1) = y_plot; group_csv.lbl{end+1} = results(i).LegendName;
        else
            group_calc.x(end+1) = x_plot; group_calc.y(end+1) = y_plot; group_calc.lbl{end+1} = results(i).LegendName;
        end
    end
end

fig3 = figure('Position', [200, 200, 800, 600]);
hold on;
h_csv = scatter(group_csv.x, group_csv.y, 80, 'filled', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', 'k');
h_calc = scatter(group_calc.x, group_calc.y, 80, 'filled', ...
    'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'k');
for k = 1:length(group_csv.lbl)
    text(group_csv.x(k), group_csv.y(k), ['  ' group_csv.lbl{k}], 'FontSize', 9, 'VerticalAlignment', 'middle');
end
for k = 1:length(group_calc.lbl)
    text(group_calc.x(k), group_calc.y(k), ['  ' group_calc.lbl{k}], 'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', [0.4 0.4 0.4]);
end
xlabel('Natural Log of Deliquescence RH (ln(DRH))');
ylabel('Standard Enthalpy of Solution (\DeltaH_{soln}^{\circ}, kJ/mol)');
title('Enthalpy of Solution vs. DRH');
grid on; set(gcf, 'color', 'w');
leg_handles = []; leg_txt = {};
if ~isempty(group_csv.x), leg_handles(end+1) = h_csv; leg_txt{end+1} = 'Lit. DRH'; end
if ~isempty(group_calc.x), leg_handles(end+1) = h_calc; leg_txt{end+1} = 'Est. DRH'; end
if ~isempty(leg_handles), legend(leg_handles, leg_txt, 'Location', 'best'); end
save_fig(fig3, output_dir, 'uptake_enthalpy_vs_drh');

%% FIGURE 4: Excess Gibbs Energy of Mixing (G_ex) Curves
fig4 = figure('Position', [300, 300, 900, 700]);
hold on;

for i = 1:length(results)
    drh = results(i).BestDRH;
    func = results(i).CalcFunc;
    func_args_i = results(i).FuncArgs;
    T_salt_i = results(i).T_salt;
    MW_salt = results(i).MW;
    name = results(i).LegendName;

    eval_mf = @(rh) eval_mf_safe(func, rh, func_args_i, T_salt_i);

    try
        ws_sat = eval_mf(drh);
        if isnan(ws_sat) || ws_sat <= 0 || ws_sat >= 1, continue; end

        ns_sat = ws_sat / MW_salt;
        nw_sat = (1 - ws_sat) / MWw;
        xs_sat = ns_sat / (ns_sat + nw_sat);
    catch
        continue;
    end

    xs_grid = linspace(0.01, 0.9*xs_sat, 20);
    xs_grid(1) = 1e-6;
    ln_aw_vec = zeros(size(xs_grid));
    valid_mask = true(size(xs_grid));

    for k = 1:length(xs_grid)
        xs_val = xs_grid(k);
        ws_val = (xs_val * MW_salt) / (xs_val * MW_salt + (1 - xs_val) * MWw);
        obj_func_sq = @(rh) (ws_val - eval_mf(rh))^2;
        try
            aw_val = fminbnd(obj_func_sq, 1e-5, 0.9999);
            if aw_val <= 0
                valid_mask(k) = false;
            else
                ln_aw_vec(k) = log(aw_val);
            end
        catch
            valid_mask(k) = false;
        end
    end

    if sum(valid_mask) > 10
        xs_valid = xs_grid(valid_mask);
        ln_aw_valid = ln_aw_vec(valid_mask);
        integrand = ln_aw_valid ./ (1 - xs_valid);
        cum_integral = cumtrapz(xs_valid, integrand);
        g_ex_curve = R * T_K * cum_integral;
        xw_plot = 1 - xs_valid;
        plot(xw_plot, g_ex_curve, 'LineWidth', 1.5, ...
             'LineStyle', results(i).LineStyle, 'Color', colors(i,:));
    end
end

xlabel('Mole Fraction of Water (x_w)');
ylabel('Excess Gibbs Energy of Mixing (G^{ex}, J/mol)');
title('Excess Gibbs Energy of Mixing vs. Water Concentration');
grid on; set(gcf, 'color', 'w');
xlim([0.8, 1.0]);
h_endo = plot(NaN,NaN,'k-','LineWidth',1.5,'DisplayName','Endothermic');
h_exo = plot(NaN,NaN,'k--','LineWidth',1.5,'DisplayName','Exothermic');
legend([h_endo, h_exo], 'Location', 'southwest');
save_fig(fig4, output_dir, 'uptake_gmix_excess_vs_mole_fraction');

disp('Isotherms screening completed. Run validate_solubility_vs_rh.m for solubility validation.');

%% Helper functions
function mf = eval_mf_safe(func, rh, func_args, T_salt)
    try
        if func_args == 1
            mf = func(rh, T_salt);
        else
            mf = func(rh);
        end
    catch
        mf = NaN;
    end
end

function save_fig(h, dir, name)
    set(h, 'Color', 'w');
    fname = fullfile(dir, [name '.png']);
    print(h, fname, '-dpng', '-r150');
end
