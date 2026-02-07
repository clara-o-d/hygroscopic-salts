close all 
clear
clc 

% --- Setup Paths ---
[filepath,~,~] = fileparts(mfilename('fullpath'));
% Ensure you have these folders or comment them out if not needed
addpath(fullfile(filepath, '..', 'calculate_mf')); 
addpath(fullfile(filepath, '..', 'util'));         

output_dir = '../figures/uptake';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% --- Constants ---
MWw = 18.015;   % Molecular weight of water
R   = 8.314;    % Gas constant (J / mol K)
T   = 298.15;   % Temperature (K)

% --- Data Input ---
% Cols: {Name, MW, RH_min, RH_max, FuncName, should_plot, Cation_n, Anion_n}
salt_data = {
    % Endothermic salts (should_plot = 0)
    {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1};
    {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1};
    {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1};
    {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1};
    {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1};
    {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1};
    {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1};
    {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
    {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
    {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
    {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1}; 
    {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};    
    {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};      
    {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};      
    {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};    
    {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};    
    {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};      
    
    % Exothermic salts (should_plot = 1)
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 0, 1, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 0, 1, 2};
    {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2};
    {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
    {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1};
    {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2};
    {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2};
    {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
    {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1};
    
    % Sulfates
    {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1}; 
    {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};    
    {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
    {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};    
    {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};    
    {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1}; 
    {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};    
    {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};    
    {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};    
    
    % Nitrates (additional)
    {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2}; 
    {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2}; 
    
    % Halides (additional)
    {'CaBr2', 199.89, 0.6405, 0.9530, 'calculate_mf_CaBr2', 0, 1, 2}; 
    {'CaI2', 293.89, 0.8331, 0.9514, 'calculate_mf_CaI2', 0, 1, 2};    
    {'SrCl2', 158.53, 0.8069, 0.9768, 'calculate_mf_SrCl2', 0, 1, 2}; 
    {'SrBr2', 247.43, 0.7786, 0.9561, 'calculate_mf_SrBr2', 0, 1, 2}; 
    {'SrI2', 341.43, 0.6795, 0.9559, 'calculate_mf_SrI2', 0, 1, 2};    
    {'BaCl2', 208.23, 0.9385, 0.9721, 'calculate_mf_BaCl2', 0, 1, 2}; 
    {'BaBr2', 297.14, 0.8231, 0.9577, 'calculate_mf_BaBr2', 0, 1, 2}; 
    
    % Chlorates
    {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1}; 
};

%% LOAD REFERENCE DATA
baselineFile = fullfile(filepath, '..', 'data', 'baseline_numeric_only.csv');
ref_map_drh = containers.Map();
ref_map_sol = containers.Map();

if isfile(baselineFile)
    opts = detectImportOptions(baselineFile);
    opts.VariableNamingRule = 'preserve'; 
    tbl = readtable(baselineFile, opts);
    varNames = tbl.Properties.VariableNames;
    
    % Identify Columns
    idx_elec = find(strcmpi(varNames, 'electrolyte'), 1);
    idx_sol  = find(ismember(lower(varNames), {'solubility','solubility_limit','solubilitylimit'}), 1);
    idx_drh  = find(ismember(lower(varNames), {'deliquescence_humidity','drh','deliquescencehumidity'}), 1);
    
    if ~isempty(idx_elec)
        raw_names = string(tbl.(varNames{idx_elec}));
        clean_names = erase(raw_names, ["(", ")"]);
        clean_names = strtrim(clean_names);
        
        % Store Solubility
        if ~isempty(idx_sol)
            sol_vals = tbl.(varNames{idx_sol});
            for k = 1:length(clean_names)
                ref_map_sol(upper(clean_names{k})) = sol_vals(k);
            end
        end
        
        % Store DRH
        if ~isempty(idx_drh)
            drh_vals = tbl.(varNames{idx_drh});
            for k = 1:length(clean_names)
                ref_map_drh(upper(clean_names{k})) = drh_vals(k);
            end
        end
    end
else
    warning('Baseline CSV not found.');
end

%% PROCESSING LOOP
results = struct();
for i = 1:length(salt_data)
    row = salt_data{i}; 
    name    = row{1};
    MW_salt = row{2};
    rh_min_calc = row{3}; 
    rh_max  = row{4};
    funcStr = row{5};
    should_plot   = row{6};
    nu      = row{7} + row{8}; 
    
    % --- Determine Best DRH ---
    best_drh = rh_min_calc;
    is_drh_from_csv = false;
    csv_solubility = NaN;
    
    % Try to find in CSV
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
    
    if isKey(ref_map_sol, upper(name))
        csv_solubility = ref_map_sol(upper(name));
    end

    % --- Calculations ---
    lineStyle = '-';
    if should_plot == 1, lineStyle = '--'; end
    
    rh_vec = linspace(rh_min_calc, rh_max, 100);
    calc_func = str2func(funcStr);
    mf_vec = zeros(size(rh_vec));
    for j = 1:length(rh_vec)
        try
            mf_vec(j) = calc_func(rh_vec(j));
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
    results(i).CSVSolubility = csv_solubility;
    results(i).CalcFunc = calc_func;
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
save_fig(fig1, output_dir, ['RH_vs_Uptake_gg']);


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
save_fig(fig2, output_dir, ['RH_vs_Uptake_molmol']);


%% ENTHALPY VS LOG(DRH) PLOT
enthalpyMap = containers.Map();
% Standard Enthalpy of Solution (kJ/mol) at 25C
enthalpyMap('NaCl') = 3.88; enthalpyMap('KCl') = 17.22; enthalpyMap('NH4Cl') = 14.77;
enthalpyMap('CsCl') = 17.78; enthalpyMap('NaNO3') = 20.50; enthalpyMap('AgNO3') = 22.59;
enthalpyMap('KI') = 20.33; enthalpyMap('KNO3') = 34.89; enthalpyMap('NaClO4') = 13.88;
enthalpyMap('KClO3') = 41.38; enthalpyMap('KBr') = 19.87; enthalpyMap('RbCl') = 17.28;
enthalpyMap('K2SO4') = 23.8; enthalpyMap('NH42SO4') = 6.6;

% Grouping
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

% --- Figure 3 ---
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
title('Enthalpy of Solution vs. DRH (Endothermic Salts)');
grid on; set(gcf, 'color', 'w');
leg_handles = []; leg_txt = {};
if ~isempty(group_csv.x), leg_handles(end+1) = h_csv; leg_txt{end+1} = 'Lit. DRH'; end
if ~isempty(group_calc.x), leg_handles(end+1) = h_calc; leg_txt{end+1} = 'Est. DRH'; end
if ~isempty(leg_handles), legend(leg_handles, leg_txt, 'Location', 'best'); end
save_fig(fig3, output_dir, ['Enthalpy_vs_DRH']);

%% VALIDATION: Solubility Comparison
val_csv.x = []; val_csv.y = []; val_csv.lbl = {};
val_calc.x = []; val_calc.y = []; val_calc.lbl = {};

for i = 1:length(results)
    s_limit = results(i).CSVSolubility;
    if ~isnan(s_limit) && s_limit > 0
        calc_drh = results(i).BestDRH;
        func = results(i).CalcFunc;
        try
            mf_at_drh = func(calc_drh);
            u_gg_at_drh = (1 / mf_at_drh) - 1;
            if u_gg_at_drh > 0
                s_calc = (1 / u_gg_at_drh) * 100;
                if results(i).IsDRHFromCSV
                    val_csv.x(end+1) = s_limit; val_csv.y(end+1) = s_calc; val_csv.lbl{end+1} = results(i).LegendName;
                else
                    val_calc.x(end+1) = s_limit; val_calc.y(end+1) = s_calc; val_calc.lbl{end+1} = results(i).LegendName;
                end
            end
        catch
        end
    end
end

% --- Figure 4 ---
fig4 = figure('Position', [250, 250, 700, 700]);
hold on;
h_v_csv = scatter(val_csv.x, val_csv.y, 60, 'filled', ...
    'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k');
h_v_calc = scatter(val_calc.x, val_calc.y, 60, 'filled', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', 'k');

for k = 1:length(val_csv.lbl)
    text(val_csv.x(k), val_csv.y(k), ['  ' val_csv.lbl{k}], 'FontSize', 8, 'Interpreter', 'tex');
end
for k = 1:length(val_calc.lbl)
    text(val_calc.x(k), val_calc.y(k), ['  ' val_calc.lbl{k}], 'FontSize', 8, 'Interpreter', 'tex', 'Color', [0.4 0.4 0.4]);
end
maxVal = max([val_csv.x, val_csv.y, val_calc.x, val_calc.y]) * 1.1;
plot([0, maxVal], [0, maxVal], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Match');
xlabel('Literature Solubility (g salt / 100 g water)');
ylabel('Calculated Solubility from DRH (g salt / 100 g water)');
title('Validation: Calculated vs Literature Solubility');
grid on; set(gcf, 'color', 'w'); axis square; xlim([0, maxVal]); ylim([0, maxVal]);
leg_handles = []; leg_txt = {};
if ~isempty(val_csv.x), leg_handles(end+1) = h_v_csv; leg_txt{end+1} = 'Lit. DRH'; end
if ~isempty(val_calc.x), leg_handles(end+1) = h_v_calc; leg_txt{end+1} = 'Est. DRH'; end
legend(leg_handles, leg_txt, 'Location', 'northwest');
save_fig(fig4, output_dir, ['Solubility_vs_Est_Solubility']);

%% FIGURE 5: Excess Gibbs Energy of Mixing (G_ex) Curves
% Plotting G_ex across the full range of water concentration (xw)
% from xw = 1 (pure water) down to xw_sat (saturation).

fig5 = figure('Position', [300, 300, 900, 700]);
hold on;

for i = 1:length(results)
    drh = results(i).BestDRH;
    func = results(i).CalcFunc;
    MW_salt = results(i).MW;
    name = results(i).LegendName;
    
    % 1. Determine saturation state
    try
        ws_sat = func(drh); 
        if isnan(ws_sat) || ws_sat <= 0 || ws_sat >= 1, continue; end
        
        % Calculate saturation mole fraction of salt
        ns_sat = ws_sat / MW_salt;
        nw_sat = (1 - ws_sat) / MWw;
        xs_sat = ns_sat / (ns_sat + nw_sat);
        
    catch
        continue;
    end
    
    % 2. Create Grid for Salt Mole Fraction (xs)
    % From 0 (pure water) to xs_sat
    xs_grid = linspace(0.01, 0.9*xs_sat, 20);
    
    % Avoid exact 0 for numerical stability during root finding if needed
    xs_grid(1) = 1e-6; 
    
    ln_aw_vec = zeros(size(xs_grid));
    valid_mask = true(size(xs_grid));
    
    % 3. Loop through grid to find ln(aw) at each xs
    for k = 1:length(xs_grid)
        xs_val = xs_grid(k);
        
        % Convert xs -> ws
        % ws = xs * MWs / (xs * MWs + (1-xs) * MWw)
        ws_val = (xs_val * MW_salt) / (xs_val * MW_salt + (1 - xs_val) * MWw);
        
        % Find RH (aw) that gives this ws
        % obj = ws_val - func(rh)
        % Use squared error for robustness
        obj_func_sq = @(rh) (ws_val - func(rh))^2;
        
        try
            % Limit search to valid RH range for the salt.
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
    
    % 4. Integrate Cumulative Trapezoidal
    % G_ex(xs) = RT * int_0^xs (ln(aw)/(1-xs')) dxs'
    if sum(valid_mask) > 10
        xs_valid = xs_grid(valid_mask);
        ln_aw_valid = ln_aw_vec(valid_mask);
        
        integrand = ln_aw_valid ./ (1 - xs_valid);
        cum_integral = cumtrapz(xs_valid, integrand);
        g_ex_curve = R * T * cum_integral;
        
        % 5. Plot vs Water Mole Fraction (xw = 1 - xs)
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
save_fig(fig5, output_dir, ['Mole_Fraction_vs_Gmix_excess']);
% 
% 
% %% FIGURE 6: Water Activity Coefficient vs Water Concentration at Saturation
% % Endothermic Salts Only.
% % Y-Axis: Gamma_w at saturation = a_w / x_w = DRH / x_w
% % X-Axis: x_w at saturation
% 
% fig6 = figure('Position', [350, 350, 800, 600]);
% hold on;
% 
% for i = 1:length(results)
%     % Filter: Endothermic salts only
%     if results(i).should_plot == 1, continue; end
% 
%     drh = results(i).BestDRH;
%     func = results(i).CalcFunc;
%     MW_salt = results(i).MW;
% 
%     try
%         % 1. Calculate Mass Fraction at Saturation (from DRH)
%         ws_sat = func(drh); 
% 
%         if isnan(ws_sat) || ws_sat <= 0 || ws_sat >= 1, continue; end
% 
%         % 2. Convert to Water Mole Fraction (x_w)
%         ns = ws_sat / MW_salt;
%         nw = (1 - ws_sat) / MWw;
% 
%         xw_sat = nw / (ns + nw);
% 
%         % 3. Calculate Activity Coefficient of Water (gamma_w)
%         % a_w = gamma_w * x_w  => gamma_w = a_w / x_w
%         % a_w is simply the DRH fraction.
%         aw_sat = drh;
%         gamma_w_sat = aw_sat / xw_sat;
% 
%         % 4. Plot
%         scatter(log(xw_sat), log(gamma_w_sat), 80, 'filled', ...
%             'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', 'k');
% 
%         text(xw_sat, gamma_w_sat, ['  ' results(i).LegendName], ...
%             'FontSize', 9, 'VerticalAlignment', 'middle');
% 
%     catch
%         continue;
%     end
% end
% 
% xlabel('Mole Fraction of Water at Saturation (x_{w,sat})');
% ylabel('Activity Coefficient of Water at Saturation (\gamma_{w,sat})');
% title('Water Activity Coefficient vs Concentration at Saturation (Endothermic)');
% grid on; set(gcf, 'color', 'w');

disp('All tasks completed.');

function save_fig(h, dir, name)
    set(h, 'Color', 'w');
    fname = fullfile(dir, name);
    print(h, fname, '-dpng', '-r150');
    % close(h);
end