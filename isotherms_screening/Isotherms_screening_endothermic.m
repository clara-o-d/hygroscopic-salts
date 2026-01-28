close all 
clear
clc 

% --- Setup Paths ---
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf')); 
addpath(fullfile(filepath, '..', 'util'));         

% --- Constants ---
MWw = 18.015; % Molecular weight of water

% --- Data Input ---
% Cols: {Name, MW, RH_min, RH_max, FuncName, IsExo, Cation_n, Anion_n}
salt_data = {
    % Endothermic salts (IsExo = 0)
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
    
    % Exothermic salts (IsExo = 1)
    {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1};
    {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
    {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
    {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
    {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2};
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

% --- Processing Loop ---
results = struct();

for i = 1:length(salt_data)
    % 1. Extract the row
    row = salt_data{i}; 
    
    % 2. Extract parameters
    name    = row{1};
    MW_salt = row{2};
    rh_min  = row{3};
    rh_max  = row{4};
    funcStr = row{5};
    isExo   = row{6}; % The 6th column: 1 if Exothermic, 0 if Endothermic
    nu      = row{7} + row{8}; 
    
    % Determine Line Style
    if isExo == 1
        lineStyle = '--';
    else
        lineStyle = '-';
    end

    % Generate RH vector
    rh_vec = linspace(rh_min, rh_max, 100);
    
    % Get function handle and calculate
    calc_func = str2func(funcStr);
    
    % Safe calculation loop
    mf_vec = zeros(size(rh_vec));
    for j = 1:length(rh_vec)
        try
            mf_vec(j) = calc_func(rh_vec(j));
        catch
            mf_vec(j) = NaN; 
        end
    end
    
    % Calculate Uptake (g/g)
    % U_gg = m_water / m_salt = (1/mf) - 1
    u_gg = (1 ./ mf_vec) - 1;
    
    % Calculate Uptake (mol/mol)
    u_mm = u_gg * (MW_salt / (nu * MWw));
    
    % Format name for Legend
    cleanName = regexprep(name, '(\d+)', '_$1');
    
    % Store in struct
    results(i).Name = name;
    results(i).LegendName = cleanName;
    results(i).MW = MW_salt; 
    results(i).RH = rh_vec;
    results(i).U_gg = u_gg;
    results(i).U_mm = u_mm;
    results(i).LineStyle = lineStyle;
    results(i).IsExo = isExo;
end

%% PLOTTING

% Define a set of colors to cycle through
colors = lines(length(results));

% --- Figure 1: g/g Uptake ---
fig1 = figure('Position', [100, 100, 900, 700]);
hold on
for i = 1:length(results)
    % Plot with specific LineStyle
    p = plot(results(i).RH, results(i).U_gg, 'LineWidth', 1.5, ...
         'LineStyle', results(i).LineStyle, 'Color', colors(i,:));
     
    % Label at start of curve (approx Deliquescence Point)
    x_label = results(i).RH(1);
    y_label = results(i).U_gg(1);
    
    text(x_label, y_label, ['  ' results(i).LegendName], ...
        'Color', colors(i,:), ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', ...
        'Interpreter', 'tex');
end

xlabel('Relative Humidity (RH)');
ylabel('Uptake (g/g)');
title('Salt Water Uptake (Mass Basis)');
grid on; set(gcf,'color','w');
xlim([0, 1]);

% Legend for Line Styles only
h1 = plot(NaN,NaN,'k-','LineWidth',1.5,'DisplayName','Endothermic');
h2 = plot(NaN,NaN,'k--','LineWidth',1.5,'DisplayName','Exothermic');
legend([h1, h2], 'Location', 'northwest');

print(fullfile(filepath, '..', 'figures', 'Uptake_All_gg'),'-dpng','-r600');

% --- Figure 2: mol/mol Uptake ---
fig2 = figure('Position', [150, 150, 900, 700]);
hold on
for i = 1:length(results)
    plot(results(i).RH * 100, results(i).U_mm, 'LineWidth', 1.5, ...
         'LineStyle', results(i).LineStyle, 'Color', colors(i,:));
     
    % Label at start
    x_label = results(i).RH(1) * 100;
    y_label = results(i).U_mm(1);
    
    text(x_label, y_label, ['  ' results(i).LegendName], ...
        'Color', colors(i,:), ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'bottom', ...
        'Interpreter', 'tex');
end

xlabel('Relative Humidity (%)');
ylabel('Uptake (mol water / mol dissociation particle)');
title('Salt Water Uptake (Molar Basis)');
grid on; set(gcf,'color','w');
xlim([0, 100]);
legend([h1, h2], 'Location', 'northwest');

print(fullfile(filepath, '..', 'figures', 'Uptake_All_molmol'),'-dpng','-r600');


%% VALIDATION: Solubility Comparison
baselineFile = fullfile(filepath, '..', 'data', 'baseline_numeric_only.csv');
fprintf('Processing validation against: %s\n', baselineFile);

if isfile(baselineFile)
    opts = detectImportOptions(baselineFile);
    opts.VariableNamingRule = 'preserve'; 
    tbl = readtable(baselineFile, opts);
    varNames = tbl.Properties.VariableNames;
    
    % 1. Identify "electrolyte" column
    elec_col_name = '';
    % Case insensitive search for 'electrolyte'
    col_idx = find(strcmpi(varNames, 'electrolyte'));
    
    if ~isempty(col_idx)
        elec_col_name = varNames{col_idx(1)};
        fprintf('Found electrolyte column: "%s"\n', elec_col_name);
        
        % 2. Identify "solubility_limit" column
        sol_col_name = '';
        possibleSolNames = {'solubility_limit', 'solubility', 'Solubility', 'SolubilityLimit'};
        for k = 1:length(possibleSolNames)
            if ismember(possibleSolNames{k}, varNames)
                sol_col_name = possibleSolNames{k};
                break;
            end
        end
        
        if ~isempty(sol_col_name)
             fprintf('Found solubility column: "%s"\n', sol_col_name);
             
             % Get CSV Data
             raw_names = string(tbl.(elec_col_name));
             sol_values = tbl.(sol_col_name);
             
             % Clean CSV names: Remove parentheses ()
             clean_csv_names = erase(raw_names, ["(", ")"]);
             % Trim whitespace
             clean_csv_names = strtrim(clean_csv_names);
             
             sol_ref_vec = [];
             sol_calc_vec = [];
             labels = {};
             
             % Loop through our salt data
             for i = 1:length(results)
                 target_name = results(i).Name;
                 
                 % Find match in cleaned CSV names (Case Insensitive)
                 idx = find(strcmpi(clean_csv_names, target_name), 1);
                 
                 if ~isempty(idx)
                     s_limit = sol_values(idx);
                     
                     if s_limit > 0 && ~isnan(s_limit)
                         u_calc_drh = results(i).U_gg(1);
                         
                         if u_calc_drh > 0
                             s_calc = (1 / u_calc_drh) * 100;
                             
                             sol_ref_vec(end+1) = s_limit;
                             sol_calc_vec(end+1) = s_calc;
                             labels{end+1} = results(i).LegendName;
                         else
                             fprintf('Skipping %s: Calculated uptake is zero/invalid.\n', target_name);
                         end
                     else
                         fprintf('Skipping %s: CSV Solubility is %f (invalid).\n', target_name, s_limit);
                     end
                 else
                     fprintf('No match found in CSV for: %s (Looked for "%s")\n', target_name, target_name);
                 end
             end
             
             % --- Plot Validation ---
             if ~isempty(sol_ref_vec)
                fig3 = figure('Position', [200, 200, 700, 700]);
                
                scatter(sol_ref_vec, sol_calc_vec, 60, 'filled', ...
                    'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k');
                hold on;
                
                % Add Labels to Scatter Points
                for k = 1:length(labels)
                    text(sol_ref_vec(k), sol_calc_vec(k), ['  ' labels{k}], ...
                        'FontSize', 8, 'Interpreter', 'tex', ...
                        'VerticalAlignment', 'middle');
                end
                
                % 1:1 Line
                maxVal = max([sol_ref_vec, sol_calc_vec]) * 1.1;
                plot([0, maxVal], [0, maxVal], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Match');
                
                xlabel('Literature Solubility (g salt / 100 g water)');
                ylabel('Calculated Solubility from DRH (g salt / 100 g water)');
                title('Validation: Calculated vs Literature Solubility');
                grid on; set(gcf, 'color', 'w');
                axis square;
                xlim([0, maxVal]);
                ylim([0, maxVal]);
                legend('Data', '1:1 Line', 'Location', 'northwest');
                
                print(fullfile(filepath, '..', 'figures', 'Solubility_Validation'), '-dpng', '-r600');
                disp('Solubility validation plot generated successfully.');
             else
                 warning('No matching valid data points found for validation plot.');
             end
             
        else
            warning('Could not find solubility column. Checked: %s', strjoin(possibleSolNames, ', '));
        end
    else
        warning('Could not find "electrolyte" column in CSV.');
    end
else
    warning('Baseline CSV file not found at: %s', baselineFile);
end

disp('All tasks completed.');