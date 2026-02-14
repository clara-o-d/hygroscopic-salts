% VALIDATE_SOLUBILITY_VS_RH
% Compares literature solubility with calculated solubility from DRH.
% Run independently or after isotherms_screening_all_salts.m.
% Output: figures/uptake/validation_solubility_lit_vs_calc.png

close all
clear
clc

[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

output_dir = fullfile(filepath, '..', 'figures', 'uptake');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Load salt data
raw_salt_data = load_salt_data();
baselineFile = fullfile(filepath, '..', 'data', 'baseline_numeric_only.csv');

ref_map_drh = containers.Map();
ref_map_sol = containers.Map();

if ~isfile(baselineFile)
    error('Baseline CSV not found: %s', baselineFile);
end

opts = detectImportOptions(baselineFile);
opts.VariableNamingRule = 'preserve';
tbl = readtable(baselineFile, opts);
varNames = tbl.Properties.VariableNames;
idx_elec = find(strcmpi(varNames, 'electrolyte_stripped'), 1);
idx_sol  = find(ismember(lower(varNames), {'solubility','solubility_limit','solubilitylimit'}), 1);
idx_drh  = find(ismember(lower(varNames), {'deliquescence_humidity','drh','deliquescencehumidity'}), 1);

if isempty(idx_elec)
    error('electrolyte_stripped column not found in baseline CSV');
end

raw_names = string(tbl.(varNames{idx_elec}));
clean_names = erase(raw_names, ["(", ")"]);
clean_names = strtrim(clean_names);

if ~isempty(idx_sol)
    sol_vals = tbl.(varNames{idx_sol});
    for k = 1:length(clean_names)
        ref_map_sol(upper(clean_names{k})) = sol_vals(k);
    end
end

if ~isempty(idx_drh)
    drh_vals = tbl.(varNames{idx_drh});
    for k = 1:length(clean_names)
        ref_map_drh(upper(clean_names{k})) = drh_vals(k);
    end
end

%% Build results for salts with solubility data
val_csv.x = []; val_csv.y = []; val_csv.lbl = {};
val_calc.x = []; val_calc.y = []; val_calc.lbl = {};

for s = 1:length(raw_salt_data)
    r = raw_salt_data{s};
    name = r{1};
    MW_salt = r{2};
    rh_min = r{3};
    func_name = r{5};
    func_args = r{6};
    nu = r{10};
    T_salt = r{9};

    s_limit = NaN;
    if isKey(ref_map_sol, upper(name))
        s_limit = ref_map_sol(upper(name));
    end
    if isnan(s_limit) || s_limit <= 0
        continue;
    end

    best_drh = rh_min;
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

    try
        if func_args == 1
            mf_at_drh = feval(func_name, best_drh, T_salt);
        else
            mf_at_drh = feval(func_name, best_drh);
        end
        u_gg_at_drh = (1 / mf_at_drh) - 1;
        if u_gg_at_drh > 0
            s_calc = (1 / u_gg_at_drh) * 100;
            cleanName = regexprep(name, '(\d+)', '_$1');
            if is_drh_from_csv
                val_csv.x(end+1) = s_limit; val_csv.y(end+1) = s_calc; val_csv.lbl{end+1} = cleanName;
            else
                val_calc.x(end+1) = s_limit; val_calc.y(end+1) = s_calc; val_calc.lbl{end+1} = cleanName;
            end
        end
    catch
    end
end

%% Plot
fig = figure('Position', [250, 250, 700, 700]);
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

all_vals = [val_csv.x, val_csv.y, val_calc.x, val_calc.y];
if isempty(all_vals)
    warning('No salts with solubility data found. Skipping plot.');
    return;
end
maxVal = max(all_vals) * 1.1;
plot([0, maxVal], [0, maxVal], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Match');

xlabel('Literature Solubility (g salt / 100 g water)');
ylabel('Calculated Solubility from DRH (g salt / 100 g water)');
title('Validation: Calculated vs Literature Solubility');
grid on; set(gcf, 'color', 'w'); axis square; xlim([0, maxVal]); ylim([0, maxVal]);

leg_handles = []; leg_txt = {};
if ~isempty(val_csv.x), leg_handles(end+1) = h_v_csv; leg_txt{end+1} = 'Lit. DRH'; end
if ~isempty(val_calc.x), leg_handles(end+1) = h_v_calc; leg_txt{end+1} = 'Est. DRH'; end
legend(leg_handles, leg_txt, 'Location', 'northwest');

fname = fullfile(output_dir, 'validation_solubility_lit_vs_calc.png');
print(fig, fname, '-dpng', '-r150');
set(fig, 'Color', 'w');

fprintf('Saved: %s\n', fname);
disp('Solubility validation completed.');
