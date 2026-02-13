% PLOT_WATER_ACTIVITY_SUBSET
% Parameterized script for activity coefficient plots by salt subset.
% Usage: plot_water_activity_subset('sulfates')
%        plot_water_activity_subset  % runs endothermic (default)
%        plot_water_activity_subset('all')  % runs all 6 subsets
%
% SALT_FILTER: 'endothermic', 'exothermic', 'sulfates', 'halides',
%              'nitrates', 'chlorates', or 'all'

function plot_water_activity_subset(SALT_FILTER)

if nargin < 1
    SALT_FILTER = 'endothermic';  % default
end

if strcmpi(SALT_FILTER, 'all')
    filters = {'endothermic', 'exothermic', 'sulfates', 'halides', 'nitrates', 'chlorates'};
    for i = 1:length(filters)
        plot_water_activity_subset(filters{i});
    end
    return;
end

close all
clc

% --- Salt subsets (from original subset scripts) ---
SUBSETS = struct();
SUBSETS.endothermic = {'KCl', 'NH4Cl', 'CsCl', 'NaNO3', 'AgNO3', 'KI', 'LiNO3', 'KNO3', ...
                       'NaClO4', 'KClO3', 'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI'};
SUBSETS.exothermic = {'LiCl', 'CaCl2', 'MgCl2', 'LiBr', 'ZnCl2', 'LiI', 'ZnBr2', 'ZnI2', ...
                      'HCl', 'MgNO32', 'LiOH', 'NaOH'};
SUBSETS.exothermic_fit = {'LiCl', 'CaCl2', 'MgCl2', 'LiBr', 'ZnCl2', 'LiI', 'ZnBr2', 'ZnI2', ...
                          'HCl', 'LiOH', 'NaOH'};  % exclude MgNO32 from fit
SUBSETS.sulfates = {'Na2SO4', 'K2SO4', 'NH42SO4', 'MgSO4', 'MnSO4', ...
                    'Li2SO4', 'NiSO4', 'CuSO4', 'ZnSO4'};
SUBSETS.halides = {'NaBr', 'NaI', 'KBr', 'RbCl', 'CsBr', 'CsI', 'CaBr2', 'CaI2', ...
                   'SrCl2', 'SrBr2', 'SrI2', 'BaCl2', 'BaBr2'};
SUBSETS.nitrates = {'NaNO3', 'AgNO3', 'LiNO3', 'NH4NO3', 'KNO3', 'BaNO3', 'CaNO3'};
SUBSETS.chlorates = {'NaClO4', 'LiClO4', 'KClO3'};

% Plot settings per subset
PLOT_OPTS = struct();
PLOT_OPTS.endothermic = struct('xlim_mf', [0.2 1.0], 'ylim_mf', [0.6 1.1], 'xlim_rh', [50 100], 'show_fit', false, 'title_suffix', '');
PLOT_OPTS.exothermic  = struct('xlim_mf', [0.5 1.0], 'ylim_mf', [0 1.2], 'xlim_rh', [0 100], 'show_fit', true, 'title_suffix', '');
PLOT_OPTS.sulfates    = struct('xlim_mf', [0.9 1.0], 'ylim_mf', [0.7 1.1], 'xlim_rh', [80 100], 'show_fit', true, 'title_suffix', 'Sulfate Solutions: ');
PLOT_OPTS.halides     = struct('xlim_mf', [0.6 1.0], 'ylim_mf', [0.6 1.1], 'xlim_rh', [55 100], 'show_fit', true, 'title_suffix', 'Halide Solutions: ');
PLOT_OPTS.nitrates    = struct('xlim_mf', [0.7 1.0], 'ylim_mf', [0.7 1.1], 'xlim_rh', [60 100], 'show_fit', true, 'title_suffix', 'Nitrate Solutions: ');
PLOT_OPTS.chlorates   = struct('xlim_mf', [0.7 1.0], 'ylim_mf', [0.8 1.1], 'xlim_rh', [75 100], 'show_fit', true, 'title_suffix', 'Chlorate Solutions: ');

filter_lower = lower(SALT_FILTER);
if ~isfield(SUBSETS, filter_lower)
    error('Unknown SALT_FILTER: %s. Use: endothermic, exothermic, sulfates, halides, nitrates, chlorates', SALT_FILTER);
end

fit_salts = SUBSETS.(filter_lower);
if strcmp(filter_lower, 'exothermic')
    fit_salts_for_avg = SUBSETS.exothermic_fit;
else
    fit_salts_for_avg = fit_salts;
end

opts = PLOT_OPTS.(filter_lower);

[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));

T = 25;
MWw = 18.015;
num_points = 100;

%% Load salt data and compute RH, mole fraction, gamma
salt_data = load_salt_data();
data = struct();

for k = 1:length(fit_salts)
    salt_name = fit_salts{k};
    idx = find(cellfun(@(r) strcmp(r{1}, salt_name), salt_data), 1);
    if isempty(idx)
        error('Salt %s not found in load_salt_data', salt_name);
    end
    r = salt_data{idx};
    MW_salt = r{2}; RH_min = r{3}; RH_max = r{4}; func_name = r{5}; func_args = r{6}; T_salt = r{9};

    RH_vec = linspace(RH_min, RH_max, num_points);
    mf_salt_vec = zeros(size(RH_vec));
    x_water_vec = zeros(size(RH_vec));

    for i = 1:num_points
        if func_args == 1
            mf_s = feval(func_name, RH_vec(i), T_salt);
        else
            mf_s = feval(func_name, RH_vec(i));
        end
        mf_w = 1 - mf_s;
        x_w = (mf_w / MWw) / ((mf_w / MWw) + (mf_s / MW_salt));
        mf_salt_vec(i) = mf_s;
        x_water_vec(i) = x_w;
    end
    gamma_vec = RH_vec ./ x_water_vec;

    data.(strrep(salt_name, '.', '_')) = struct('RH', RH_vec, 'x_water', x_water_vec, 'gamma', gamma_vec);
end

%% Calculate constrained weighted polynomial fit
all_RH_fit = [];
all_gamma_fit = [];
all_weights = [];

for k = 1:length(fit_salts_for_avg)
    salt_name = fit_salts_for_avg{k};
    fn = strrep(salt_name, '.', '_');
    if ~isfield(data, fn), continue; end
    rh_vec = data.(fn).RH;
    gamma_vec = data.(fn).gamma;
    rh_range = max(rh_vec) - min(rh_vec);
    w_vec = ones(size(rh_vec)) * rh_range;
    all_RH_fit = [all_RH_fit, rh_vec * 100];
    all_gamma_fit = [all_gamma_fit, gamma_vec];
    all_weights = [all_weights, w_vec];
end

Y_vec = all_gamma_fit' - 1;
X_mat = [(all_RH_fit.^2 - 10000)', (all_RH_fit - 100)'];
coeffs = lscov(X_mat, Y_vec, all_weights');
a_fit = coeffs(1); b_fit = coeffs(2); c_fit = 1 - 10000*a_fit - 100*b_fit;

x_fit_line = linspace(opts.xlim_rh(1), 100, 200);
y_fit_line = a_fit * x_fit_line.^2 + b_fit * x_fit_line + c_fit;

%% Exothermic-specific: R^2 and RMSE
if strcmp(filter_lower, 'exothermic')
    gamma_pred = a_fit * all_RH_fit.^2 + b_fit * all_RH_fit + c_fit;
    residuals = all_gamma_fit - gamma_pred;
    RMSE = sqrt(mean(residuals.^2));
    SS_res = sum(residuals.^2);
    SS_tot = sum((all_gamma_fit - mean(all_gamma_fit)).^2);
    R_squared = 1 - SS_res / SS_tot;
end

%% Colors and display names (from original scripts)
colors = lines(length(fit_salts));
disp_names = get_display_names(fit_salts, filter_lower);

%% FIGURE 1: Activity Coefficient vs Mole Fraction
fig1 = figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

for k = 1:length(fit_salts)
    fn = strrep(fit_salts{k}, '.', '_');
    if isfield(data, fn)
        plot(data.(fn).x_water, data.(fn).gamma, 'LineWidth', 2.5, ...
             'DisplayName', disp_names{k}, 'Color', colors(k,:));
    end
end
plot([0.5 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold');
title([opts.title_suffix 'Activity Coefficient vs Mole Fraction'], 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', min(2, ceil(length(fit_salts)/8)));
xlim(opts.xlim_mf); ylim(opts.ylim_mf);
set(gca, 'FontSize', 12); set(gcf, 'color', 'w');

out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient', 'activity_coefficient_subsets');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_name = sprintf('activity_coefficient_gamma_vs_mole_fraction_%s', filter_lower);
print(fig1, fullfile(out_dir, [out_name '.png']), '-dpng', '-r600');

%% FIGURE 2: Activity Coefficient vs Relative Humidity
fig2 = figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

for k = 1:length(fit_salts)
    fn = strrep(fit_salts{k}, '.', '_');
    if isfield(data, fn)
        plot(data.(fn).RH*100, data.(fn).gamma, 'LineWidth', 2.5, ...
             'DisplayName', disp_names{k}, 'Color', colors(k,:));
    end
end
if opts.show_fit
    plot(x_fit_line, y_fit_line, 'k:', 'LineWidth', 3, 'DisplayName', [filter_lower ' Average']);
end
plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');
xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold');
title([opts.title_suffix 'Activity Coefficient vs RH'], 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10, 'NumColumns', min(2, ceil(length(fit_salts)/8)));
xlim(opts.xlim_rh); ylim(opts.ylim_mf);
set(gca, 'FontSize', 12); set(gcf, 'color', 'w');

if strcmp(filter_lower, 'exothermic')
    eqn_str = sprintf(['$\\gamma_w = %.6f\\,RH^2 %+.6f\\,RH %+.6f$\n' ...
                       '$R^2 = %.4f$,   RMSE = %.4f'], a_fit, b_fit, c_fit, R_squared, RMSE);
    annotation('textbox', [0.10 0.70 0.45 0.08], 'String', eqn_str, 'Interpreter', 'latex', ...
        'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1.1);
    savefig(fig2, fullfile(out_dir, 'activity_coefficient_gamma_vs_rh_exothermic.fig'));
    savefig(fig1, fullfile(out_dir, 'activity_coefficient_gamma_vs_mole_fraction_exothermic.fig'));
end

out_name2 = sprintf('activity_coefficient_gamma_vs_rh_%s', filter_lower);
print(fig2, fullfile(out_dir, [out_name2 '.png']), '-dpng', '-r600');

fprintf('Plots saved for subset: %s\n', filter_lower);
end

function disp_names = get_display_names(salts, filter_lower)
% Return TeX display names for legend
disp_names = cell(size(salts));
map = containers.Map(...
    {'Na2SO4','K2SO4','NH42SO4','MgSO4','MnSO4','Li2SO4','NiSO4','CuSO4','ZnSO4', ...
     'CaBr2','CaI2','SrCl2','SrBr2','SrI2','BaCl2','BaBr2', ...
     'NaNO3','AgNO3','LiNO3','NH4NO3','KNO3','BaNO3','CaNO3', ...
     'NaClO4','LiClO4','KClO3', ...
     'NH4Cl','CaCl2','MgCl2','ZnCl2','ZnBr2','ZnI2','MgNO32'}, ...
    {'Na_2SO_4','K_2SO_4','(NH_4)_2SO_4','MgSO_4','MnSO_4','Li_2SO_4','NiSO_4','CuSO_4','ZnSO_4', ...
     'CaBr_2','CaI_2','SrCl_2','SrBr_2','SrI_2','BaCl_2','BaBr_2', ...
     'NaNO_3','AgNO_3','LiNO_3','NH_4NO_3','KNO_3','Ba(NO_3)_2','Ca(NO_3)_2', ...
     'NaClO_4','LiClO_4','KClO_3', ...
     'NH_4Cl','CaCl_2','MgCl_2','ZnCl_2','ZnBr_2','ZnI_2','Mg(NO_3)_2'});
for i = 1:length(salts)
    if isKey(map, salts{i})
        disp_names{i} = map(salts{i});
    else
        disp_names{i} = salts{i};
    end
end
end
