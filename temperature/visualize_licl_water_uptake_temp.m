close all
clear
clc

% Script: visualize_licl_water_uptake_temp.m
%
% Visualizes the molar water uptake  x_w  and the water activity coefficient
% gamma_w for LiCl solutions across ambient RH and temperature.
%
% Physical basis
% ──────────────
%   aw = gamma_w * x_w
%
%   where  x_w = n_water / (n_water + n_ions)
%              = 55.51   / (55.51   + 2*m_eq)          [LiCl → Li⁺ + Cl⁻, ν=2]
%
%   The solution equilibrates to the molality where
%       calculate_activity_temperature_LiCl(mf, T) = RH_ambient
%   which is solved numerically with fzero.  The resulting m_eq is used to
%   compute both x_w and  gamma_w = RH / x_w.
%
% Plots produced
% ──────────────
%   1.  3D surface  – x_w vs (RH, T)
%   2.  Contour map – x_w vs (RH, T)
%   3.  2D lines    – x_w vs T at fixed ambient RH   (blue → red, low→high RH)
%   4.  2D lines    – x_w vs RH at fixed T            (blue → red, cold → hot)
%   5.  3D surface  – gamma_w vs (RH, T)
%   6.  2D lines    – gamma_w vs RH at fixed T         (blue → red)
%   7.  2D lines    – gamma_w vs T at fixed ambient RH (blue → red, low→high RH)

%% ── Paths ───────────────────────────────────────────────────────────────────
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_temperature'));
addpath(fullfile(filepath, '..', 'util'));

fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end

%% ── Constants ───────────────────────────────────────────────────────────────
MW_LiCl = 42.394;   % g/mol
MWw     = 18.015;   % g/mol
n_w0    = 1000 / MWw;   % mol water per kg H₂O  ≈ 55.51

%% ── Calibration limits ──────────────────────────────────────────────────────
mf_min_cal = 0.041481;
mf_max_cal = 0.440682;
T_min_cal  = 25.0;
T_max_cal  = 100.0;

%% ── Define grids ─────────────────────────────────────────────────────────────
RH_vec = linspace(0.15, 0.99, 60);
T_vec  = linspace(25,  100,   60);

[RH_grid, T_grid] = meshgrid(RH_vec, T_vec);

fprintf('Solving for equilibrium LiCl concentration (x_w analysis)...\n');
fprintf('RH range:          %.2f to %.2f\n', min(RH_vec), max(RH_vec));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(T_vec), max(T_vec));
fprintf('Grid size:         %d × %d = %d points\n', length(RH_vec), length(T_vec), numel(RH_grid));

%% ── Invert activity polynomial ───────────────────────────────────────────────
mf_eq_grid = nan(size(RH_grid));
m_eq_grid  = nan(size(RH_grid));

fzero_opts = optimset('Display', 'off', 'TolX', 1e-7);

for i = 1:numel(RH_grid)
    RH_t = RH_grid(i);
    T    = T_grid(i);

    if T < T_min_cal || T > T_max_cal, continue; end

    residual = @(mf) calculate_activity_temperature_LiCl(mf, T) - RH_t;

    try
        r_lo = residual(mf_min_cal);
        r_hi = residual(mf_max_cal);

        if sign(r_lo) ~= sign(r_hi)
            mf_sol = fzero(residual, [mf_min_cal, mf_max_cal], fzero_opts);
            if mf_sol >= mf_min_cal && mf_sol <= mf_max_cal
                mf_eq_grid(i) = mf_sol;
                m_eq_grid(i)  = 1000 * mf_sol / (MW_LiCl * (1 - mf_sol));
            end
        end
    catch
    end
end

n_valid = sum(~isnan(m_eq_grid(:)));
fprintf('Converged: %d / %d points\n\n', n_valid, numel(RH_grid));

%% ── Derived quantities ───────────────────────────────────────────────────────
%   x_w = n_w / (n_w + 2*m_eq)       [LiCl → 2 ions]
%   gamma_w = aw / x_w  =  RH / x_w
x_w_grid    = n_w0 ./ (n_w0 + 2.*m_eq_grid);
gamma_w_grid = RH_grid ./ x_w_grid;

% Mask unphysical gamma values (NaN or gamma > 1.5)
gamma_w_grid(gamma_w_grid > 1.5 | gamma_w_grid <= 0) = NaN;

fprintf('x_w range:     %.4f to %.4f\n', min(x_w_grid(:),[],'omitnan'), max(x_w_grid(:),[],'omitnan'));
fprintf('gamma_w range: %.4f to %.4f\n', min(gamma_w_grid(:),[],'omitnan'), max(gamma_w_grid(:),[],'omitnan'));

%% ── Plot 1: 3D surface – x_w vs (RH, T) ─────────────────────────────────────
figure('Position', [100 100 1200 900]);
surf(RH_grid*100, T_grid, x_w_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);

xlabel('Ambient RH (%)',   'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Molar Water Uptake  x_w', 'FontSize', 16, 'FontWeight', 'bold');
title('LiCl Equilibrium Water Uptake  x_w  vs Ambient RH and Temperature', ...
      'FontSize', 18, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = 'x_w';
cb.Label.FontSize = 14; cb.Label.FontWeight = 'bold';
colormap(parula);
view(-45, 30);
grid on; box on;
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_xw_3D_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_xw_3D_RH_temp.fig'));
fprintf('Saved: LiCl_xw_3D_RH_temp.png\n');

%% ── Plot 2: Contour – x_w vs (RH, T) ────────────────────────────────────────
figure('Position', [100 100 1200 900]);
[C, h] = contourf(RH_grid*100, T_grid, x_w_grid, 20);
clabel(C, h, 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');

xlabel('Ambient RH (%)',   'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
title('LiCl Molar Water Uptake  x_w  Contours', 'FontSize', 20, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = 'Molar Water Uptake  x_w';
cb.Label.FontSize = 14; cb.Label.FontWeight = 'bold';
colormap(parula);
grid on; box on;
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_xw_contour_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_xw_contour_RH_temp.fig'));
fprintf('Saved: LiCl_xw_contour_RH_temp.png\n');

%% ── Plot 3: x_w vs T at fixed ambient RH ─────────────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

RH_fixed = [0.20, 0.35, 0.50, 0.65, 0.80, 0.95];
n_RH     = length(RH_fixed);

% Dry/low-RH → blue  ;  humid/high-RH → red
colors_RH = interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], (1:n_RH)');
markers   = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_RH
    [~, idx] = min(abs(RH_vec - RH_fixed(i)));
    xw   = x_w_grid(:, idx);        % vary T
    valid = ~isnan(xw);
    if ~any(valid), continue; end

    mk_idx   = round(linspace(1, sum(valid), 8));
    t_valid  = T_vec(valid);
    xw_valid = xw(valid);

    plot(t_valid, xw_valid, 'LineWidth', 2.5, 'Color', colors_RH(i,:), ...
         'DisplayName', sprintf('RH = %.0f%%', RH_fixed(i)*100));

    mk_idx = min(mk_idx, length(t_valid));
    plot(t_valid(mk_idx), xw_valid(mk_idx), markers{i}, 'MarkerSize', 7, ...
         'Color', colors_RH(i,:), 'MarkerFaceColor', colors_RH(i,:), ...
         'HandleVisibility', 'off');
end

colormap(gca, interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_RH,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks      = linspace(0, 1, n_RH);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f%%', x*100), RH_fixed, 'UniformOutput', false);
cb.Label.String = 'Ambient RH'; cb.Label.FontSize = 13; cb.Label.FontWeight = 'bold';

xlabel('Temperature (°C)',         'FontSize', 16, 'FontWeight', 'bold');
ylabel('Molar Water Uptake  x_w',  'FontSize', 16, 'FontWeight', 'bold');
title({'LiCl Equilibrium Water Uptake vs Temperature at Fixed Ambient RH', ...
       '(concentration adjusts to equilibrium at each T)'}, ...
      'FontSize', 17, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_xw_vs_T_fixed_RH.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_xw_vs_T_fixed_RH.fig'));
fprintf('Saved: LiCl_xw_vs_T_fixed_RH.png\n');

%% ── Plot 4: x_w vs RH at fixed temperatures ──────────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

T_fixed = [25, 40, 55, 70, 85, 100];
n_T     = length(T_fixed);

% Cold → blue,  hot → red
colors_T = interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], (1:n_T)');

for i = 1:n_T
    [~, idx] = min(abs(T_vec - T_fixed(i)));
    xw   = x_w_grid(idx, :)';        % vary RH
    valid = ~isnan(xw);
    if ~any(valid), continue; end

    rh_valid = RH_vec(valid)*100;
    xw_valid = xw(valid);

    mk_idx = round(linspace(1, sum(valid), 8));
    mk_idx = min(mk_idx, length(rh_valid));

    plot(rh_valid, xw_valid, 'LineWidth', 2.5, 'Color', colors_T(i,:), ...
         'DisplayName', sprintf('T = %.0f°C', T_fixed(i)));

    plot(rh_valid(mk_idx), xw_valid(mk_idx), markers{i}, 'MarkerSize', 7, ...
         'Color', colors_T(i,:), 'MarkerFaceColor', colors_T(i,:), ...
         'HandleVisibility', 'off');
end

colormap(gca, interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_T,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks      = linspace(0, 1, n_T);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f°C', x), T_fixed, 'UniformOutput', false);
cb.Label.String = 'Temperature'; cb.Label.FontSize = 13; cb.Label.FontWeight = 'bold';

xlabel('Ambient RH (%)',           'FontSize', 16, 'FontWeight', 'bold');
ylabel('Molar Water Uptake  x_w',  'FontSize', 16, 'FontWeight', 'bold');
title({'LiCl Equilibrium Water Uptake vs Ambient RH at Fixed Temperatures', ...
       '(concentration adjusts to equilibrium at each condition)'}, ...
      'FontSize', 17, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_xw_vs_RH_fixed_T.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_xw_vs_RH_fixed_T.fig'));
fprintf('Saved: LiCl_xw_vs_RH_fixed_T.png\n');

%% ── Plot 5: 3D surface – gamma_w vs (RH, T) ──────────────────────────────────
figure('Position', [100 100 1200 900]);
surf(RH_grid*100, T_grid, gamma_w_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);

xlabel('Ambient RH (%)',   'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('\gamma_w  =  a_w / x_w', 'FontSize', 16, 'FontWeight', 'bold');
title('LiCl Water Activity Coefficient  \gamma_w  vs Ambient RH and Temperature', ...
      'FontSize', 18, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = '\gamma_w'; cb.Label.FontSize = 14; cb.Label.FontWeight = 'bold';
colormap(turbo);
view(-45, 30);
grid on; box on;
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_gammaw_3D_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_gammaw_3D_RH_temp.fig'));
fprintf('Saved: LiCl_gammaw_3D_RH_temp.png\n');

%% ── Plot 6: gamma_w vs RH at fixed temperatures ──────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for i = 1:n_T
    [~, idx] = min(abs(T_vec - T_fixed(i)));
    gw    = gamma_w_grid(idx, :)';
    valid = ~isnan(gw);
    if ~any(valid), continue; end

    rh_valid = RH_vec(valid)*100;
    gw_valid = gw(valid);

    mk_idx = round(linspace(1, sum(valid), 8));
    mk_idx = min(mk_idx, length(rh_valid));

    plot(rh_valid, gw_valid, 'LineWidth', 2.5, 'Color', colors_T(i,:), ...
         'DisplayName', sprintf('T = %.0f°C', T_fixed(i)));

    plot(rh_valid(mk_idx), gw_valid(mk_idx), markers{i}, 'MarkerSize', 7, ...
         'Color', colors_T(i,:), 'MarkerFaceColor', colors_T(i,:), ...
         'HandleVisibility', 'off');
end

yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', '\gamma_w = 1 (ideal)');

colormap(gca, interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_T,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks      = linspace(0, 1, n_T);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f°C', x), T_fixed, 'UniformOutput', false);
cb.Label.String = 'Temperature'; cb.Label.FontSize = 13; cb.Label.FontWeight = 'bold';

xlabel('Ambient RH (%)',               'FontSize', 16, 'FontWeight', 'bold');
ylabel('\gamma_w  =  a_w / x_w',       'FontSize', 16, 'FontWeight', 'bold');
title('LiCl Water Activity Coefficient  \gamma_w  vs Ambient RH at Fixed Temperatures', ...
      'FontSize', 17, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_gammaw_vs_RH_fixed_T.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_gammaw_vs_RH_fixed_T.fig'));
fprintf('Saved: LiCl_gammaw_vs_RH_fixed_T.png\n');

%% ── Plot 7: gamma_w vs T at fixed ambient RH ─────────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for i = 1:n_RH
    [~, idx] = min(abs(RH_vec - RH_fixed(i)));
    gw    = gamma_w_grid(:, idx);
    valid = ~isnan(gw);
    if ~any(valid), continue; end

    t_valid  = T_vec(valid);
    gw_valid = gw(valid);

    mk_idx = round(linspace(1, sum(valid), 8));
    mk_idx = min(mk_idx, length(t_valid));

    plot(t_valid, gw_valid, 'LineWidth', 2.5, 'Color', colors_RH(i,:), ...
         'DisplayName', sprintf('RH = %.0f%%', RH_fixed(i)*100));

    plot(t_valid(mk_idx), gw_valid(mk_idx), markers{i}, 'MarkerSize', 7, ...
         'Color', colors_RH(i,:), 'MarkerFaceColor', colors_RH(i,:), ...
         'HandleVisibility', 'off');
end

yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', '\gamma_w = 1 (ideal)');

colormap(gca, interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_RH,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks      = linspace(0, 1, n_RH);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f%%', x*100), RH_fixed, 'UniformOutput', false);
cb.Label.String = 'Ambient RH'; cb.Label.FontSize = 13; cb.Label.FontWeight = 'bold';

xlabel('Temperature (°C)',             'FontSize', 16, 'FontWeight', 'bold');
ylabel('\gamma_w  =  a_w / x_w',       'FontSize', 16, 'FontWeight', 'bold');
title('LiCl Water Activity Coefficient  \gamma_w  vs Temperature at Fixed Ambient RH', ...
      'FontSize', 17, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_gammaw_vs_T_fixed_RH.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_gammaw_vs_T_fixed_RH.fig'));
fprintf('Saved: LiCl_gammaw_vs_T_fixed_RH.png\n');

%% ── Summary ──────────────────────────────────────────────────────────────────
fprintf('\n=== Summary Statistics ===\n');
fprintf('x_w range:     %.4f – %.4f\n', min(x_w_grid(:),[],'omitnan'), max(x_w_grid(:),[],'omitnan'));
fprintf('gamma_w range: %.4f – %.4f\n', min(gamma_w_grid(:),[],'omitnan'), max(gamma_w_grid(:),[],'omitnan'));

fprintf('\nTemperature effect on water uptake:\n');
fprintf('(At constant ambient RH, how does x_w change from 25°C to 100°C?)\n');
for rh_test = [0.35, 0.50, 0.65, 0.80]
    [~, ri]    = min(abs(RH_vec - rh_test));
    [~, ti_lo] = min(abs(T_vec - 25));
    [~, ti_hi] = min(abs(T_vec - 100));
    xw_lo = x_w_grid(ti_lo, ri);
    xw_hi = x_w_grid(ti_hi, ri);
    if ~isnan(xw_lo) && ~isnan(xw_hi)
        fprintf('  RH = %.0f%%:  x_w(25°C) = %.4f,  x_w(100°C) = %.4f,  Δx_w = %+.4f\n', ...
            rh_test*100, xw_lo, xw_hi, xw_hi - xw_lo);
    end
end

fprintf('\nAll figures saved to:\n%s\n', fig_out_dir);
fprintf('Done.\n');
