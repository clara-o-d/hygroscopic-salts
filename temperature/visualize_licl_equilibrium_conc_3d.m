close all
clear
clc

% Script: visualize_licl_equilibrium_conc_3d.m
%
% Unlike visualize_licl_temperature_activity_3d.m (which fixes LiCl molality),
% this script accounts for the fact that a real LiCl solution equilibrates
% with its surroundings: at each (ambient RH, temperature) pair the solution
% absorbs or releases water until its water activity matches the ambient RH.
%
% The equilibrium molality at (RH_target, T) is found by inverting
%     calculate_activity_temperature_LiCl(mf, T) = RH_target
% using fzero, then converting the equilibrium mf to molality.
%
% Plots produced:
%   1.  3D surface  – equilibrium molality vs (RH, T)
%   2.  Contour map – equilibrium molality vs (RH, T)
%   3.  2D lines    – equilibrium molality vs T, one curve per ambient RH
%         (cool→warm colour scale; unique markers)
%   4.  2D lines    – equilibrium molality vs ambient RH, one curve per T
%         (light→dark colour scale; unique markers)
%   5.  3D surface  – equilibrium mass fraction vs (RH, T)
%   6.  2D lines    – aw of equilibrated solution vs T (verification: should
%         equal the target RH line, confirming the inversion is correct)

%% Setup paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_temperature'));
addpath(fullfile(filepath, '..', 'util'));

fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% Constants
MW_LiCl = 42.394;  % g/mol

%% Calibration limits of the fitted polynomial
mf_min_cal = 0.041481;
mf_max_cal = 0.440682;
T_min_cal  = 25.0;
T_max_cal  = 100.0;

%% Define ranges
RH_vec = linspace(0.15, 0.99, 60);   % ambient water activity
T_vec  = linspace(25, 100, 60);       % temperature (°C)

[RH_grid, T_grid] = meshgrid(RH_vec, T_vec);

fprintf('Solving for equilibrium LiCl concentration...\n');
fprintf('RH range:          %.2f to %.2f\n', min(RH_vec), max(RH_vec));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(T_vec), max(T_vec));
fprintf('Grid size:         %d x %d = %d points\n', ...
    length(RH_vec), length(T_vec), numel(RH_grid));

%% Solve for equilibrium mass fraction at every (RH, T) point
mf_eq_grid = nan(size(RH_grid));
m_eq_grid  = nan(size(RH_grid));   % mol / kg H2O

fzero_opts = optimset('Display', 'off', 'TolX', 1e-7);

for i = 1:numel(RH_grid)
    RH_target = RH_grid(i);
    T         = T_grid(i);

    % Residual function: aw(mf,T) - RH_target = 0
    residual = @(mf) calculate_activity_temperature_LiCl(mf, T) - RH_target;

    % Skip if T is outside the calibration range of the polynomial
    if T < T_min_cal || T > T_max_cal
        continue;
    end

    % Bracket search: very high RH → dilute solution (small mf);
    %                 very low  RH → concentrated solution (large mf).
    % Use [mf_min_cal, mf_max_cal] as the search bracket.
    try
        % Check that a sign change exists within the calibrated mf range
        r_lo = residual(mf_min_cal);
        r_hi = residual(mf_max_cal);

        if sign(r_lo) ~= sign(r_hi)
            mf_sol = fzero(residual, [mf_min_cal, mf_max_cal], fzero_opts);
        else
            % No bracket found – target RH is outside achievable range at
            % this temperature; leave as NaN
            continue;
        end

        if mf_sol >= mf_min_cal && mf_sol <= mf_max_cal
            mf_eq_grid(i) = mf_sol;
            % Convert to molality: mf = m*MW / (1000 + m*MW)  →  m = 1000*mf / (MW*(1-mf))
            m_eq_grid(i) = 1000 * mf_sol / (MW_LiCl * (1 - mf_sol));
        end
    catch
        % fzero failed; leave as NaN
    end
end

n_valid = sum(~isnan(m_eq_grid(:)));
fprintf('Converged:  %d / %d points\n\n', n_valid, numel(RH_grid));

%% ── Plot 1: 3D surface – equilibrium molality ─────────────────────────────
figure('Position', [100 100 1200 900]);
surf(RH_grid*100, T_grid, m_eq_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);

xlabel('Ambient RH (%)',             'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)',           'FontSize', 16, 'FontWeight', 'bold');
zlabel('Equilibrium LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium LiCl Molality vs Ambient RH and Temperature', ...
      'FontSize', 20, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = 'Molality (mol/kg)';
cb.Label.FontSize = 14;  cb.Label.FontWeight = 'bold';
colormap(jet);
view(-45, 30);
grid on; box on;
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_molality_3D_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_molality_3D_RH_temp.fig'));
fprintf('Saved: LiCl_eq_molality_3D_RH_temp.png\n');

%% ── Plot 2: Contour map – equilibrium molality ────────────────────────────
figure('Position', [100 100 1200 900]);
[C, h] = contourf(RH_grid*100, T_grid, m_eq_grid, 20);
clabel(C, h, 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');

xlabel('Ambient RH (%)',   'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium LiCl Molality Contours', 'FontSize', 20, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = 'Equilibrium Molality (mol/kg)';
cb.Label.FontSize = 14;  cb.Label.FontWeight = 'bold';
colormap(jet);
grid on; box on;
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_molality_contour_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_molality_contour_RH_temp.fig'));
fprintf('Saved: LiCl_eq_molality_contour_RH_temp.png\n');

%% ── Plot 3: 2D cross-sections – equilibrium molality vs T at fixed RH ────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

RH_fixed = [0.20, 0.35, 0.50, 0.65, 0.80, 0.95];
n_RH     = length(RH_fixed);

% Cool (dry/concentrated) → warm (humid/dilute) colour ramp
colors_RH = interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], (1:n_RH)');
markers   = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_RH
    [~, idx] = min(abs(RH_vec - RH_fixed(i)));
    y = m_eq_grid(:, idx);           % molality at this RH, varying T

    mk_idx = round(linspace(1, length(T_vec), 8));
    valid  = ~isnan(y);

    if any(valid)
        plot(T_vec(valid), y(valid), 'LineWidth', 2.5, ...
             'Color', colors_RH(i,:), ...
             'DisplayName', sprintf('RH = %.0f%%', RH_vec(idx)*100));

        mk_valid = mk_idx(~isnan(y(mk_idx)));
        if ~isempty(mk_valid)
            plot(T_vec(mk_valid), y(mk_valid), markers{i}, ...
                 'MarkerSize', 7, 'Color', colors_RH(i,:), ...
                 'MarkerFaceColor', colors_RH(i,:), 'HandleVisibility', 'off');
        end
    end
end

% Colorbar indicating RH level
colormap(gca, interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_RH,256)'));
cb = colorbar;
cb.Ticks      = linspace(0,1,n_RH);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f%%', x*100), RH_fixed, 'UniformOutput', false);
cb.Label.String    = 'Ambient RH';
cb.Label.FontSize  = 13;  cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('Temperature (°C)',                       'FontSize', 16, 'FontWeight', 'bold');
ylabel('Equilibrium LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium LiCl Molality vs Temperature at Fixed Ambient RH', ...
      'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12);
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_molality_vs_T_fixed_RH.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_molality_vs_T_fixed_RH.fig'));
fprintf('Saved: LiCl_eq_molality_vs_T_fixed_RH.png\n');

%% ── Plot 4: 2D cross-sections – equilibrium molality vs RH at fixed T ────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

T_fixed = [25, 40, 55, 70, 85, 100];
n_T     = length(T_fixed);

% Cool-to-warm ramp: blue (cold) → red (hot)
colors_T = interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], (1:n_T)');

for i = 1:n_T
    [~, idx] = min(abs(T_vec - T_fixed(i)));
    y = m_eq_grid(idx, :)';          % molality at this T, varying RH

    mk_idx   = round(linspace(1, length(RH_vec), 8));
    valid    = ~isnan(y);

    if any(valid)
        plot(RH_vec(valid)*100, y(valid), 'LineWidth', 2.5, ...
             'Color', colors_T(i,:), ...
             'DisplayName', sprintf('T = %.0f°C', T_vec(idx)));

        mk_valid = mk_idx(~isnan(y(mk_idx)));
        if ~isempty(mk_valid)
            plot(RH_vec(mk_valid)*100, y(mk_valid), markers{i}, ...
                 'MarkerSize', 7, 'Color', colors_T(i,:), ...
                 'MarkerFaceColor', colors_T(i,:), 'HandleVisibility', 'off');
        end
    end
end

colormap(gca, interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_T,256)'));
cb = colorbar;
cb.Ticks      = linspace(0,1,n_T);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f°C', x), T_fixed, 'UniformOutput', false);
cb.Label.String   = 'Temperature';
cb.Label.FontSize = 13;  cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('Ambient RH (%)',                          'FontSize', 16, 'FontWeight', 'bold');
ylabel('Equilibrium LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium LiCl Molality vs Ambient RH at Fixed Temperatures', ...
      'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_molality_vs_RH_fixed_T.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_molality_vs_RH_fixed_T.fig'));
fprintf('Saved: LiCl_eq_molality_vs_RH_fixed_T.png\n');

%% ── Plot 5: 3D surface – equilibrium mass fraction ────────────────────────
figure('Position', [100 100 1200 900]);
surf(RH_grid*100, T_grid, mf_eq_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);

xlabel('Ambient RH (%)',         'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)',       'FontSize', 16, 'FontWeight', 'bold');
zlabel('Equilibrium Mass Fraction of LiCl', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium LiCl Mass Fraction vs Ambient RH and Temperature', ...
      'FontSize', 20, 'FontWeight', 'bold');

cb = colorbar;
cb.Label.String = 'Mass Fraction';
cb.Label.FontSize = 14;  cb.Label.FontWeight = 'bold';
colormap(parula);
view(-45, 30);
grid on; box on;
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_massfrac_3D_RH_temp.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_massfrac_3D_RH_temp.fig'));
fprintf('Saved: LiCl_eq_massfrac_3D_RH_temp.png\n');

%% ── Plot 6: Verification – aw of equilibrated solution vs T ──────────────
% At equilibrium, calculate_activity_temperature_LiCl(mf_eq, T) should equal
% the target RH. Plot a few lines to confirm visually.
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

% Reuse the same RH_fixed colours from Plot 3
for i = 1:n_RH
    [~, idx] = min(abs(RH_vec - RH_fixed(i)));

    aw_verify = nan(length(T_vec), 1);
    for ti = 1:length(T_vec)
        mf_v = mf_eq_grid(ti, idx);
        T_v  = T_vec(ti);
        if ~isnan(mf_v)
            aw_verify(ti) = calculate_activity_temperature_LiCl(mf_v, T_v);
        end
    end

    valid = ~isnan(aw_verify);
    if any(valid)
        plot(T_vec(valid), aw_verify(valid)*100, 'LineWidth', 2.5, ...
             'Color', colors_RH(i,:), ...
             'DisplayName', sprintf('Target RH = %.0f%%', RH_fixed(i)*100));
        % Dashed horizontal line showing the target
        yline(RH_fixed(i)*100, '--', 'Color', colors_RH(i,:), ...
              'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

colormap(gca, interp1([1; n_RH], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_RH,256)'));
cb = colorbar;
cb.Ticks      = linspace(0,1,n_RH);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f%%', x*100), RH_fixed, 'UniformOutput', false);
cb.Label.String   = 'Ambient RH';
cb.Label.FontSize = 13;  cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Computed Water Activity of Equilibrated Solution (%)', ...
       'FontSize', 16, 'FontWeight', 'bold');
title({'Verification: a_w of Equilibrated LiCl Solution vs Temperature', ...
       '(solid = computed, dashed = target RH; should overlap)'}, ...
      'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 12);
ylim([0 100]);
set(gca, 'FontSize', 14);  set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'LiCl_eq_aw_verify_vs_T.png'));
savefig(     fullfile(fig_out_dir, 'LiCl_eq_aw_verify_vs_T.fig'));
fprintf('Saved: LiCl_eq_aw_verify_vs_T.png\n');

%% ── Summary ───────────────────────────────────────────────────────────────
fprintf('\n=== Summary Statistics ===\n');
fprintf('Equilibrium molality range: %.3f to %.3f mol/kg\n', ...
    min(m_eq_grid(:), [], 'omitnan'), max(m_eq_grid(:), [], 'omitnan'));

fprintf('\nTemperature effect on equilibrium concentration:\n');
fprintf('(At constant ambient RH, how much does equilibrium molality change?)\n');
for rh_test = [0.35, 0.50, 0.65, 0.80]
    [~, ri] = min(abs(RH_vec - rh_test));
    [~, ti_lo] = min(abs(T_vec - 25));
    [~, ti_hi] = min(abs(T_vec - 100));
    m_lo = m_eq_grid(ti_lo, ri);
    m_hi = m_eq_grid(ti_hi, ri);
    if ~isnan(m_lo) && ~isnan(m_hi)
        fprintf('  RH = %.0f%%:  m(25°C) = %.3f, m(100°C) = %.3f, Δm = %+.3f mol/kg\n', ...
            rh_test*100, m_lo, m_hi, m_hi - m_lo);
    end
end

fprintf('\nAll figures saved to:\n%s\n', fig_out_dir);
fprintf('\nVisualization complete!\n');
