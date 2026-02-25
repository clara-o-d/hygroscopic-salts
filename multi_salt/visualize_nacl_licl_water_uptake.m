close all
clear
clc

% Script: visualize_nacl_licl_water_uptake.m
%
% Visualizes the molar water uptake  x_w  and water activity coefficient
% gamma_w for the NaCl + LiCl system at 25°C.
%
% Physical basis
% ──────────────
%   aw = gamma_w * x_w
%
%   where  x_w = n_water / (n_water + n_ions)
%              = 55.51   / (55.51   + 2*m_NaCl + 2*m_LiCl)
%
%   Both NaCl and LiCl are 1:1 electrolytes (nu = 2 each), so for a
%   mixture with total molality m_total and LiCl mole fraction y_LiCl:
%       n_ions = 2*m_total  (independent of the NaCl/LiCl split)
%
%   At a given ambient RH the solution equilibrates to the molality where
%   aw(mf_NaCl, mf_LiCl) = RH.  The inversion is performed with fzero
%   using the bivariate polynomial in calculate_activity_NaCl_LiCl.
%
%   Pure NaCl is inverted via calculate_mf_NaCl (valid 0.762–0.993).
%   Pure LiCl is inverted via calculate_mf_LiCl (T = 25°C).
%
% Plots
% ─────
%   1. x_w vs RH  – one curve per NaCl:LiCl ratio, pure salts included
%   2. gamma_w vs RH – same curves
%   3. ionic mole fraction (1 - x_w) vs RH – shows "salting" effect
%   4. equilibrium molality vs RH (log scale)
%   5. water uptake per dollar of salt – g H₂O / $ vs RH
%        LiCl: $0.25–0.50 / g   |   NaCl: $0.01–0.03 / g
%        midpoint line + shaded cost-uncertainty band for each ratio

%% ── Paths ──────────────────────────────────────────────────────────────────
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_mixtures'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

fig_out_dir = fullfile(filepath, '..', 'figures', 'multi_salt');
if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end

%% ── Constants ───────────────────────────────────────────────────────────────
MWw     = 18.015;   % g/mol water
MW_NaCl = 58.443;   % g/mol
MW_LiCl = 42.394;   % g/mol
n_w0    = 1000 / MWw;   % mol water per kg  ≈ 55.51

T25 = 25;  % °C

% Salt cost ranges ($ per gram of dry salt)
cost_NaCl_lo  = 0.01;   cost_NaCl_hi  = 0.03;
cost_LiCl_lo  = 0.25;   cost_LiCl_hi  = 0.50;
cost_NaCl_mid = (cost_NaCl_lo + cost_NaCl_hi) / 2;   % $0.02 / g
cost_LiCl_mid = (cost_LiCl_lo + cost_LiCl_hi) / 2;   % $0.375 / g

%% ── Mixture configurations ──────────────────────────────────────────────────
% y_LiCl = mole fraction of LiCl in the dissolved-salt mixture
%   0   → pure NaCl    (inverted via calculate_mf_NaCl)
%   1   → pure LiCl    (inverted via calculate_mf_LiCl)
%   else → mixture polynomial

y_LiCl_list  = [0, 0.25, 0.50, 0.75, 1.0];
labels = {'Pure NaCl', 'NaCl:LiCl = 3:1', 'NaCl:LiCl = 1:1 (equal)', ...
          'NaCl:LiCl = 1:3', 'Pure LiCl'};

% Colour: orange-red (NaCl-rich) → navy (LiCl-rich)
n_cases  = length(y_LiCl_list);
colors   = interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], (1:n_cases)');
markers  = {'o', 's', '^', 'd', 'v'};

%% ── Calibration limits of mixture polynomial ────────────────────────────────
mf_NaCl_min = 0.005745;  mf_NaCl_max = 0.236635;
mf_LiCl_min = 0.004704;  mf_LiCl_max = 0.154883;

% RH vector – mixture polynomial covers ~ 0.76 – 0.98
RH_vec_mix  = linspace(0.762, 0.982, 120);
% Pure LiCl covers down to very low aw
RH_vec_licl = linspace(0.15,  0.98,  120);

fzero_opts = optimset('Display', 'off', 'TolX', 1e-8);

%% ── Helper: mass fractions from total molality and y_LiCl ──────────────────
mf_of_m = @(m, y) deal( ...
    (1-y)*m * MW_NaCl ./ ((1-y)*m*MW_NaCl + y*m*MW_LiCl + 1000), ...
    y    *m * MW_LiCl ./ ((1-y)*m*MW_NaCl + y*m*MW_LiCl + 1000));

%% ── Solve for equilibrium molality at each (RH, y) ─────────────────────────
fprintf('Inverting activity polynomials for each NaCl:LiCl ratio...\n');

results = struct();

for k = 1:n_cases
    y    = y_LiCl_list(k);
    tag  = sprintf('y%02d', round(y*100));

    if y == 0
        % Pure NaCl – use existing calculate_mf_NaCl (valid RH 0.762–0.993)
        RH_vec = linspace(0.763, 0.993, 120);
        m_eq  = nan(size(RH_vec));
        for j = 1:length(RH_vec)
            try
                mf_s     = calculate_mf_NaCl(RH_vec(j));
                m_eq(j)  = 1000 * mf_s / (MW_NaCl * (1 - mf_s));
            catch; end
        end

    elseif y == 1
        % Pure LiCl – use existing calculate_mf_LiCl (T = 25°C)
        RH_vec = RH_vec_licl;
        m_eq  = nan(size(RH_vec));
        for j = 1:length(RH_vec)
            try
                mf_s     = calculate_mf_LiCl(RH_vec(j), T25);
                m_eq(j)  = 1000 * mf_s / (MW_LiCl * (1 - mf_s));
            catch; end
        end

    else
        % Mixture – invert bivariate polynomial parameterised by m_total
        RH_vec = RH_vec_mix;
        m_eq  = nan(size(RH_vec));

        for j = 1:length(RH_vec)
            RH_t = RH_vec(j);
            residual = @(m) activity_mixture(m, y, MW_NaCl, MW_LiCl) - RH_t;
            try
                r_lo = residual(0.05);
                r_hi = residual(15);
                if sign(r_lo) ~= sign(r_hi)
                    m_sol = fzero(residual, [0.05, 15], fzero_opts);
                    if m_sol > 0
                        m_eq(j) = m_sol;
                    end
                end
            catch; end
        end
    end

    % Compute x_w and gamma_w from equilibrium molality
    % n_ions = 2*m_total (both salts contribute 2 ions per formula unit)
    x_w    = n_w0 ./ (n_w0 + 2.*m_eq);
    gamma_w = RH_vec ./ x_w;

    results.(tag) = struct('RH', RH_vec, 'm_eq', m_eq, 'x_w', x_w, ...
                            'gamma_w', gamma_w, 'y', y, 'label', labels{k}, ...
                            'color', colors(k,:), 'marker', markers{k});

    n_ok = sum(~isnan(m_eq));
    fprintf('  y_LiCl = %.2f  (%s): %d / %d RH points converged\n', ...
        y, labels{k}, n_ok, length(RH_vec));
end
fprintf('\n');

tags = fieldnames(results);

%% ── Plot 1: Molar water uptake  x_w  vs  RH ────────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for k = 1:length(tags)
    d   = results.(tags{k});
    ok  = ~isnan(d.x_w);
    if ~any(ok), continue; end

    rh  = d.RH(ok)*100;
    xw  = d.x_w(ok);

    mk_idx = round(linspace(1, sum(ok), 7));
    plot(rh, xw, 'LineWidth', 2.5, 'Color', d.color, 'DisplayName', d.label);
    plot(rh(mk_idx), xw(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

xlabel('Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Molar Water Uptake  x_w', 'FontSize', 16, 'FontWeight', 'bold');
title('Molar Water Uptake of NaCl, LiCl, and Mixtures vs RH  (25°C)', ...
      'FontSize', 19, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12);
ylim([0 1]);

% Colorbar indicating NaCl → LiCl gradient
colormap(gca, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks = linspace(0,1,n_cases);
cb.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb.Label.String = 'LiCl mole fraction in salt';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';

set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_xw_vs_RH.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_xw_vs_RH.fig'));
fprintf('Saved: NaCl_LiCl_xw_vs_RH.png\n');

%% ── Plot 2: Activity coefficient of water  gamma_w  vs  RH ─────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for k = 1:length(tags)
    d  = results.(tags{k});
    ok = ~isnan(d.gamma_w) & d.gamma_w > 0 & d.gamma_w <= 1.5;
    if ~any(ok), continue; end

    rh = d.RH(ok)*100;
    gw = d.gamma_w(ok);

    mk_idx = round(linspace(1, sum(ok), 7));
    plot(rh, gw, 'LineWidth', 2.5, 'Color', d.color, 'DisplayName', d.label);
    plot(rh(mk_idx), gw(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', '\gamma_w = 1 (ideal)');

xlabel('Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('\gamma_w  =  a_w / x_w', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity Coefficient  \gamma_w  for NaCl + LiCl  (25°C)', ...
      'FontSize', 19, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);

colormap(gca, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks = linspace(0,1,n_cases);
cb.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb.Label.String = 'LiCl mole fraction in salt';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';

set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_gammaw_vs_RH.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_gammaw_vs_RH.fig'));
fprintf('Saved: NaCl_LiCl_gammaw_vs_RH.png\n');

%% ── Plot 3: Ionic mole fraction  (1 - x_w)  vs  RH ─────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for k = 1:length(tags)
    d  = results.(tags{k});
    ok = ~isnan(d.x_w);
    if ~any(ok), continue; end

    rh   = d.RH(ok)*100;
    ions = 1 - d.x_w(ok);

    mk_idx = round(linspace(1, sum(ok), 7));
    plot(rh, ions, 'LineWidth', 2.5, 'Color', d.color, 'DisplayName', d.label);
    plot(rh(mk_idx), ions(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

xlabel('Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Ionic Mole Fraction  (1 - x_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Ionic Mole Fraction in Equilibrated NaCl + LiCl Solutions  (25°C)', ...
      'FontSize', 19, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
ylim([0 0.6]);

colormap(gca, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks = linspace(0,1,n_cases);
cb.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb.Label.String = 'LiCl mole fraction in salt';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';

set(gca, 'FontSize', 14); set(gcf, 'color', 'w');

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_ionic_fraction_vs_RH.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_ionic_fraction_vs_RH.fig'));
fprintf('Saved: NaCl_LiCl_ionic_fraction_vs_RH.png\n');

%% ── Plot 4: Equilibrium molality vs RH ─────────────────────────────────────
figure('Position', [100 100 1200 900]);
hold on; grid on; box on;

for k = 1:length(tags)
    d  = results.(tags{k});
    ok = ~isnan(d.m_eq);
    if ~any(ok), continue; end

    rh = d.RH(ok)*100;
    m  = d.m_eq(ok);

    mk_idx = round(linspace(1, sum(ok), 7));
    plot(rh, m, 'LineWidth', 2.5, 'Color', d.color, 'DisplayName', d.label);
    plot(rh(mk_idx), m(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

xlabel('Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Equilibrium Molality  m_{total}  (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium Molality of NaCl + LiCl Solutions vs RH  (25°C)', ...
      'FontSize', 19, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);

colormap(gca, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb = colorbar; clim([0 1]);
cb.Ticks = linspace(0,1,n_cases);
cb.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb.Label.String = 'LiCl mole fraction in salt';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';

set(gca, 'FontSize', 14); set(gcf, 'color', 'w');
set(gca, 'YScale', 'log');

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_molality_vs_RH.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_molality_vs_RH.fig'));
fprintf('Saved: NaCl_LiCl_molality_vs_RH.png\n');

%% ── Plot 5: Water uptake per dollar of salt  ────────────────────────────────
%
% Metric: grams of water absorbed per dollar of dry salt invested
%
%   cost_per_kg_water($) = (1-y)*m_eq*MW_NaCl*cost_NaCl
%                        + y    *m_eq*MW_LiCl*cost_LiCl
%
%   water_per_dollar (g H₂O / $) = 1000 / cost_per_kg_water
%
% "Dry-salt" basis: we charge for the salt initially loaded; the water it
% absorbs from the air is then the output.  Lines show nominal (midpoint)
% costs; shaded bands span the lo–hi price ranges.

figure('Position', [100 100 1400 950]);
ax_main = axes('Position', [0.08 0.13 0.78 0.78]);
hold(ax_main, 'on'); grid(ax_main, 'on'); box(ax_main, 'on');

for k = 1:length(tags)
    d   = results.(tags{k});
    y   = d.y;
    ok  = ~isnan(d.m_eq);
    if ~any(ok), continue; end

    rh  = d.RH(ok) * 100;
    m   = d.m_eq(ok);

    % Cost per kg water  ($ / kg H₂O)
    cost_mid  = (1-y)*m*MW_NaCl*cost_NaCl_mid + y*m*MW_LiCl*cost_LiCl_mid;
    cost_best = (1-y)*m*MW_NaCl*cost_NaCl_lo  + y*m*MW_LiCl*cost_LiCl_lo;
    cost_worst= (1-y)*m*MW_NaCl*cost_NaCl_hi  + y*m*MW_LiCl*cost_LiCl_hi;

    % Guard against zero cost (pure water edge)
    cost_mid (cost_mid  <= 0) = NaN;
    cost_best(cost_best <= 0) = NaN;
    cost_worst(cost_worst <= 0) = NaN;

    wpd_mid   = 1000 ./ cost_mid;
    wpd_best  = 1000 ./ cost_best;
    wpd_worst = 1000 ./ cost_worst;

    % Shaded cost-range band
    rh_patch  = [rh, fliplr(rh)];
    wpd_patch = [wpd_best, fliplr(wpd_worst)];
    valid_patch = ~any(isnan(wpd_patch));
    if valid_patch
        fill(ax_main, rh_patch, wpd_patch, d.color, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Midpoint line with markers
    mk_idx = round(linspace(1, sum(ok), 7));
    plot(ax_main, rh, wpd_mid, 'LineWidth', 2.5, 'Color', d.color, ...
         'DisplayName', d.label);
    plot(ax_main, rh(mk_idx), wpd_mid(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

xlabel(ax_main, 'Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(ax_main, 'Water Uptake per Dollar  (g H_2O / $)', 'FontSize', 16, 'FontWeight', 'bold');
title(ax_main, {'Cost-Effectiveness of Water Uptake: NaCl vs LiCl Mixtures  (25°C)', ...
      sprintf('NaCl: $%.2f–$%.2f/g  |  LiCl: $%.2f–$%.2f/g   (shaded = cost range, line = midpoint)', ...
      cost_NaCl_lo, cost_NaCl_hi, cost_LiCl_lo, cost_LiCl_hi)}, ...
      'FontSize', 15, 'FontWeight', 'bold');
legend(ax_main, 'Location', 'northwest', 'FontSize', 12);
set(ax_main, 'YScale', 'log', 'FontSize', 14);
set(gcf, 'color', 'w');

% Colorbar
colormap(ax_main, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb5 = colorbar(ax_main); clim(ax_main, [0 1]);
cb5.Ticks = linspace(0,1,n_cases);
cb5.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb5.Label.String = 'LiCl mole fraction in salt';
cb5.Label.FontSize = 12; cb5.Label.FontWeight = 'bold';

% Annotation box: reference prices
annotation('textbox', [0.10 0.14 0.30 0.10], 'String', ...
    {sprintf('NaCl cost: $%.2f – $%.2f / g (mid $%.3f/g)', cost_NaCl_lo, cost_NaCl_hi, cost_NaCl_mid), ...
     sprintf('LiCl cost: $%.2f – $%.2f / g (mid $%.3f/g)', cost_LiCl_lo, cost_LiCl_hi, cost_LiCl_mid)}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.6 0.6 0.6], ...
    'FitBoxToText', 'on');

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_per_dollar.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_water_per_dollar.fig'));
fprintf('Saved: NaCl_LiCl_water_per_dollar.png\n');

% ── Bonus panel: same data but overlaid on a linear scale inset ─────────────
% Create a second figure showing just the high-RH region where NaCl competes
figure('Position', [100 100 1400 950]);
ax2 = axes('Position', [0.08 0.13 0.78 0.78]);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

for k = 1:length(tags)
    d   = results.(tags{k});
    y   = d.y;
    ok  = ~isnan(d.m_eq) & d.RH >= 0.75;
    if ~any(ok), continue; end

    rh  = d.RH(ok) * 100;
    m   = d.m_eq(ok);

    cost_mid  = (1-y)*m*MW_NaCl*cost_NaCl_mid + y*m*MW_LiCl*cost_LiCl_mid;
    cost_best = (1-y)*m*MW_NaCl*cost_NaCl_lo  + y*m*MW_LiCl*cost_LiCl_lo;
    cost_worst= (1-y)*m*MW_NaCl*cost_NaCl_hi  + y*m*MW_LiCl*cost_LiCl_hi;

    cost_mid (cost_mid  <= 0) = NaN;
    cost_best(cost_best <= 0) = NaN;
    cost_worst(cost_worst <= 0) = NaN;

    wpd_mid   = 1000 ./ cost_mid;
    wpd_best  = 1000 ./ cost_best;
    wpd_worst = 1000 ./ cost_worst;

    rh_patch  = [rh, fliplr(rh)];
    wpd_patch = [wpd_best, fliplr(wpd_worst)];
    if ~any(isnan(wpd_patch))
        fill(ax2, rh_patch, wpd_patch, d.color, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    mk_idx = round(linspace(1, sum(ok), 7));
    plot(ax2, rh, wpd_mid, 'LineWidth', 2.5, 'Color', d.color, ...
         'DisplayName', d.label);
    plot(ax2, rh(mk_idx), wpd_mid(mk_idx), d.marker, 'MarkerSize', 7, ...
         'Color', d.color, 'MarkerFaceColor', d.color, 'HandleVisibility', 'off');
end

xlabel(ax2, 'Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(ax2, 'Water Uptake per Dollar  (g H_2O / $)', 'FontSize', 16, 'FontWeight', 'bold');
title(ax2, {'Cost-Effectiveness of Water Uptake: High-RH Region (≥75%)  (25°C)', ...
      sprintf('NaCl: $%.2f–$%.2f/g  |  LiCl: $%.2f–$%.2f/g', ...
      cost_NaCl_lo, cost_NaCl_hi, cost_LiCl_lo, cost_LiCl_hi)}, ...
      'FontSize', 15, 'FontWeight', 'bold');
legend(ax2, 'Location', 'northwest', 'FontSize', 12);
set(ax2, 'FontSize', 14);   % linear scale to show crossover clearly
set(gcf, 'color', 'w');
xlim(ax2, [75 100]);

colormap(ax2, interp1([1; n_cases], [0.85 0.20 0.05; 0.05 0.20 0.70], linspace(1,n_cases,256)'));
cb6 = colorbar(ax2); clim(ax2, [0 1]);
cb6.Ticks = linspace(0,1,n_cases);
cb6.TickLabels = arrayfun(@(y) sprintf('y_{LiCl}=%.2f',y), y_LiCl_list, 'UniformOutput',false);
cb6.Label.String = 'LiCl mole fraction in salt';
cb6.Label.FontSize = 12; cb6.Label.FontWeight = 'bold';

saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_per_dollar_highRH.png'));
savefig(     fullfile(fig_out_dir, 'NaCl_LiCl_water_per_dollar_highRH.fig'));
fprintf('Saved: NaCl_LiCl_water_per_dollar_highRH.png\n');

fprintf('\nAll figures saved to:\n%s\n', fig_out_dir);
fprintf('Done.\n');

%% ── Local helper functions ──────────────────────────────────────────────────

function aw = activity_mixture(m_total, y_LiCl, MW_NaCl, MW_LiCl)
% Evaluate calculate_activity_NaCl_LiCl at a given total molality and
% LiCl mole fraction  y_LiCl  (mole fraction of LiCl among the two salts).
    denom   = (1-y_LiCl)*m_total*MW_NaCl + y_LiCl*m_total*MW_LiCl + 1000;
    mf_NaCl = (1-y_LiCl)*m_total*MW_NaCl / denom;
    mf_LiCl =    y_LiCl *m_total*MW_LiCl / denom;
    aw = calculate_activity_NaCl_LiCl(mf_NaCl, mf_LiCl);
end
