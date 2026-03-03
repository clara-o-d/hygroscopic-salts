close all
clear
clc

% ─────────────────────────────────────────────────────────────────────────────
%  licl_bacl2_case_study.m
%
%  Case study for the LiCl + BaCl₂ + H₂O system at 25 °C.
%
%  DATA SOURCE
%    data/osmotic_activity_coefficients.xlsx  —  Table 3 from the
%    literature: osmotic coefficient φ and mean activity coefficients
%    γ_LiCl, γ_BaCl₂ as functions of ionic strength I and the BaCl₂
%    mole fraction y_B among the dissolved salts.
%
%  WHAT THIS SCRIPT DOES
%    1. Load & parse Table 3 (I vs y_B grid).
%    2. Fit a 2-D polynomial (degree 2 in I, degree 2 in y_B) to each of
%       φ, log₁₀γ_LiCl, log₁₀γ_BaCl₂.
%    3. Plot: data vs polynomial fit — so you can verify accuracy.
%    4. Plot: water activity aᵥ vs total molality for
%         • pure LiCl       (y_B = 0)
%         • pure BaCl₂      (y_B = 1)
%         • 50:50 molar mix (y_B = 0.5, equal molalities)
%    5. Plot: molar water uptake x_w vs RH  (ionic mole-fraction basis)
%    6. Plot: mass-based water uptake mf_w vs RH
%
%  NOTATION
%    y_B  = m_BaCl₂ / (m_LiCl + m_BaCl₂)   salt mole fraction of BaCl₂
%    I    = m_LiCl + 3·m_BaCl₂              ionic strength
%    m_total = m_LiCl + m_BaCl₂             total dissolved-salt molality
%    m_total = I / (1 + 2·y_B)
%    ν_ions  = 2·m_LiCl + 3·m_BaCl₂ = (2 + y_B)·m_total
%    aᵥ      = exp(−φ · ν_ions · Mw / 1000)
% ─────────────────────────────────────────────────────────────────────────────

%% ── Paths ────────────────────────────────────────────────────────────────────
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'calculate_activity_mixtures'));
addpath(fullfile(filepath, '..', 'util'));

fig_out_dir = fullfile(filepath, '..', 'figures', 'multi_salt');
if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end

%% ── Constants ────────────────────────────────────────────────────────────────
MWw      = 18.015;     % g/mol  water
MW_LiCl  = 42.394;    % g/mol  LiCl
MW_BaCl2 = 208.23;    % g/mol  BaCl₂
n_w0     = 1000/MWw;  % mol water per kg  (≈ 55.51)
T25      = 25;         % °C

% colours used in every plot
col_LiCl  = [0.10 0.35 0.85];   % blue
col_BaCl2 = [0.85 0.15 0.10];   % red
col_mix   = [0.10 0.65 0.15];   % green

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 1 – Load & parse Table 3
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Loading Table 3 from osmotic_activity_coefficients.xlsx ...\n');

data_path = fullfile(filepath, '..', 'data', 'osmotic_activity_coefficients.xlsx');
raw = readcell(data_path, ...
    'Sheet',  'Osmotic Activity Coeffs', ...
    'Range',  'A4:E33');          % 30 data rows (skip header rows 1-3)

% Column layout:  A = I (mol/kg),  B = y_B,  C = φ,
%                 D = log₁₀γ_LiCl,  E = log₁₀γ_BaCl₂
% I is printed only in the first of each 6-row block; fill forward.

n_rows   = size(raw, 1);    % 30
I_data   = nan(n_rows,1);
yB_data  = nan(n_rows,1);
phi_data = nan(n_rows,1);
lgA_data = nan(n_rows,1);   % log₁₀γ_LiCl  (salt A)
lgB_data = nan(n_rows,1);   % log₁₀γ_BaCl₂ (salt B)

I_cur = NaN;
for k = 1:n_rows
    v = raw{k,1};
    if isnumeric(v) && ~isempty(v) && isfinite(v)
        I_cur = v;            % new block starts
    end
    I_data(k)   = I_cur;
    yB_data(k)  = raw{k,2};
    phi_data(k) = raw{k,3};
    lgA_data(k) = raw{k,4};
    lgB_data(k) = raw{k,5};
end

I_uniq  = unique(I_data);    % [0.5 1 2 3 4]
yB_uniq = unique(yB_data);   % [0 0.2 0.4 0.6 0.8 1]

fprintf('  Loaded %d points.  I = %s mol/kg\n', ...
    n_rows, mat2str(I_uniq'));
fprintf('                      y_B = %s\n\n', mat2str(yB_uniq'));

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 2 – 2-D polynomial fit
% ═══════════════════════════════════════════════════════════════════════════
% Model (total degree ≤ 2 in each variable, plus cross term):
%   f(I,y) = p1 + p2·I + p3·I² + p4·y + p5·y² + p6·I·y
% 6 coefficients, 30 data points → well-determined.

build_X    = @(Iv, yv) [ones(numel(Iv),1), Iv(:), Iv(:).^2, ...
                         yv(:), yv(:).^2,  Iv(:).*yv(:)];
eval_2D    = @(p, Iv, yv) p(1) + p(2)*Iv + p(3)*Iv.^2 + ...
                            p(4)*yv + p(5)*yv.^2 + p(6)*Iv.*yv;

X = build_X(I_data, yB_data);       % 30 × 6 design matrix

p_phi  = X \ phi_data;
p_lgA  = X \ lgA_data;
p_lgB  = X \ lgB_data;

rmse_phi = sqrt(mean((phi_data  - eval_2D(p_phi, I_data,yB_data)).^2));
rmse_lgA = sqrt(mean((lgA_data  - eval_2D(p_lgA, I_data,yB_data)).^2));
rmse_lgB = sqrt(mean((lgB_data  - eval_2D(p_lgB, I_data,yB_data)).^2));

fprintf('2-D polynomial fit RMSE:\n');
fprintf('   phi (osmotic coeff):   %.5f\n',  rmse_phi);
fprintf('   log10(gamma_LiCl):     %.5f\n',  rmse_lgA);
fprintf('   log10(gamma_BaCl2):    %.5f\n\n',rmse_lgB);

%% ── Helper: water activity from (I, y_B) via fitted φ ─────────────────────
% aᵥ = exp( −φ · ν_ions · Mw/1000 )
% ν_ions = (2+y_B)·m_total,  m_total = I/(1+2·y_B)
aw_IyB = @(Iv, yv) exp( ...
    -eval_2D(p_phi, Iv, yv) .* (2 + yv) .* (Iv ./ (1 + 2*yv)) .* MWw/1000 );

%% ═══════════════════════════════════════════════════════════════════════════
%  PLOT 1 – Polynomial fit validation
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Generating Plot 1: fit validation ...\n');

colors_yB = lines(numel(yB_uniq));
I_fine    = linspace(0.35, 4.3, 300);

fig1 = figure('Position', [60 60 1450 500], 'Color', 'w');

quants  = {phi_data,  lgA_data,           lgB_data};
p_coefs = {p_phi,     p_lgA,              p_lgB};
ylab    = {'\phi  (osmotic coefficient)', ...
           'log_{10} \gamma_{LiCl}', ...
           'log_{10} \gamma_{BaCl_2}'};
ttl     = {'Osmotic Coefficient  \phi', ...
           'LiCl Mean Activity Coefficient', ...
           'BaCl_2 Mean Activity Coefficient'};
rmses   = [rmse_phi, rmse_lgA, rmse_lgB];

for qi = 1:3
    subplot(1,3,qi);
    hold on; grid on; box on;
    q_data = quants{qi};
    p_q    = p_coefs{qi};

    for k = 1:numel(yB_uniq)
        yB   = yB_uniq(k);
        mask = abs(yB_data - yB) < 0.01;

        % --- data points ---
        plot(I_data(mask), q_data(mask), 'o', ...
            'Color', colors_yB(k,:), 'MarkerFaceColor', colors_yB(k,:), ...
            'MarkerSize', 8,  'DisplayName', sprintf('y_B=%.1f', yB));

        % --- fitted curve ---
        plot(I_fine, eval_2D(p_q, I_fine, yB*ones(size(I_fine))), '-', ...
            'Color', colors_yB(k,:), 'LineWidth', 1.6, ...
            'HandleVisibility', 'off');
    end

    xlabel('I  (mol kg^{-1} H_2O)', 'FontSize', 11);
    ylabel(ylab{qi},  'FontSize', 11);
    title(sprintf('%s\nRMSE = %.5f', ttl{qi}, rmses(qi)), 'FontSize', 11);
    if qi == 1
        legend('Location', 'northwest', 'FontSize', 8, 'NumColumns', 2);
    end
    set(gca, 'FontSize', 10);
end

sgtitle('LiCl + BaCl_2 + H_2O  —  Polynomial Fits vs. Table 3 Data  (25 °C)', ...
    'FontSize', 13, 'FontWeight', 'bold');

out = fullfile(fig_out_dir, 'LiCl_BaCl2_polynomial_fits');
saveas(fig1, [out '.png']);
savefig(fig1, [out '.fig']);
fprintf('   Saved: LiCl_BaCl2_polynomial_fits.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  PLOT 2 – Water activity vs total molality
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Generating Plot 2: aᵥ vs molality ...\n');

% Fine molality vectors
m_L  = linspace(0.05, 4.0,  400);   % pure LiCl  :  I = m,   y_B = 0
m_B  = linspace(0.05, 1.32, 400);   % pure BaCl₂ :  I = 3m,  y_B = 1
m_M  = linspace(0.05, 2.0,  400);   % 50:50 mix  :  I = 2m,  y_B = 0.5

aw_L = aw_IyB(m_L,       zeros(size(m_L)));
aw_B = aw_IyB(3*m_B,     ones(size(m_B)));
aw_M = aw_IyB(2*m_M,  0.5*ones(size(m_M)));

% Table 3 raw data converted to aᵥ
aw_raw = aw_IyB(I_data, yB_data);
m_raw  = I_data ./ (1 + 2*yB_data);   % m_total per data point

% Extract pure-salt and interpolated 50:50 markers
mk0  = abs(yB_data)       < 0.01;   % y_B = 0  (pure LiCl)
mk1  = abs(yB_data - 1.0) < 0.01;   % y_B = 1  (pure BaCl₂)
mk4  = abs(yB_data - 0.4) < 0.01;
mk6  = abs(yB_data - 0.6) < 0.01;
m_interp  = (m_raw(mk4)  + m_raw(mk6))  / 2;
aw_interp = (aw_raw(mk4) + aw_raw(mk6)) / 2;

fig2 = figure('Position', [100 100 1200 720], 'Color', 'w');
hold on; grid on; box on;

% Fitted lines
plot(m_L, aw_L, '-',  'Color', col_LiCl,  'LineWidth', 2.5, ...
    'DisplayName', 'Pure LiCl (poly fit)');
plot(m_B, aw_B, '-',  'Color', col_BaCl2, 'LineWidth', 2.5, ...
    'DisplayName', 'Pure BaCl_2 (poly fit)');
plot(m_M, aw_M, '-',  'Color', col_mix,   'LineWidth', 2.5, ...
    'DisplayName', 'LiCl + BaCl_2, 50:50 molar (poly fit)');

% Table 3 data markers (filled circles = pure salts, triangles = interp 50:50)
plot(m_raw(mk0), aw_raw(mk0), 'o', ...
    'Color', col_LiCl,  'MarkerFaceColor', col_LiCl,  'MarkerSize', 9, ...
    'HandleVisibility', 'off');
plot(m_raw(mk1), aw_raw(mk1), 'o', ...
    'Color', col_BaCl2, 'MarkerFaceColor', col_BaCl2, 'MarkerSize', 9, ...
    'HandleVisibility', 'off');
plot(m_interp, aw_interp, '^', ...
    'Color', col_mix,   'MarkerFaceColor', col_mix,   'MarkerSize', 9, ...
    'DisplayName', '50:50 interp. (y_B=0.4 & 0.6)');

xlabel('Total Molality  m_{total}  (mol kg^{-1} H_2O)', ...
    'FontSize', 14, 'FontWeight', 'bold');
ylabel('Water Activity  a_w', 'FontSize', 14, 'FontWeight', 'bold');
title({'Water Activity vs. Molality: LiCl, BaCl_2, and 50:50 Molar Mixture  (25 °C)', ...
       'm_{total} = m_{LiCl}+m_{BaCl_2}; 50:50 means m_{LiCl}=m_{BaCl_2}'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
ylim([0.70 1.01]);  xlim([0 4.3]);
set(gca, 'FontSize', 13);

annotation('textbox', [0.135 0.14 0.31 0.13], 'FitBoxToText', 'off', ...
    'String', {'Markers = Table 3 data'; ...
               'Lines = 2D polynomial fits'; ...
               '50:50 markers = interp. y_B=0.4 & 0.6'}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.6 0.6 0.6]);

out = fullfile(fig_out_dir, 'LiCl_BaCl2_aw_vs_molality');
saveas(fig2, [out '.png']);
savefig(fig2, [out '.fig']);
fprintf('   Saved: LiCl_BaCl2_aw_vs_molality.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  SECTION 3 – Invert activity functions to get equilibrium molality vs RH
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Inverting activity functions for equilibrium molality ...\n');

fz_opts = optimset('Display', 'off', 'TolX', 1e-9, 'TolFun', 1e-9);

% ── Pure LiCl: calculate_mf_LiCl covers RH 0.12 – 0.97 ─────────────────
RH_L = linspace(0.15, 0.97, 250);
m_L_eq = nan(size(RH_L));
for k = 1:numel(RH_L)
    try
        mf_s      = calculate_mf_LiCl(RH_L(k), T25);
        m_L_eq(k) = 1000 * mf_s / (MW_LiCl * (1 - mf_s));
    catch, end
end

% ── Pure BaCl₂: calculate_mf_BaCl2 covers RH 0.908 – 0.960 ─────────────
RH_B = linspace(0.909, 0.959, 80);
m_B_eq = nan(size(RH_B));
for k = 1:numel(RH_B)
    try
        mf_s      = calculate_mf_BaCl2(RH_B(k));
        m_B_eq(k) = 1000 * mf_s / (MW_BaCl2 * (1 - mf_s));
    catch, end
end

% ── 50:50 molar mix: invert calculate_activity_LiCl_BaCl2 ───────────────
% Valid mass-fraction range for the bivariate polynomial:
%   mf_LiCl  ≤ 0.162,  mf_BaCl₂ ≤ 0.2104
% For equal molalities (m_each), the limiting constraint is mf_BaCl₂ ≤ 0.2104
%   → m_each ≤ 1.35 mol/kg
% The minimum is set by mf_LiCl ≥ 0.0051  → m_each ≥ 0.124
%
% Corresponding aᵥ range is roughly 0.73 – 0.989 for this 50:50 path.

RH_M   = linspace(0.735, 0.989, 200);
m_M_eq = nan(size(RH_M));   % m_each = m_LiCl = m_BaCl₂

for k = 1:numel(RH_M)
    RH_t = RH_M(k);
    obj  = @(me) mix_activity_50(me, MW_LiCl, MW_BaCl2) - RH_t;
    try
        r0 = obj(0.124);   r1 = obj(1.35);
        if isfinite(r0) && isfinite(r1) && sign(r0) ~= sign(r1)
            m_M_eq(k) = fzero(obj, [0.124, 1.35], fz_opts);
        end
    catch, end
end

fprintf('   LiCl   converged at %d / %d RH points\n', sum(~isnan(m_L_eq)), numel(RH_L));
fprintf('   BaCl₂  converged at %d / %d RH points\n', sum(~isnan(m_B_eq)), numel(RH_B));
fprintf('   Mix    converged at %d / %d RH points\n\n', sum(~isnan(m_M_eq)), numel(RH_M));

%% ── Derived quantities ───────────────────────────────────────────────────
% Molar water uptake  x_w = n_water / (n_water + n_ions)
%   Pure LiCl:  n_ions = 2·m
%   Pure BaCl₂: n_ions = 3·m
%   50:50 mix:  n_ions = 2·m_each + 3·m_each = 5·m_each
xw_L = n_w0 ./ (n_w0 + 2*m_L_eq);
xw_B = n_w0 ./ (n_w0 + 3*m_B_eq);
xw_M = n_w0 ./ (n_w0 + 5*m_M_eq);

% Mass fraction of water  mf_w = 1000 / (1000 + mass_salt_per_kg_water)
%   50:50 mix: mass_salt = m_each·MW_LiCl + m_each·MW_BaCl₂ per kg water
mfw_L = 1000 ./ (1000 + m_L_eq * MW_LiCl);
mfw_B = 1000 ./ (1000 + m_B_eq * MW_BaCl2);
mfw_M = 1000 ./ (1000 + m_M_eq * (MW_LiCl + MW_BaCl2));

%% ═══════════════════════════════════════════════════════════════════════════
%  PLOT 3 – Molar water uptake  x_w  vs  RH
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Generating Plot 3: molar uptake vs RH ...\n');

fig3 = figure('Position', [120 120 1200 720], 'Color', 'w');
hold on; grid on; box on;

ok = ~isnan(xw_L);
plot(RH_L(ok)*100, xw_L(ok), '-', 'Color', col_LiCl,  'LineWidth', 2.5, ...
    'DisplayName', 'Pure LiCl');

ok = ~isnan(xw_B);
plot(RH_B(ok)*100, xw_B(ok), '-', 'Color', col_BaCl2, 'LineWidth', 2.5, ...
    'DisplayName', 'Pure BaCl_2');

ok = ~isnan(xw_M);
plot(RH_M(ok)*100, xw_M(ok), '-', 'Color', col_mix,   'LineWidth', 2.5, ...
    'DisplayName', 'LiCl + BaCl_2  (50:50 molar)');

xlabel('Relative Humidity (%)',    'FontSize', 14, 'FontWeight', 'bold');
ylabel('x_w  =  n_{H_2O} / (n_{H_2O} + n_{ions})', ...
    'FontSize', 14, 'FontWeight', 'bold');
title({'Molar Water Uptake vs. RH  —  LiCl, BaCl_2, and 50:50 Molar Mixture  (25 °C)', ...
       'LiCl: 2 ions/formula unit  |  BaCl_2: 3 ions/formula unit  |  mix: 5 ions per pair'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 13);
ylim([0 1]);
set(gca, 'FontSize', 13);

annotation('textbox', [0.56 0.17 0.36 0.12], 'FitBoxToText', 'off', ...
    'String', {'BaCl_2 limited to deliquescence range (0.91–0.96 RH)'; ...
               'Mix uses bivariate poly (calc\_activity\_LiCl\_BaCl2)'}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.6 0.6 0.6]);

out = fullfile(fig_out_dir, 'LiCl_BaCl2_molar_uptake');
saveas(fig3, [out '.png']);
savefig(fig3, [out '.fig']);
fprintf('   Saved: LiCl_BaCl2_molar_uptake.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  PLOT 4 – Mass-based water uptake  mf_w  vs  RH
% ═══════════════════════════════════════════════════════════════════════════
fprintf('Generating Plot 4: mass uptake vs RH ...\n');

fig4 = figure('Position', [140 140 1200 720], 'Color', 'w');
hold on; grid on; box on;

ok = ~isnan(mfw_L);
plot(RH_L(ok)*100, mfw_L(ok), '-', 'Color', col_LiCl,  'LineWidth', 2.5, ...
    'DisplayName', 'Pure LiCl');

ok = ~isnan(mfw_B);
plot(RH_B(ok)*100, mfw_B(ok), '-', 'Color', col_BaCl2, 'LineWidth', 2.5, ...
    'DisplayName', 'Pure BaCl_2');

ok = ~isnan(mfw_M);
plot(RH_M(ok)*100, mfw_M(ok), '-', 'Color', col_mix,   'LineWidth', 2.5, ...
    'DisplayName', 'LiCl + BaCl_2  (50:50 molar)');

xlabel('Relative Humidity (%)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mf_w  =  m_{H_2O} / (m_{H_2O} + m_{salt})', ...
    'FontSize', 14, 'FontWeight', 'bold');
title({'Mass-Based Water Uptake vs. RH  —  LiCl, BaCl_2, and 50:50 Molar Mixture  (25 °C)', ...
       'MW_{LiCl}=42.4 g/mol  |  MW_{BaCl_2}=208.2 g/mol  |  mix charged per pair'}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 13);
ylim([0 1]);
set(gca, 'FontSize', 13);

annotation('textbox', [0.56 0.17 0.36 0.12], 'FitBoxToText', 'off', ...
    'String', {'BaCl_2: high MW → less water by mass at same molality'; ...
               '50:50 mix blends both effects'}, ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.95], 'EdgeColor', [0.6 0.6 0.6]);

out = fullfile(fig_out_dir, 'LiCl_BaCl2_mass_uptake');
saveas(fig4, [out '.png']);
savefig(fig4, [out '.fig']);
fprintf('   Saved: LiCl_BaCl2_mass_uptake.png\n');

fprintf('\nAll four figures saved to:\n  %s\n', fig_out_dir);
fprintf('Done.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER FUNCTIONS
% ═══════════════════════════════════════════════════════════════════════════

function aw = mix_activity_50(m_each, MW_LiCl, MW_BaCl2)
% Water activity of 50:50 molar LiCl+BaCl₂ at equal molality m_each
% (m_LiCl = m_BaCl₂ = m_each).  Uses the bivariate polynomial
% calculate_activity_LiCl_BaCl2(mf1, mf2) where mf1=LiCl, mf2=BaCl₂.
    denom = 1000 + m_each*(MW_LiCl + MW_BaCl2);
    mf1   = m_each * MW_LiCl  / denom;
    mf2   = m_each * MW_BaCl2 / denom;
    aw    = calculate_activity_LiCl_BaCl2(mf1, mf2);
end
