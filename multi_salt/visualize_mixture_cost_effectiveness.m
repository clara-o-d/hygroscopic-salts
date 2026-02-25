close all
clear
clc

% Script: visualize_mixture_cost_effectiveness.m
%
% Compares the cost-effectiveness of water uptake for all binary salt
% mixtures AND the individual pure component salts at 25°C.
%
% Metric: grams of water held per $ of dry salt charged
%   cost_per_kg_water ($) = 0.5 * m_total * (MW1*c1 + MW2*c2)   [mixtures]
%                         = m * MW * cost_per_g                   [pure salts]
%   water_per_dollar (g/$) = 1000 / cost_per_kg_water
%
% Methods:
%   Mixtures  – forward sweep over m_total (no inversion, equal mole frac)
%   Pure salts – calculate_mf_* functions called over valid RH range
%
% Colour schemes (both colorblind-safe, Paul Tol):
%   Mixtures   – "bright"  palette, solid/dashed lines, filled/open markers
%   Pure salts – "muted"   palette, dash-dot lines,    diamond markers
%
% Plots:
%   1. All mixtures + pure salts, log y-scale, direct curve labels
%   2. Affordable subset (no CsCl), linear y-scale, direct curve labels
%   3. 2×2 panel – cost-band sensitivity for 4 key pairs

%% ── Paths ────────────────────────────────────────────────────────────────────
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_mixtures'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));

fig_out_dir = fullfile(filepath, '..', 'figures', 'multi_salt');
if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end

%% ── Colour palettes (both from Paul Tol, colorblind-safe) ───────────────────
% "bright" — for mixtures
CB_bright = [
     68 119 170;   % 1  blue
    102 204 238;   % 2  cyan
     34 136  51;   % 3  green
    204 187  68;   % 4  yellow
    238 102 119;   % 5  rose
    170  51 119;   % 6  purple
    187 187 187;   % 7  grey
] / 255;

% "muted" — for pure salts
CB_muted = [
     51  34 136;   % 1  indigo
    136 204 238;   % 2  light blue
     68 170 153;   % 3  teal
     17 119  51;   % 4  green
    153 153  51;   % 5  olive
    221 204 119;   % 6  sand
    204 102 119;   % 7  rose
    136  34  85;   % 8  wine
] / 255;

% 13 mixtures → 7 solid lines then 6 dashed (bright colours cycle)
n_mix = 13;
mix_colors = [CB_bright; CB_bright(1:6,:)];   % 13×3
mix_styles  = [repmat({'-'},7,1); repmat({'--'},6,1)];
mix_lw      = 2.5;
mix_markers = {'o','s','^','d','v','p','h','o','s','^','d','v','p'};
mix_filled  = [true(7,1); false(6,1)];

% 8 pure salts — dash-dot lines, diamond markers, muted colours
pure_lw = 2.0;

%% ── Salt property maps ───────────────────────────────────────────────────────
MW_map = containers.Map( ...
    {'LiCl','NaCl','CaCl2','MgCl2','BaCl2','KCl','CsCl', ...
     'NH4Cl','NaNO3','LiNO3','NaAc','LiAc','LiClO4','NaClO4'}, ...
    {42.394, 58.443, 110.983, 95.211, 208.233, 74.551, 168.358, ...
     53.491, 84.994, 68.946, 82.034, 65.989, 106.392, 122.440});

NU_map = containers.Map( ...
    {'LiCl','NaCl','CaCl2','MgCl2','BaCl2','KCl','CsCl', ...
     'NH4Cl','NaNO3','LiNO3','NaAc','LiAc','LiClO4','NaClO4'}, ...
    {2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2});

C_map = containers.Map( ...
    {'LiCl','NaCl','CaCl2','MgCl2','BaCl2','KCl','CsCl', ...
     'NH4Cl','NaNO3','LiNO3','NaAc','LiAc','LiClO4','NaClO4'}, ...
    {[0.25 0.375 0.50], [0.01 0.02 0.03], [0.005 0.013 0.02], ...
     [0.003 0.007 0.01], [0.05 0.10 0.15], [0.01 0.02 0.03], ...
     [5.00 12.5 20.0], [0.01 0.03 0.05], [0.005 0.013 0.02], ...
     [0.10 0.20 0.30], [0.01 0.03 0.05], [0.10 0.20 0.30], ...
     [0.30 0.65 1.00], [0.05 0.125 0.20]});

%% ── Mixture definitions ──────────────────────────────────────────────────────
mix_defs = { ...
    @calculate_activity_MgCl2_NaCl,         'MgCl_2 + NaCl',        'MgCl2',  'NaCl',   0.009431,0.103303,0.008690,0.174613; ...
    @calculate_activity_MgCl2_CaCl2,        'MgCl_2 + CaCl_2',      'MgCl2',  'CaCl2',  0.009431,0.222171,0.010976,0.249779; ...
    @calculate_activity_NH4Cl_LiCl,         'NH_4Cl + LiCl',         'NH4Cl',  'LiCl',   0.005321,0.176250,0.004222,0.144989; ...
    @calculate_activity_NaNO3_LiNO3,        'NaNO_3 + LiNO_3',       'NaNO3',  'LiNO3',  0.017575,0.359893,0.017757,0.247139; ...
    @calculate_activity_NaCl_LiCl,          'NaCl + LiCl',           'NaCl',   'LiCl',   0.005745,0.236635,0.004704,0.154883; ...
    @calculate_activity_NaC2H3O2_LiC2H3O2,  'NaAc + LiAc',           'NaAc',   'LiAc',   0.018404,0.243676,0.014785,0.224370; ...
    @calculate_activity_LiCl_KCl,           'LiCl + KCl',            'LiCl',   'KCl',    0.004222,0.144989,0.007400,0.229703; ...
    @calculate_activity_LiCl_MgCl2,         'LiCl + MgCl_2',         'LiCl',   'MgCl2',  0.007990,0.276173,0.004738,0.186003; ...
    @calculate_activity_LiCl_CaCl2,         'LiCl + CaCl_2',         'LiCl',   'CaCl2',  0.003801,0.385054,0.001109,0.387891; ...
    @calculate_activity_LiCl_BaCl2,         'LiCl + BaCl_2',         'LiCl',   'BaCl2',  0.005062,0.161999,0.010304,0.210444; ...
    @calculate_activity_LiClO4_NaClO4,      'LiClO_4 + NaClO_4',     'LiClO4', 'NaClO4', 0.020940,0.272325,0.023213,0.350219; ...
    @calculate_activity_NaCl_CsCl,          'NaCl + CsCl',           'NaCl',   'CsCl',   0.005810,0.237506,0.032575,0.472952; ...
    @calculate_activity_LiCl_CsCl,          'LiCl + CsCl',           'LiCl',   'CsCl',   0.004222,0.202783,0.032575,0.502527; ...
};

%% ── Pure salt definitions ────────────────────────────────────────────────────
% {display_label, calculate_mf_suffix, MW_key, cost_key,
%  RH_lo, RH_hi, T_arg ([] = no T argument)}
% RH bounds are slightly inside the valid range to avoid edge errors.
pure_defs = { ...
    'MgCl_2',  'MgCl',   'MgCl2',  'MgCl2',  0.34, 0.985, []; ...
    'CaCl_2',  'CaCl',   'CaCl2',  'CaCl2',  0.10, 0.985, 25; ...
    'NaCl',    'NaCl',   'NaCl',   'NaCl',   0.77, 0.990, []; ...
    'KCl',     'KCl',    'KCl',    'KCl',    0.86, 0.990, []; ...
    'NH_4Cl',  'NH4Cl',  'NH4Cl',  'NH4Cl',  0.82, 0.990, []; ...
    'LiNO_3',  'LiNO3',  'LiNO3',  'LiNO3',  0.74, 0.994, []; ...
    'LiCl',    'LiCl',   'LiCl',   'LiCl',   0.12, 0.985, 25; ...
    'CsCl',    'CsCl',   'CsCl',   'CsCl',   0.82, 0.990, []; ...
};
n_pure = size(pure_defs, 1);

% Label x-positions for pure salts (chosen to avoid overlapping mixture labels)
pure_label_rh = [50, 55, 92, 94, 96, 77, 40, 93];

%% ── Compute mixture curves ───────────────────────────────────────────────────
m_sweep = linspace(0.05, 30, 800);
n_w0    = 1000 / 18.015;

fprintf('Computing mixture curves...\n');
mix = struct('aw',[],'wpd_mid',[],'wpd_best',[],'wpd_worst',[], ...
             'label','','s1','','s2','','c1',[],'c2',[]);

for k = 1:n_mix
    func    = mix_defs{k,1};   label   = mix_defs{k,2};
    s1      = mix_defs{k,3};   s2      = mix_defs{k,4};
    mf1_min = mix_defs{k,5};   mf1_max = mix_defs{k,6};
    mf2_min = mix_defs{k,7};   mf2_max = mix_defs{k,8};

    mw1 = MW_map(s1);  mw2 = MW_map(s2);
    c1  = C_map(s1);   c2  = C_map(s2);

    denom = 0.5*m_sweep*mw1 + 0.5*m_sweep*mw2 + 1000;
    mf1   = 0.5*m_sweep.*mw1 ./ denom;
    mf2   = 0.5*m_sweep.*mw2 ./ denom;

    valid = mf1>=mf1_min & mf1<=mf1_max & mf2>=mf2_min & mf2<=mf2_max;
    if ~any(valid)
        fprintf('  WARNING: no valid range for %s\n', strrep(label,'_',''));
        mix(k).label=label; mix(k).s1=s1; mix(k).s2=s2; mix(k).c1=c1; mix(k).c2=c2;
        continue
    end

    mv=m_sweep(valid);  mf1v=mf1(valid);  mf2v=mf2(valid);

    warning('off','all');
    aw_v = arrayfun(@(a,b) func(a,b), mf1v, mf2v);
    warning('on','all');

    ok = aw_v>0 & aw_v<=1 & ~isnan(aw_v);
    aw_v=aw_v(ok);  mv=mv(ok);
    [aw_v,srt]=sort(aw_v);  mv=mv(srt);

    mix(k).aw        = aw_v;
    mix(k).wpd_mid   = 1000 ./ (0.5*mv.*(mw1*c1(2)+mw2*c2(2)));
    mix(k).wpd_best  = 1000 ./ (0.5*mv.*(mw1*c1(1)+mw2*c2(1)));
    mix(k).wpd_worst = 1000 ./ (0.5*mv.*(mw1*c1(3)+mw2*c2(3)));
    mix(k).label=label; mix(k).s1=s1; mix(k).s2=s2; mix(k).c1=c1; mix(k).c2=c2;

    fprintf('  [mix] %-32s aw %.2f–%.2f  wpd %.0f–%.0f g/$\n', ...
        strrep(label,'_',''), min(aw_v), max(aw_v), ...
        min(mix(k).wpd_mid), max(mix(k).wpd_mid));
end

%% ── Compute pure salt curves ─────────────────────────────────────────────────
fprintf('Computing pure salt curves...\n');
pure = struct('aw',[],'wpd_mid',[],'wpd_best',[],'wpd_worst',[],'label','','c',[]);

for k = 1:n_pure
    label   = pure_defs{k,1};
    suf     = pure_defs{k,2};
    mw_key  = pure_defs{k,3};
    c_key   = pure_defs{k,4};
    rh_lo   = pure_defs{k,5};
    rh_hi   = pure_defs{k,6};
    T_arg   = pure_defs{k,7};

    MW_s = MW_map(mw_key);
    c    = C_map(c_key);          % [lo, mid, hi]
    func_h = str2func(['calculate_mf_' suf]);

    RH_vec = linspace(rh_lo, rh_hi, 300);
    mf_v   = nan(size(RH_vec));

    for j = 1:length(RH_vec)
        try
            if isempty(T_arg)
                mf_v(j) = func_h(RH_vec(j));
            else
                mf_v(j) = func_h(RH_vec(j), T_arg);
            end
        catch
        end
    end

    ok  = ~isnan(mf_v) & mf_v > 0 & mf_v < 1;
    rh_ok = RH_vec(ok);
    mf_ok = mf_v(ok);

    % cost_per_kg_water = m * MW * cost_per_g
    %   where m = 1000*mf / (MW*(1-mf))
    %   → cost = 1000*mf*cost_per_g / (1-mf)
    denom_cost = 1 - mf_ok;
    pure(k).aw        = rh_ok;
    pure(k).wpd_mid   = denom_cost ./ (mf_ok * c(2));
    pure(k).wpd_best  = denom_cost ./ (mf_ok * c(1));
    pure(k).wpd_worst = denom_cost ./ (mf_ok * c(3));
    pure(k).label = label;
    pure(k).c     = c;

    fprintf('  [pure] %-12s  aw %.2f–%.2f  wpd %.0f–%.0f g/$\n', ...
        strrep(label,'_',''), min(rh_ok), max(rh_ok), ...
        min(pure(k).wpd_mid), max(pure(k).wpd_mid));
end
fprintf('\n');

% Label placement for mixture lines (staggered RH positions)
mix_label_rh = [88, 86, 84, 82, 80, 87, 85, 83, 81, 79, 78, 76, 74];

%% ── Plot 1: All mixtures + pure salts, log y-scale ───────────────────────────
cost_str = build_cost_str(C_map);

fig1 = figure('Position', [30 30 1650 1000]);
ax1  = axes(fig1, 'Position', [0.06 0.13 0.68 0.80]);
hold(ax1,'on');  grid(ax1,'on');  box(ax1,'on');

% ── Mixture curves ────────────────────────────────────────────────────────────
for k = 1:n_mix
    d   = mix(k);
    col = mix_colors(k,:);
    lst = mix_styles{k};
    mrk = mix_markers{k};
    if isempty(d.aw), continue; end

    rh  = d.aw * 100;
    wpd = d.wpd_mid;

    % Shaded cost band
    rh_p  = [rh, fliplr(rh)];
    wpd_p = [d.wpd_best, fliplr(d.wpd_worst)];
    if ~any(isnan(wpd_p))
        fill(ax1, rh_p, wpd_p, col, 'FaceAlpha', 0.10, ...
             'EdgeColor','none', 'HandleVisibility','off');
    end

    plot(ax1, rh, wpd, lst, 'LineWidth', mix_lw, 'Color', col, ...
         'HandleVisibility','off');

    n  = length(rh);
    mi = unique(round(linspace(max(1,round(n*0.07)), n, 5)));
    mi = mi(mi <= n);
    if mix_filled(k)
        plot(ax1, rh(mi), wpd(mi), mrk, 'MarkerSize',6, ...
             'Color',col, 'MarkerFaceColor',col, 'HandleVisibility','off');
    else
        plot(ax1, rh(mi), wpd(mi), mrk, 'MarkerSize',6, ...
             'Color',col, 'MarkerFaceColor','w', 'HandleVisibility','off');
    end

    % Direct label
    tgt = mix_label_rh(k);
    [~,ti] = min(abs(rh - tgt));
    ti = min(ti, length(rh));
    text(ax1, rh(ti), wpd(ti), [' ' strrep(d.label,'_','')], ...
         'Color', min(col*0.80, ones(1,3)), 'FontSize', 8, 'FontWeight','bold', ...
         'VerticalAlignment','middle', 'Interpreter','tex');
end

% ── Pure salt curves ──────────────────────────────────────────────────────────
for k = 1:n_pure
    d   = pure(k);
    col = CB_muted(k,:);
    if isempty(d.aw), continue; end

    rh  = d.aw * 100;
    wpd = d.wpd_mid;

    % Shaded cost band
    rh_p  = [rh, fliplr(rh)];
    wpd_p = [d.wpd_best, fliplr(d.wpd_worst)];
    if ~any(isnan(wpd_p))
        fill(ax1, rh_p, wpd_p, col, 'FaceAlpha', 0.13, ...
             'EdgeColor','none', 'HandleVisibility','off');
    end

    plot(ax1, rh, wpd, '-.', 'LineWidth', pure_lw, 'Color', col, ...
         'HandleVisibility','off');

    n  = length(rh);
    mi = unique(round(linspace(max(1,round(n*0.07)), n, 5)));
    mi = mi(mi <= n);
    plot(ax1, rh(mi), wpd(mi), 'd', 'MarkerSize', 7, ...
         'Color', col, 'MarkerFaceColor', col, 'HandleVisibility','off');

    % Direct label — placed at pure_label_rh target
    tgt = pure_label_rh(k);
    [~,ti] = min(abs(rh - tgt));
    ti = min(ti, length(rh));
    text(ax1, rh(ti), wpd(ti), [' ' strrep(d.label,'_','')], ...
         'Color', min(col*0.80, ones(1,3)), 'FontSize', 8.5, 'FontWeight','bold', ...
         'VerticalAlignment','middle', 'Interpreter','tex', ...
         'BackgroundColor',[1 1 1 0.6], 'Margin', 1);
end

set(ax1, 'YScale','log', 'FontSize',13);
xlabel(ax1, 'Equilibrium Relative Humidity  (%)', 'FontSize',15,'FontWeight','bold');
ylabel(ax1, 'Water Uptake per Dollar   (g H_2O / $)', 'FontSize',15,'FontWeight','bold');
title(ax1, {'Cost-Effectiveness of Salt Mixtures and Pure Salts  (25°C)', ...
    'Shaded bands = lo–hi price range   ·   Mixtures at equal mole fractions'}, ...
    'FontSize', 14, 'FontWeight','bold');

% ── Legend: dummy lines for mixtures, then pure salts ────────────────────────
% Mixtures section
h_mix = gobjects(n_mix,1);
for k = 1:n_mix
    d = mix(k);
    if isempty(d.aw), continue; end
    col = mix_colors(k,:);
    h_mix(k) = plot(ax1, NaN, NaN, [mix_styles{k} mix_markers{k}], ...
                    'LineWidth',2.0, 'Color',col, 'MarkerSize',6, ...
                    'DisplayName', strrep(d.label,' ',' '));
    if mix_filled(k), h_mix(k).MarkerFaceColor = col;
    else,             h_mix(k).MarkerFaceColor = 'w'; end
end

% Pure salts section — labelled with (pure) to distinguish
h_pure = gobjects(n_pure,1);
for k = 1:n_pure
    d = pure(k);
    if isempty(d.aw), continue; end
    col = CB_muted(k,:);
    h_pure(k) = plot(ax1, NaN, NaN, '-. d', 'LineWidth',2.0, 'Color',col, ...
                     'MarkerSize',7, 'MarkerFaceColor',col, ...
                     'DisplayName', [strrep(d.label,' ',' ') ' (pure)']);
end

% Build legend with section headers using invisible placeholder entries
h_sep1 = plot(ax1, NaN,NaN,'w.','DisplayName','— Mixtures (equal mol fractions) —');
h_sep2 = plot(ax1, NaN,NaN,'w.','DisplayName','— Pure salts —');

valid_mix  = h_mix(arrayfun(@(h) ishandle(h) && h~=0, h_mix));
valid_pure = h_pure(arrayfun(@(h) ishandle(h) && h~=0, h_pure));
leg1 = legend(ax1, [h_sep1; valid_mix; h_sep2; valid_pure], ...
              'Location','eastoutside', 'FontSize',8.5, 'Interpreter','tex');
leg1.Title.String = 'Salt / mixture';
leg1.Title.FontSize = 9;

% Cost annotation
annotation(fig1,'textbox',[0.06 0.005 0.88 0.060],'String',cost_str, ...
    'FontSize',7,'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.7 0.7 0.7], ...
    'FitBoxToText','off','HorizontalAlignment','left');

set(fig1,'color','w');
saveas(fig1, fullfile(fig_out_dir,'mixture_cost_effectiveness_all.png'));
savefig(fig1, fullfile(fig_out_dir,'mixture_cost_effectiveness_all.fig'));
fprintf('Saved: mixture_cost_effectiveness_all.png\n');

%% ── Plot 2: Exclude CsCl, linear y-scale ────────────────────────────────────
no_CsCl_mix  = ~cellfun(@(s) contains(s,'CsCl'), mix_defs(:,3)) & ...
               ~cellfun(@(s) contains(s,'CsCl'), mix_defs(:,4));
no_CsCl_pure = ~cellfun(@(s) contains(s,'CsCl'), pure_defs(:,4));

fig2 = figure('Position', [30 30 1650 1000]);
ax2  = axes(fig2, 'Position', [0.06 0.13 0.68 0.80]);
hold(ax2,'on');  grid(ax2,'on');  box(ax2,'on');

for k = 1:n_mix
    if ~no_CsCl_mix(k), continue; end
    d=mix(k); col=mix_colors(k,:); lst=mix_styles{k}; mrk=mix_markers{k};
    if isempty(d.aw), continue; end
    rh=d.aw*100; wpd=d.wpd_mid;
    rh_p=[rh,fliplr(rh)]; wpd_p=[d.wpd_best,fliplr(d.wpd_worst)];
    if ~any(isnan(wpd_p))
        fill(ax2,rh_p,wpd_p,col,'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax2,rh,wpd,lst,'LineWidth',mix_lw,'Color',col,'HandleVisibility','off');
    n=length(rh); mi=unique(round(linspace(max(1,round(n*0.07)),n,5))); mi=mi(mi<=n);
    if mix_filled(k), plot(ax2,rh(mi),wpd(mi),mrk,'MarkerSize',6,'Color',col,'MarkerFaceColor',col,'HandleVisibility','off');
    else,             plot(ax2,rh(mi),wpd(mi),mrk,'MarkerSize',6,'Color',col,'MarkerFaceColor','w','HandleVisibility','off'); end
    tgt=mix_label_rh(k); [~,ti]=min(abs(rh-tgt)); ti=min(ti,length(rh));
    text(ax2,rh(ti),wpd(ti),[' ' strrep(d.label,'_','')],'Color',min(col*0.80,ones(1,3)), ...
         'FontSize',8,'FontWeight','bold','VerticalAlignment','middle','Interpreter','tex');
end

for k = 1:n_pure
    if ~no_CsCl_pure(k), continue; end
    d=pure(k); col=CB_muted(k,:);
    if isempty(d.aw), continue; end
    rh=d.aw*100; wpd=d.wpd_mid;
    rh_p=[rh,fliplr(rh)]; wpd_p=[d.wpd_best,fliplr(d.wpd_worst)];
    if ~any(isnan(wpd_p))
        fill(ax2,rh_p,wpd_p,col,'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax2,rh,wpd,'-.','LineWidth',pure_lw,'Color',col,'HandleVisibility','off');
    n=length(rh); mi=unique(round(linspace(max(1,round(n*0.07)),n,5))); mi=mi(mi<=n);
    plot(ax2,rh(mi),wpd(mi),'d','MarkerSize',7,'Color',col,'MarkerFaceColor',col,'HandleVisibility','off');
    tgt=pure_label_rh(k); [~,ti]=min(abs(rh-tgt)); ti=min(ti,length(rh));
    text(ax2,rh(ti),wpd(ti),[' ' strrep(d.label,'_','')],'Color',min(col*0.80,ones(1,3)), ...
         'FontSize',8.5,'FontWeight','bold','VerticalAlignment','middle','Interpreter','tex', ...
         'BackgroundColor',[1 1 1 0.6],'Margin',1);
end

set(ax2,'FontSize',13);
xlabel(ax2,'Equilibrium Relative Humidity  (%)','FontSize',15,'FontWeight','bold');
ylabel(ax2,'Water Uptake per Dollar   (g H_2O / $)','FontSize',15,'FontWeight','bold');
title(ax2,{'Cost-Effectiveness of Salt Mixtures and Pure Salts — Excluding CsCl  (25°C)', ...
    'Shaded bands = lo–hi price range   ·   Mixtures at equal mole fractions'}, ...
    'FontSize',14,'FontWeight','bold');

% Legend
h_m2 = gobjects(n_mix,1);
for k=1:n_mix
    if ~no_CsCl_mix(k)||isempty(mix(k).aw), continue; end
    col=mix_colors(k,:);
    h_m2(k)=plot(ax2,NaN,NaN,[mix_styles{k} mix_markers{k}],'LineWidth',2.0,'Color',col, ...
                 'MarkerSize',6,'DisplayName',strrep(mix(k).label,' ',' '));
    if mix_filled(k), h_m2(k).MarkerFaceColor=col;
    else,             h_m2(k).MarkerFaceColor='w'; end
end
h_p2 = gobjects(n_pure,1);
for k=1:n_pure
    if ~no_CsCl_pure(k)||isempty(pure(k).aw), continue; end
    col=CB_muted(k,:);
    h_p2(k)=plot(ax2,NaN,NaN,'-.d','LineWidth',2.0,'Color',col,'MarkerSize',7, ...
                 'MarkerFaceColor',col,'DisplayName',[strrep(pure(k).label,' ',' ') ' (pure)']);
end
hs1=plot(ax2,NaN,NaN,'w.','DisplayName','— Mixtures (equal mol fractions) —');
hs2=plot(ax2,NaN,NaN,'w.','DisplayName','— Pure salts —');
vm2=h_m2(arrayfun(@(h)ishandle(h)&&h~=0,h_m2));
vp2=h_p2(arrayfun(@(h)ishandle(h)&&h~=0,h_p2));
leg2=legend(ax2,[hs1;vm2;hs2;vp2],'Location','eastoutside','FontSize',8.5,'Interpreter','tex');
leg2.Title.String='Salt / mixture'; leg2.Title.FontSize=9;

annotation(fig2,'textbox',[0.06 0.005 0.88 0.060],'String',cost_str, ...
    'FontSize',7,'BackgroundColor',[0.97 0.97 0.97],'EdgeColor',[0.7 0.7 0.7], ...
    'FitBoxToText','off','HorizontalAlignment','left');

set(fig2,'color','w');
saveas(fig2,fullfile(fig_out_dir,'mixture_cost_effectiveness_noCsCl.png'));
savefig(fig2,fullfile(fig_out_dir,'mixture_cost_effectiveness_noCsCl.fig'));
fprintf('Saved: mixture_cost_effectiveness_noCsCl.png\n');

%% ── Plot 3: 2×2 panel – cost-band sensitivity for 4 key pairs ───────────────
key_pairs = [1, 2, 5, 8];
panel_pos = {[0.07 0.55 0.42 0.39],[0.55 0.55 0.42 0.39], ...
             [0.07 0.08 0.42 0.39],[0.55 0.08 0.42 0.39]};

fig3 = figure('Position', [30 30 1400 980]);
for p = 1:4
    k=key_pairs(p); d=mix(k); col=mix_colors(k,:);
    ax=axes(fig3,'Position',panel_pos{p}); %#ok<LAXES>
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    if ~isempty(d.aw)
        rh=d.aw*100;
        fill(ax,[rh fliplr(rh)],[d.wpd_best fliplr(d.wpd_worst)],col, ...
             'FaceAlpha',0.25,'EdgeColor',col*0.7,'LineStyle',':', ...
             'LineWidth',0.8,'DisplayName','Cost range');
        plot(ax,rh,d.wpd_mid,  '-', 'LineWidth',3,  'Color',col,    'DisplayName','Midpoint cost');
        plot(ax,rh,d.wpd_best, '--','LineWidth',1.4,'Color',col*0.6,'DisplayName','Best (lo prices)');
        plot(ax,rh,d.wpd_worst,':' ,'LineWidth',1.8,'Color',col*0.6,'DisplayName','Worst (hi prices)');
    end
    xlabel(ax,'RH  (%)','FontSize',12,'FontWeight','bold');
    ylabel(ax,'g H_2O / $','FontSize',12,'FontWeight','bold');
    title(ax,strrep(d.label,'_',' '),'FontSize',13,'FontWeight','bold', ...
          'Color',min(col*0.65,ones(1,3)));
    legend(ax,'Location','northwest','FontSize',8);
    set(ax,'FontSize',11);
end
sgtitle(fig3,{'Water Uptake per Dollar — Price-Sensitivity Panels  (25°C, equal mole fractions)', ...
    cost_str},'FontSize',12,'FontWeight','bold');
set(fig3,'color','w');
saveas(fig3,fullfile(fig_out_dir,'mixture_cost_bands_panels.png'));
savefig(fig3,fullfile(fig_out_dir,'mixture_cost_bands_panels.fig'));
fprintf('Saved: mixture_cost_bands_panels.png\n');

fprintf('\nAll figures saved to:\n%s\n', fig_out_dir);
fprintf('Done.\n');

%% ══ Local functions (must be at end of script) ═══════════════════════════════

function s = build_cost_str(C_map)
    pairs = { ...
        'MgCl_2','MgCl2'; 'NaCl','NaCl'; 'CaCl_2','CaCl2'; 'KCl','KCl'; ...
        'NH_4Cl','NH4Cl'; 'LiCl','LiCl'; 'LiNO_3','LiNO3'; 'NaNO_3','NaNO3'; ...
        'NaAc','NaAc'; 'LiAc','LiAc'; 'BaCl_2','BaCl2'; ...
        'LiClO_4','LiClO4'; 'NaClO_4','NaClO4'; 'CsCl','CsCl'};
    parts = cell(size(pairs,1),1);
    for i = 1:size(pairs,1)
        c = C_map(pairs{i,2});
        parts{i} = sprintf('%s: $%.3g–%.3g/g', pairs{i,1}, c(1), c(3));
    end
    s = ['Salt costs:   ' strjoin(parts,'    ')];
end
