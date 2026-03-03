function plot_mixture_3d_facets()
%PLOT_MIXTURE_3D_FACETS  3-D faceted comparison: polynomial fit vs measured data
%
%   For every salt mixture in mixture_sources_list.xlsx, draws:
%     - SURFACE  : water activity predicted by the bivariate polynomial fit
%                  (the calculate_activity_<Name> functions in this folder)
%     - POINTS   : measured (mf_salt1, mf_salt2, a_w) data from the Excel file
%
%   One individual figure is saved per mixture, plus one summary facet grid.
%
%   Outputs (written to figures/mixture_3d_facets/):
%     all_mixtures_3d_facets.png / .fig   –– summary grid, all mixtures
%     <FuncName>_3d.png / .fig            –– one figure per mixture

close all; clc;
warning('off', 'all');

%% ---- Paths ---------------------------------------------------------------
[filepath, ~, ~] = fileparts(mfilename('fullpath'));
addpath(filepath);  % make calculate_activity_* functions callable

excel_file  = fullfile(filepath, '..', 'data', 'mixture_sources_list.xlsx');
fig_out_dir = fullfile(filepath, '..', 'figures', 'mixture_3d_facets');
if ~exist(fig_out_dir, 'dir'), mkdir(fig_out_dir); end

%% ---- Read Excel ----------------------------------------------------------
fprintf('Reading %s ...\n', excel_file);
raw = readcell(excel_file);
[n_rows_raw, n_cols_raw] = size(raw);

%% ---- Locate mixture column blocks ----------------------------------------
% Row 1 contains headers like "Paper: NaCl + LiCl".  Each mixture occupies
% a block of up to 8 columns starting at that header column.
mix_col_starts = [];
mix_name_strs  = {};

for c = 1:n_cols_raw
    val = raw{1, c};
    if ischar(val) || isstring(val)
        val = char(val);
        if contains(val, 'Paper') && contains(val, ':')
            mix_col_starts(end+1) = c;                         %#ok<AGROW>
            colon_pos = strfind(val, ':');
            mix_name_strs{end+1} = strtrim(val(colon_pos(1)+1:end));  %#ok<AGROW>
        end
    end
end
fprintf('Found %d mixture blocks.\n\n', numel(mix_col_starts));

%% ---- Parse each mixture block --------------------------------------------
mix_data = [];   % struct array, grows below

for mi = 1:numel(mix_col_starts)
    c0       = mix_col_starts(mi);
    c_end    = min(c0 + 7, n_cols_raw);
    name_str = mix_name_strs{mi};

    % ---- split into component names ("Salt1 + Salt2") --------------------
    parts = strsplit(name_str, '+');
    if numel(parts) ~= 2, continue; end
    comp1 = strtrim(parts{1});
    comp2 = strtrim(parts{2});

    % ---- derive expected function name -----------------------------------
    func_name = build_func_name(name_str);

    % ---- locate header row (contains 'aw' or 'a_w') within rows 2..15 ---
    hdr_row = -1;
    for r = 2 : min(15, n_rows_raw)
        for c = c0:c_end
            v = raw{r, c};
            if ischar(v) || isstring(v)
                vs = lower(strtrim(char(v)));
                if strcmp(vs, 'aw') || strcmp(vs, 'a_w')
                    hdr_row = r;
                    break;
                end
            end
        end
        if hdr_row > 0, break; end
    end

    if hdr_row < 0
        fprintf('  [%2d] %-35s : header row not found – skipping.\n', mi, name_str);
        continue;
    end

    % ---- identify aw column and mass-fraction columns --------------------
    aw_col  = -1;
    mf_cols = [];
    for c = c0:c_end
        v = raw{hdr_row, c};
        if ischar(v) || isstring(v)
            vs = lower(strtrim(char(v)));
            if strcmp(vs, 'aw') || strcmp(vs, 'a_w')
                aw_col = c;
            elseif contains(vs, 'mass fraction')
                mf_cols(end+1) = c;  %#ok<AGROW>
            end
        end
    end

    if aw_col < 0 || numel(mf_cols) < 2
        fprintf('  [%2d] %-35s : aw/mf columns not found – skipping.\n', mi, name_str);
        continue;
    end

    % ---- extract numeric data rows (stop at first missing/non-numeric) ---
    mf1 = [];  mf2 = [];  aw = [];
    for r = hdr_row+1 : n_rows_raw
        v_aw = raw{r, aw_col};
        if ~isnumeric(v_aw) || ~isscalar(v_aw) || ~isfinite(v_aw)
            break;   % end of data block
        end
        v1 = raw{r, mf_cols(1)};
        v2 = raw{r, mf_cols(2)};
        if ~isnumeric(v1) || ~isscalar(v1) || ~isfinite(v1), v1 = 0; end
        if ~isnumeric(v2) || ~isscalar(v2) || ~isfinite(v2), v2 = 0; end
        mf1(end+1) = v1;    %#ok<AGROW>
        mf2(end+1) = v2;    %#ok<AGROW>
        aw(end+1)  = v_aw;  %#ok<AGROW>
    end

    if isempty(aw)
        fprintf('  [%2d] %-35s : no data rows found – skipping.\n', mi, name_str);
        continue;
    end

    fprintf('  [%2d] %-35s : %d points\n', mi, name_str, numel(aw));

    entry.name      = func_name;
    entry.comp1     = comp1;
    entry.comp2     = comp2;
    entry.func_name = func_name;
    entry.mf1       = mf1(:);
    entry.mf2       = mf2(:);
    entry.aw        = aw(:);

    if isempty(mix_data)
        mix_data = entry;
    else
        mix_data(end+1) = entry;  %#ok<AGROW>
    end
end

n_mix = numel(mix_data);
fprintf('\nSuccessfully parsed %d mixtures.\n\n', n_mix);

if n_mix == 0
    warning('No mixture data found – check Excel file path / format.');
    return;
end

%% ---- Build summary facet figure ------------------------------------------
n_cols_fig = 4;
n_rows_fig = ceil(n_mix / n_cols_fig);
fig_w = n_cols_fig * 370;
fig_h = n_rows_fig * 330;

fig_all = figure('Position', [20 20 fig_w fig_h], 'Color', 'w', 'Visible', 'off');

%% ---- Loop: individual figures + summary panels ---------------------------
fprintf('Generating figures...\n');

for mi = 1:n_mix
    md = mix_data(mi);

    % Check the .m function file is available
    if exist(md.func_name, 'file') ~= 2
        fprintf('  WARNING: %s.m not found – skipping plot.\n', md.func_name);
        continue;
    end

    % ---- individual 3-D figure -------------------------------------------
    fig_ind = figure('Position', [50 50 660 540], 'Color', 'w', 'Visible', 'off');
    ax_ind  = axes(fig_ind);
    draw_3d(ax_ind, md, true);

    png_ind = fullfile(fig_out_dir, [md.name '_3d.png']);
    fig_ind_path = fullfile(fig_out_dir, [md.name '_3d.fig']);
    exportgraphics(fig_ind, png_ind, 'Resolution', 150);
    savefig(fig_ind, fig_ind_path);
    close(fig_ind);
    fprintf('  Saved: %s_3d.png\n', md.name);

    % ---- panel in summary figure -----------------------------------------
    ax_sum = subplot(n_rows_fig, n_cols_fig, mi, 'Parent', fig_all);
    draw_3d(ax_sum, md, false);
end

sgtitle(fig_all, ...
    'Water Activity: Polynomial Fit (surface) vs Measured Data (black dots)', ...
    'FontSize', 12, 'FontWeight', 'bold');

out_png = fullfile(fig_out_dir, 'all_mixtures_3d_facets.png');
out_fig = fullfile(fig_out_dir, 'all_mixtures_3d_facets.fig');
exportgraphics(fig_all, out_png, 'Resolution', 150);
savefig(fig_all, out_fig);
close(fig_all);

fprintf('\nSummary figure saved:\n  %s\n', out_png);
fprintf('Individual figures saved to:\n  %s\n', fig_out_dir);

warning('on', 'all');
fprintf('\nDone.\n');

end  % --- main function ---


%% =========================================================================
function draw_3d(ax, md, show_full_title)
%DRAW_3D  Render one mixture's 3-D surface + scatter of measured data.
%
%   SURFACE  – polynomial fit evaluated on a meshgrid of (mf1, mf2)
%   SCATTER  – actual (mf1, mf2, aw) data points (black filled dots)

    mf1_data = md.mf1(:);
    mf2_data = md.mf2(:);
    aw_data  = md.aw(:);

    % Grid extents: slightly beyond the data maximum
    mf1_hi = min(max(mf1_data) * 1.15, 0.98);
    mf2_hi = min(max(mf2_data) * 1.15, 0.98);

    % Build surface
    n_g  = 35;
    mf1v = linspace(0, mf1_hi, n_g);
    mf2v = linspace(0, mf2_hi, n_g);
    [MF1, MF2] = meshgrid(mf1v, mf2v);
    AW = NaN(size(MF1));

    for ii = 1:numel(MF1)
        if MF1(ii) + MF2(ii) < 0.99
            try
                v = feval(md.func_name, MF1(ii), MF2(ii));
                if isfinite(v) && v >= 0 && v <= 1.05
                    AW(ii) = v;
                end
            catch
                % leave as NaN
            end
        end
    end

    % Surface
    h_surf = surf(ax, MF1, MF2, AW, ...
        'FaceAlpha', 0.60, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap(ax, parula);
    hold(ax, 'on');

    % Data points (black filled dots with light edge for visibility)
    scatter3(ax, mf1_data, mf2_data, aw_data, ...
        35, 'k', 'filled', ...
        'MarkerEdgeColor', [0.85 0.85 0.85], 'LineWidth', 0.5);

    % Axis labels
    xlabel(ax, sprintf('mf  %s', md.comp1), 'FontSize', 7, 'Interpreter', 'none');
    ylabel(ax, sprintf('mf  %s', md.comp2), 'FontSize', 7, 'Interpreter', 'none');
    zlabel(ax, 'a_w', 'FontSize', 7);

    % Z limits: hug the data range
    z_lo = max(0,    min(aw_data) * 0.97);
    z_hi = min(1.01, max(aw_data) * 1.02);
    zlim(ax, [z_lo, z_hi]);

    % Title
    if show_full_title
        title(ax, sprintf('%s  +  %s', md.comp1, md.comp2), ...
            'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
        subtitle(ax, sprintf('%d data points  |  RMSE printed in console', ...
            numel(aw_data)), 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    else
        title(ax, sprintf('%s + %s', md.comp1, md.comp2), ...
            'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none');
    end

    view(ax, 35, 28);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 7);

    % Print per-point RMSE in console for the individual figure call
    if show_full_title
        aw_pred = NaN(size(aw_data));
        for ii = 1:numel(aw_data)
            try
                aw_pred(ii) = feval(md.func_name, mf1_data(ii), mf2_data(ii));
            catch
            end
        end
        valid = isfinite(aw_pred);
        if any(valid)
            rmse = sqrt(mean((aw_pred(valid) - aw_data(valid)).^2));
            fprintf('    RMSE  %-35s : %.5f\n', [md.comp1 ' + ' md.comp2], rmse);
        end
    end
end  % draw_3d


%% =========================================================================
function fn = build_func_name(name_str)
%BUILD_FUNC_NAME  "NaCl + LiCl"  -->  "calculate_activity_NaCl_LiCl"
%  Mirrors the sanitize_filename logic in generate_mixture_activity_files.py:
%    1. Remove characters that are neither word chars, spaces, nor '+'
%    2. Replace runs of spaces / '+' with a single '_'
    s = strtrim(name_str);
    s = regexprep(s, '[^\w\s+]', '');   % strip parens, dots, etc.
    s = regexprep(s, '[\s+]+',   '_');  % spaces and + → underscore
    fn = ['calculate_activity_' s];
end
