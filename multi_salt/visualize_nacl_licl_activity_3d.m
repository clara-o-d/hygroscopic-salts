close all
clear
clc

% Script to visualize water activity for NaCl + LiCl mixtures as a 3D surface
% Shows how water activity varies with different molalities of both salts

%% Setup paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_mixtures'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory
fig_out_dir = fullfile(filepath, '..', 'figures', 'multi_salt');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% Constants
MWw = 18.015;     % Molecular weight of water (g/mol)
MW_NaCl = 58.443; % Molecular weight of NaCl (g/mol)
MW_LiCl = 42.394; % Molecular weight of LiCl (g/mol)

%% Define molality ranges
% Create grid of molalities
m_NaCl_vec = linspace(0.01, 6.0, 50);  % NaCl molality (mol/kg)
m_LiCl_vec = linspace(0.01, 6.0, 50);  % LiCl molality (mol/kg)

[m_NaCl_grid, m_LiCl_grid] = meshgrid(m_NaCl_vec, m_LiCl_vec);

fprintf('Calculating water activity for NaCl + LiCl mixtures...\n');
fprintf('NaCl molality range: %.2f to %.2f mol/kg\n', min(m_NaCl_vec), max(m_NaCl_vec));
fprintf('LiCl molality range: %.2f to %.2f mol/kg\n', min(m_LiCl_vec), max(m_LiCl_vec));
fprintf('Grid size: %d x %d = %d points\n', length(m_NaCl_vec), length(m_LiCl_vec), ...
    length(m_NaCl_vec) * length(m_LiCl_vec));

%% Calculate mass fractions and water activity
aw_grid = zeros(size(m_NaCl_grid));

for i = 1:numel(m_NaCl_grid)
    m_NaCl = m_NaCl_grid(i);
    m_LiCl = m_LiCl_grid(i);
    
    % Convert molality to mass fraction
    % mass of water = 1000 g
    % mass of NaCl = m_NaCl * MW_NaCl
    % mass of LiCl = m_LiCl * MW_LiCl
    
    mass_water = 1000; % g
    mass_NaCl = m_NaCl * MW_NaCl; % g
    mass_LiCl = m_LiCl * MW_LiCl; % g
    total_mass = mass_water + mass_NaCl + mass_LiCl;
    
    mf_NaCl = mass_NaCl / total_mass;
    mf_LiCl = mass_LiCl / total_mass;
    
    % Calculate water activity using the fitted polynomial
    try
        aw = calculate_activity_NaCl_LiCl(mf_NaCl, mf_LiCl);
        
        % Ensure physically reasonable values
        if aw < 0 || aw > 1 || isnan(aw) || isinf(aw)
            aw = NaN;
        end
        
        aw_grid(i) = aw;
    catch
        aw_grid(i) = NaN;
    end
end

fprintf('Calculation complete!\n');
fprintf('Valid data points: %d / %d\n', sum(~isnan(aw_grid(:))), numel(aw_grid));

%% Plot 1: 3D Surface plot
figure('Position', [100, 100, 1200, 900]);
surf(m_NaCl_grid, m_LiCl_grid, aw_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

xlabel('NaCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of NaCl + LiCl Mixtures (25°C)', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

% Use a nice colormap
colormap(jet);
caxis([min(aw_grid(:)), max(aw_grid(:))]);

% Set view angle
view(-45, 30);

% Improve appearance
grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Set z-limits
zlim([0 1]);

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_3D_surface.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_3D_surface.fig'));

fprintf('Saved: NaCl_LiCl_water_activity_3D_surface.png\n');

%% Plot 2: Contour plot
figure('Position', [100, 100, 1200, 900]);
[C, h] = contourf(m_NaCl_grid, m_LiCl_grid, aw_grid, 20);
hold on;

xlabel('NaCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity Contours: NaCl + LiCl (25°C)', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity (a_w)';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(jet);

% Add contour labels
clabel(C, h, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_contour.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_contour.fig'));

fprintf('Saved: NaCl_LiCl_water_activity_contour.png\n');

%% Plot 3: 3D Surface with different viewing angle
figure('Position', [100, 100, 1200, 900]);
surf(m_NaCl_grid, m_LiCl_grid, aw_grid, 'EdgeColor', 'interp', 'FaceAlpha', 0.85);
hold on;

xlabel('NaCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of NaCl + LiCl Mixtures (Top View)', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(parula);
caxis([min(aw_grid(:)), max(aw_grid(:))]);

% Top-down view
view(0, 90);

% Improve appearance
grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_topview.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_topview.fig'));

fprintf('Saved: NaCl_LiCl_water_activity_topview.png\n');

%% Plot 4: Multiple cross-sections at fixed NaCl molalities
figure('Position', [100, 100, 1200, 900]);
hold on; grid on; box on;

% Select several NaCl molality values
m_NaCl_fixed = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0];
n_NaCl = length(m_NaCl_fixed);

% Sequential colormap: light orange (dilute) -> dark red (concentrated) for NaCl
colors_NaCl = interp1([1; n_NaCl], [0.98 0.78 0.50; 0.60 0.05 0.05], (1:n_NaCl)');

% Distinct markers for each line
markers = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_NaCl
    % Find closest index
    [~, idx] = min(abs(m_NaCl_vec - m_NaCl_fixed(i)));
    
    n_pts = length(m_LiCl_vec);
    mk_idx = round(linspace(1, n_pts, 8));
    
    % Plot water activity vs LiCl molality at this NaCl molality
    plot(m_LiCl_vec, aw_grid(:, idx), 'LineWidth', 2.5, ...
        'Color', colors_NaCl(i, :), ...
        'DisplayName', sprintf('NaCl = %.1f mol/kg', m_NaCl_vec(idx)));
    plot(m_LiCl_vec(mk_idx), aw_grid(mk_idx, idx), ...
        markers{i}, 'MarkerSize', 7, ...
        'Color', colors_NaCl(i, :), ...
        'MarkerFaceColor', colors_NaCl(i, :), ...
        'HandleVisibility', 'off');
end

% Add colorbar to indicate NaCl concentration gradient
colormap(gca, interp1([1; n_NaCl], [0.98 0.78 0.50; 0.60 0.05 0.05], linspace(1,n_NaCl,256)'));
cb = colorbar;
cb.Ticks = linspace(0,1,n_NaCl);
cb.TickLabels = arrayfun(@(x) sprintf('%.1f mol/kg', x), m_NaCl_fixed, 'UniformOutput', false);
cb.Label.String = 'NaCl Molality';
cb.Label.FontSize = 13;
cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity vs LiCl Molality at Fixed NaCl Concentrations', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_crosssections_fixed_NaCl.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_crosssections_fixed_NaCl.fig'));

fprintf('Saved: NaCl_LiCl_crosssections_fixed_NaCl.png\n');

%% Plot 5: Multiple cross-sections at fixed LiCl molalities
figure('Position', [100, 100, 1200, 900]);
hold on; grid on; box on;

% Select several LiCl molality values
m_LiCl_fixed = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0];
n_LiCl = length(m_LiCl_fixed);

% Sequential colormap: light blue (dilute) -> dark navy (concentrated) for LiCl
colors_LiCl = interp1([1; n_LiCl], [0.50 0.80 0.97; 0.03 0.15 0.55], (1:n_LiCl)');

% Distinct markers for each line
markers = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_LiCl
    % Find closest index
    [~, idx] = min(abs(m_LiCl_vec - m_LiCl_fixed(i)));
    
    n_pts = length(m_NaCl_vec);
    mk_idx = round(linspace(1, n_pts, 8));
    
    % Plot water activity vs NaCl molality at this LiCl molality
    plot(m_NaCl_vec, aw_grid(idx, :), 'LineWidth', 2.5, ...
        'Color', colors_LiCl(i, :), ...
        'DisplayName', sprintf('LiCl = %.1f mol/kg', m_LiCl_vec(idx)));
    plot(m_NaCl_vec(mk_idx), aw_grid(idx, mk_idx), ...
        markers{i}, 'MarkerSize', 7, ...
        'Color', colors_LiCl(i, :), ...
        'MarkerFaceColor', colors_LiCl(i, :), ...
        'HandleVisibility', 'off');
end

% Add colorbar to indicate LiCl concentration gradient
colormap(gca, interp1([1; n_LiCl], [0.50 0.80 0.97; 0.03 0.15 0.55], linspace(1,n_LiCl,256)'));
cb = colorbar;
cb.Ticks = linspace(0,1,n_LiCl);
cb.TickLabels = arrayfun(@(x) sprintf('%.1f mol/kg', x), m_LiCl_fixed, 'UniformOutput', false);
cb.Label.String = 'LiCl Molality';
cb.Label.FontSize = 13;
cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('NaCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity vs NaCl Molality at Fixed LiCl Concentrations', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_crosssections_fixed_LiCl.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_crosssections_fixed_LiCl.fig'));

fprintf('Saved: NaCl_LiCl_crosssections_fixed_LiCl.png\n');

%% Plot 6: Interactive 3D surface with mesh
figure('Position', [100, 100, 1400, 900]);

% Create surface
surf(m_NaCl_grid, m_LiCl_grid, aw_grid, 'FaceAlpha', 0.8);
hold on;

% Add mesh overlay
mesh(m_NaCl_grid, m_LiCl_grid, aw_grid, 'EdgeColor', 'k', 'FaceColor', 'none', 'EdgeAlpha', 0.15);

xlabel('NaCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of NaCl + LiCl Mixtures (Interactive View)', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(turbo);
caxis([min(aw_grid(:)), max(aw_grid(:))]);

% Set viewing angle
view(-135, 25);

% Add lighting for better 3D perception
lighting gouraud;
camlight('headlight');
material shiny;

% Improve appearance
grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');
zlim([0 1]);

% Save
saveas(gcf, fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_3D_mesh.png'));
savefig(fullfile(fig_out_dir, 'NaCl_LiCl_water_activity_3D_mesh.fig'));

fprintf('Saved: NaCl_LiCl_water_activity_3D_mesh.png\n');

%% Summary statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Water activity range: %.4f to %.4f\n', min(aw_grid(:)), max(aw_grid(:)));
fprintf('Mean water activity: %.4f\n', mean(aw_grid(:), 'omitnan'));
fprintf('Std water activity: %.4f\n', std(aw_grid(:), 'omitnan'));

% Find minimum water activity point
[min_aw, min_idx] = min(aw_grid(:));
[min_i, min_j] = ind2sub(size(aw_grid), min_idx);
fprintf('\nLowest water activity: %.4f\n', min_aw);
fprintf('  at NaCl = %.2f mol/kg, LiCl = %.2f mol/kg\n', ...
    m_NaCl_grid(min_i, min_j), m_LiCl_grid(min_i, min_j));

% Find maximum water activity point (excluding edges)
[max_aw, max_idx] = max(aw_grid(:));
[max_i, max_j] = ind2sub(size(aw_grid), max_idx);
fprintf('\nHighest water activity: %.4f\n', max_aw);
fprintf('  at NaCl = %.2f mol/kg, LiCl = %.2f mol/kg\n', ...
    m_NaCl_grid(max_i, max_j), m_LiCl_grid(max_i, max_j));

fprintf('\n=== All figures saved to: ===\n%s\n', fig_out_dir);

fprintf('\nVisualization complete!\n');
