close all
clear
clc

% Script to visualize water activity for LiCl as a function of molality and temperature
% Creates 3D surface plots and cross-sections showing temperature dependence

%% Setup paths
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_activity_temperature'));
addpath(fullfile(filepath, '..', 'util'));

% Define output directory
fig_out_dir = fullfile(filepath, '..', 'figures', 'temperature');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

%% Constants
MWw = 18.015;     % Molecular weight of water (g/mol)
MW_LiCl = 42.394; % Molecular weight of LiCl (g/mol)

%% Define ranges
% Molality range (mol/kg H2O)
m_LiCl_vec = linspace(0.5, 15.0, 50);

% Temperature range (°C)
T_vec = linspace(25, 100, 50);

[m_grid, T_grid] = meshgrid(m_LiCl_vec, T_vec);

fprintf('Calculating water activity for LiCl at different temperatures and molalities...\n');
fprintf('Molality range: %.2f to %.2f mol/kg\n', min(m_LiCl_vec), max(m_LiCl_vec));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min(T_vec), max(T_vec));
fprintf('Grid size: %d x %d = %d points\n', length(m_LiCl_vec), length(T_vec), ...
    length(m_LiCl_vec) * length(T_vec));

%% Calculate mass fractions and water activity
aw_grid = zeros(size(m_grid));

for i = 1:numel(m_grid)
    m_LiCl = m_grid(i);
    T = T_grid(i);
    
    % Convert molality to mass fraction
    % mass of water = 1000 g
    % mass of LiCl = m_LiCl * MW_LiCl
    
    mass_water = 1000; % g
    mass_LiCl = m_LiCl * MW_LiCl; % g
    total_mass = mass_water + mass_LiCl;
    
    mf_LiCl = mass_LiCl / total_mass;
    
    % Calculate water activity using the temperature-dependent fitted polynomial
    try
        aw = calculate_activity_temperature_LiCl(mf_LiCl, T);
        
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
surf(m_grid, T_grid, aw_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of LiCl Solutions vs Temperature and Molality', 'FontSize', 20, 'FontWeight', 'bold');

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
saveas(gcf, fullfile(fig_out_dir, 'LiCl_water_activity_3D_temp_molality.png'));
savefig(fullfile(fig_out_dir, 'LiCl_water_activity_3D_temp_molality.fig'));

fprintf('Saved: LiCl_water_activity_3D_temp_molality.png\n');

%% Plot 2: Contour plot
figure('Position', [100, 100, 1200, 900]);
[C, h] = contourf(m_grid, T_grid, aw_grid, 20);
hold on;

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity Contours: LiCl Solutions', 'FontSize', 20, 'FontWeight', 'bold');

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
saveas(gcf, fullfile(fig_out_dir, 'LiCl_water_activity_contour_temp_molality.png'));
savefig(fullfile(fig_out_dir, 'LiCl_water_activity_contour_temp_molality.fig'));

fprintf('Saved: LiCl_water_activity_contour_temp_molality.png\n');

%% Plot 3: 3D Surface with different viewing angle (side view)
figure('Position', [100, 100, 1200, 900]);
surf(m_grid, T_grid, aw_grid, 'EdgeColor', 'interp', 'FaceAlpha', 0.85);
hold on;

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of LiCl Solutions (Side View)', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(parula);
caxis([min(aw_grid(:)), max(aw_grid(:))]);

% Side view to see temperature dependence
view(90, 0);

% Improve appearance
grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');
zlim([0 1]);

% Save
saveas(gcf, fullfile(fig_out_dir, 'LiCl_water_activity_sideview_temp.png'));
savefig(fullfile(fig_out_dir, 'LiCl_water_activity_sideview_temp.fig'));

fprintf('Saved: LiCl_water_activity_sideview_temp.png\n');

%% Plot 4: Cross-sections at fixed temperatures
figure('Position', [100, 100, 1200, 900]);
hold on; grid on; box on;

% Select several temperatures
T_fixed = [25, 40, 55, 70, 85, 100];
n_T = length(T_fixed);

% Cool-to-warm colormap: blue (cold) -> red (hot)
colors_T = interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], (1:n_T)');

% Distinct markers for each line
markers = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_T
    % Find closest index
    [~, idx] = min(abs(T_vec - T_fixed(i)));
    
    % Plot water activity vs molality at this temperature
    % Use marker spacing so markers don't crowd the line
    n_pts = length(m_LiCl_vec);
    mk_idx = round(linspace(1, n_pts, 8));
    
    plot(m_LiCl_vec, aw_grid(idx, :), 'LineWidth', 2.5, ...
        'Color', colors_T(i, :), ...
        'DisplayName', sprintf('T = %.0f°C', T_vec(idx)));
    plot(m_LiCl_vec(mk_idx), aw_grid(idx, mk_idx), ...
        markers{i}, 'MarkerSize', 7, ...
        'Color', colors_T(i, :), ...
        'MarkerFaceColor', colors_T(i, :), ...
        'HandleVisibility', 'off');
end

% Add colorbar to indicate temperature gradient
colormap(gca, interp1([1; n_T], [0.08 0.40 0.85; 0.88 0.10 0.10], linspace(1,n_T,256)'));
cb = colorbar;
cb.Ticks = linspace(0,1,n_T);
cb.TickLabels = arrayfun(@(x) sprintf('%.0f°C', x), T_fixed, 'UniformOutput', false);
cb.Label.String = 'Temperature';
cb.Label.FontSize = 13;
cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity vs Molality at Fixed Temperatures', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'LiCl_crosssections_fixed_temp.png'));
savefig(fullfile(fig_out_dir, 'LiCl_crosssections_fixed_temp.fig'));

fprintf('Saved: LiCl_crosssections_fixed_temp.png\n');

%% Plot 5: Cross-sections at fixed molalities
figure('Position', [100, 100, 1200, 900]);
hold on; grid on; box on;

% Select several molalities
m_fixed = [1.0, 3.0, 5.0, 7.0, 9.0, 12.0];
n_m = length(m_fixed);

% Sequential colormap: light blue (dilute) -> dark navy (concentrated)
colors_m = interp1([1; n_m], [0.40 0.76 0.96; 0.05 0.10 0.45], (1:n_m)');

% Distinct markers for each line
markers = {'o', 's', '^', 'd', 'v', 'p'};

for i = 1:n_m
    % Find closest index
    [~, idx] = min(abs(m_LiCl_vec - m_fixed(i)));
    
    n_pts = length(T_vec);
    mk_idx = round(linspace(1, n_pts, 8));
    
    plot(T_vec, aw_grid(:, idx), 'LineWidth', 2.5, ...
        'Color', colors_m(i, :), ...
        'DisplayName', sprintf('m = %.1f mol/kg', m_LiCl_vec(idx)));
    plot(T_vec(mk_idx), aw_grid(mk_idx, idx), ...
        markers{i}, 'MarkerSize', 7, ...
        'Color', colors_m(i, :), ...
        'MarkerFaceColor', colors_m(i, :), ...
        'HandleVisibility', 'off');
end

% Add colorbar to indicate concentration gradient
colormap(gca, interp1([1; n_m], [0.40 0.76 0.96; 0.05 0.10 0.45], linspace(1,n_m,256)'));
cb = colorbar;
cb.Ticks = linspace(0,1,n_m);
cb.TickLabels = arrayfun(@(x) sprintf('%.1f mol/kg', x), m_fixed, 'UniformOutput', false);
cb.Label.String = 'LiCl Molality';
cb.Label.FontSize = 13;
cb.Label.FontWeight = 'bold';
clim([0 1]);

xlabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity vs Temperature at Fixed Molalities', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'LiCl_crosssections_fixed_molality.png'));
savefig(fullfile(fig_out_dir, 'LiCl_crosssections_fixed_molality.fig'));

fprintf('Saved: LiCl_crosssections_fixed_molality.png\n');

%% Plot 6: Interactive 3D surface with mesh and lighting
figure('Position', [100, 100, 1400, 900]);

% Create surface
surf(m_grid, T_grid, aw_grid, 'FaceAlpha', 0.8);
hold on;

% Add mesh overlay
mesh(m_grid, T_grid, aw_grid, 'EdgeColor', 'k', 'FaceColor', 'none', 'EdgeAlpha', 0.15);

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Water Activity (a_w)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity of LiCl Solutions (Interactive View)', 'FontSize', 20, 'FontWeight', 'bold');

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
saveas(gcf, fullfile(fig_out_dir, 'LiCl_water_activity_3D_mesh_temp.png'));
savefig(fullfile(fig_out_dir, 'LiCl_water_activity_3D_mesh_temp.fig'));

fprintf('Saved: LiCl_water_activity_3D_mesh_temp.png\n');

%% Plot 7: Heatmap with better visibility
figure('Position', [100, 100, 1200, 900]);
imagesc(m_LiCl_vec, T_vec, aw_grid);
set(gca, 'YDir', 'normal');

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
title('Water Activity Heatmap: LiCl Solutions', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'Water Activity (a_w)';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(jet);
caxis([min(aw_grid(:)), max(aw_grid(:))]);

% Add contour lines
hold on;
[C, h] = contour(m_grid, T_grid, aw_grid, 10, 'k', 'LineWidth', 0.5);
clabel(C, h, 'FontSize', 8, 'Color', 'k');

grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');

% Save
saveas(gcf, fullfile(fig_out_dir, 'LiCl_water_activity_heatmap_temp.png'));
savefig(fullfile(fig_out_dir, 'LiCl_water_activity_heatmap_temp.fig'));

fprintf('Saved: LiCl_water_activity_heatmap_temp.png\n');

%% Plot 8: Relative Humidity equivalent
% Convert water activity to RH percentage
RH_grid = aw_grid * 100;

figure('Position', [100, 100, 1200, 900]);
surf(m_grid, T_grid, RH_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;

xlabel('LiCl Molality (mol/kg H_2O)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 16, 'FontWeight', 'bold');
zlabel('Relative Humidity (%)', 'FontSize', 16, 'FontWeight', 'bold');
title('Equilibrium RH of LiCl Solutions vs Temperature and Molality', 'FontSize', 20, 'FontWeight', 'bold');

% Add colorbar
c = colorbar;
c.Label.String = 'RH (%)';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

colormap(jet);
caxis([min(RH_grid(:)), max(RH_grid(:))]);

view(-45, 30);

grid on;
box on;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');
zlim([0 100]);

% Save
saveas(gcf, fullfile(fig_out_dir, 'LiCl_RH_3D_temp_molality.png'));
savefig(fullfile(fig_out_dir, 'LiCl_RH_3D_temp_molality.fig'));

fprintf('Saved: LiCl_RH_3D_temp_molality.png\n');

%% Summary statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('Water activity range: %.4f to %.4f\n', min(aw_grid(:)), max(aw_grid(:)));
fprintf('Mean water activity: %.4f\n', mean(aw_grid(:), 'omitnan'));
fprintf('Std water activity: %.4f\n', std(aw_grid(:), 'omitnan'));

% Find minimum water activity point
[min_aw, min_idx] = min(aw_grid(:));
[min_i, min_j] = ind2sub(size(aw_grid), min_idx);
fprintf('\nLowest water activity: %.4f\n', min_aw);
fprintf('  at molality = %.2f mol/kg, T = %.1f°C\n', ...
    m_grid(min_i, min_j), T_grid(min_i, min_j));

% Find maximum water activity point
[max_aw, max_idx] = max(aw_grid(:));
[max_i, max_j] = ind2sub(size(aw_grid), max_idx);
fprintf('\nHighest water activity: %.4f\n', max_aw);
fprintf('  at molality = %.2f mol/kg, T = %.1f°C\n', ...
    m_grid(max_i, max_j), T_grid(max_i, max_j));

% Temperature dependence analysis
fprintf('\n=== Temperature Dependence ===\n');
fprintf('Effect of temperature on water activity:\n');
test_molalities = [2.0, 5.0, 10.0];
for m_test = test_molalities
    [~, idx] = min(abs(m_LiCl_vec - m_test));
    aw_25 = aw_grid(1, idx);  % First temperature (25°C)
    aw_100 = aw_grid(end, idx);  % Last temperature (100°C)
    delta_aw = aw_100 - aw_25;
    fprintf('  At m = %.1f mol/kg: Δaw = %.4f (from 25°C to 100°C)\n', ...
        m_LiCl_vec(idx), delta_aw);
end

fprintf('\n=== All figures saved to: ===\n%s\n', fig_out_dir);
fprintf('\nVisualization complete!\n');
