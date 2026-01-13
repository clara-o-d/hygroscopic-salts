close all
clear
clc

% This script plots water activity coefficient versus cation charge density
% Charge density is calculated as z/r^3 where z is charge and r is ionic radius

T = 25; 
MWw = 18;

% Define RH value at which to evaluate the activity coefficient (can be adjusted)
RH_eval = 0.7;  % 70% relative humidity

%% Define salt properties
% Format: [Name, MW, Cation, Charge, Ionic_Radius_pm, calculate_function, RH_range_min, RH_range_max]

% Ionic radii (in picometers) from Shannon (1976) for 6-coordinate ions
% Li+: 76 pm, Na+: 102 pm, Ca2+: 100 pm, Mg2+: 72 pm, Zn2+: 74 pm, H+: ~10 pm (bare proton, unrealistic)
% For H+, using hydronium H3O+ effective radius ~140 pm

salts = {
    % Salt name, MW, Cation, Charge, Radius(pm), Calc_func, RH_min, RH_max, Color
    'LiCl',    42.4,    'Li^+',     1,  76,  @calculate_mf_LiCl,   0.12, 0.9, [0, 0.5, 0];
    'LiOH',    24,      'Li^+',     1,  76,  @calculate_mf_LiOH,   0.85, 0.9, [1, 0, 1];
    'LiBr',    86.85,   'Li^+',     1,  76,  @calculate_mf_LiBr,   0.07, 0.9, [0.3010, 0.7450, 0.9330];
    'LiI',     133.85,  'Li^+',     1,  76,  @calculate_mf_LiI,    0.18, 0.9, [0.4940, 0.1840, 0.5560];
    'NaOH',    40,      'Na^+',     1,  102, @calculate_mf_NaOH,   0.23, 0.9, [0.75, 0, 0.75];
    'CaCl_2',  111,     'Ca^{2+}',  2,  100, @calculate_mf_CaCl,   0.31, 0.9, [0.9290, 0.6940, 0.1250];
    'MgCl_2',  95.2,    'Mg^{2+}',  2,  72,  @calculate_mf_MgCl,   0.33, 0.9, [0.8500, 0.3250, 0.0980];
    'Mg(NO_3)_2', 148.3, 'Mg^{2+}', 2,  72,  @calculate_mf_MgNO3,  0.55, 0.9, [0, 1, 1];
    'ZnCl_2',  136.3,   'Zn^{2+}',  2,  74,  @calculate_mf_ZnCl,   0.07, 0.8, [0.6350, 0.0780, 0.1840];
    'ZnBr_2',  225.2,   'Zn^{2+}',  2,  74,  @calculate_mf_ZnBr,   0.08, 0.85, [0.4660, 0.6740, 0.1880];
    'ZnI_2',   319.18,  'Zn^{2+}',  2,  74,  @calculate_mf_ZnI,    0.25, 0.9, [0.9290, 0.6940, 0.1250];
    'HCl',     36.5,    'H^+',      1,  140, @calculate_mf_HCl,    0.17, 0.9, [0, 0.4470, 0.7410];
};

n_salts = size(salts, 1);

%% Calculate charge density and activity coefficients for each salt
charge_density = zeros(n_salts, 1);
gamma_w = zeros(n_salts, 1);
cation_names = cell(n_salts, 1);
salt_names = cell(n_salts, 1);
colors = zeros(n_salts, 3);

for i = 1:n_salts
    salt_name = salts{i, 1};
    MW_salt = salts{i, 2};
    cation = salts{i, 3};
    charge = salts{i, 4};
    radius_pm = salts{i, 5};
    calc_func = salts{i, 6};
    RH_min = salts{i, 7};
    RH_max = salts{i, 8};
    color = salts{i, 9};
    
    % Calculate charge density (charge per unit volume)
    % z / r^3, where r is in pm, result in units of e/pm^3
    charge_density(i) = charge / (radius_pm^3);
    
    % Check if RH_eval is within valid range for this salt
    if RH_eval < RH_min || RH_eval > RH_max
        % Use midpoint of valid range instead
        RH_use = (RH_min + RH_max) / 2;
        fprintf('Warning: RH_eval = %.2f is outside valid range [%.2f, %.2f] for %s. Using RH = %.2f instead.\n', ...
                RH_eval, RH_min, RH_max, salt_name, RH_use);
    else
        RH_use = RH_eval;
    end
    
    % Calculate mass fraction and mole fraction of water
    if strcmp(salt_name, 'LiCl') || strcmp(salt_name, 'CaCl_2')
        mf_salt = calc_func(RH_use, T);
    else
        mf_salt = calc_func(RH_use);
    end
    
    mf_water = 1 - mf_salt;
    x_water = (mf_water / MWw) / ((mf_water / MWw) + (mf_salt / MW_salt));
    
    % Calculate water activity coefficient: gamma_w = a_w / x_w = RH / x_w
    gamma_w(i) = RH_use / x_water;
    
    cation_names{i} = cation;
    salt_names{i} = salt_name;
    colors(i, :) = color;
end

%% Create the plot
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% Plot each salt with its specific color and label
for i = 1:n_salts
    scatter(charge_density(i), gamma_w(i), 150, colors(i, :), 'filled', ...
            'DisplayName', sprintf('%s (%s)', salt_names{i}, cation_names{i}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% Add ideal line (gamma_w = 1)
plot([0, max(charge_density)*1.1], [1, 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');

% Labels and formatting
xlabel('Cation Charge Density (e/pm^3) = z/r^3', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('Water Activity Coefficient vs Cation Charge Density (RH \approx %.0f%%)', RH_eval*100), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 10)
xlim([0, max(charge_density)*1.1])
ylim([0, max(gamma_w)*1.15])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save figure
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print('Activity_Coefficient_vs_Charge_Density', '-dpng', '-r600')
savefig('Activity_Coefficient_vs_Charge_Density.fig')

%% Display results table
fprintf('\n=== Water Activity Coefficient vs Cation Charge Density ===\n');
fprintf('%-15s %-10s %-10s %-15s %-15s\n', 'Salt', 'Cation', 'Charge', 'Radius (pm)', 'z/r^3 (e/pm^3)');
fprintf('%-15s %-10s %-10s %-15s %-15s\n', repmat('-', 1, 15), repmat('-', 1, 10), ...
        repmat('-', 1, 10), repmat('-', 1, 15), repmat('-', 1, 15));

% Sort by charge density for display
[charge_density_sorted, sort_idx] = sort(charge_density);
for i = 1:n_salts
    idx = sort_idx(i);
    fprintf('%-15s %-10s %-10d %-15.1f %-15.6e\n', ...
            salt_names{idx}, strrep(cation_names{idx}, '^', ''), ...
            salts{idx, 4}, salts{idx, 5}, charge_density_sorted(i));
end

fprintf('\n%-15s %-20s\n', 'Salt', 'gamma_w');
fprintf('%-15s %-20s\n', repmat('-', 1, 15), repmat('-', 1, 20));
for i = 1:n_salts
    idx = sort_idx(i);
    fprintf('%-15s %-20.4f\n', salt_names{idx}, gamma_w(idx));
end

fprintf('\nPlot generated successfully!\n');
fprintf('  - Activity_Coefficient_vs_Charge_Density.png\n');
fprintf('  - Activity_Coefficient_vs_Charge_Density.fig\n');
