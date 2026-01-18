close all
clear
clc

% Add calculate_mf folder to path
addpath('calculate_mf');

% This script plots water activity coefficient versus mean charge density
% Charge density is calculated as z/r^3 where z is charge and r is ionic radius
% Mean charge density is the average of cation and anion charge densities

T = 25; 
MWw = 18;

% Define RH value at which to evaluate the activity coefficient (can be adjusted)
RH_eval = 0.70;  % 70% relative humidity

%% Define salt properties
% Format: [Name, MW, Cation, Cation_Charge, Cation_Radius_pm, Cation_MW, Anion_Charge, Anion_Radius_pm, Anion_MW, calculate_function, RH_range_min, RH_range_max]

% Ionic radii (in picometers) from Shannon (1976) for 6-coordinate ions
% Cations: Li+: 76 pm, Na+: 102 pm, Ca2+: 100 pm, Mg2+: 72 pm, Zn2+: 74 pm, H+: ~10 pm (bare proton, unrealistic)
% For H+, using hydronium H3O+ effective radius ~140 pm
% Anions: Cl-: 181 pm, Br-: 196 pm, I-: 220 pm, OH-: ~140 pm, NO3-: ~200 pm (approximate for polyatomic)
% Molar masses (g/mol): Li+: 6.94, Na+: 22.99, Ca2+: 40.08, Mg2+: 24.31, Zn2+: 65.38, H+: 1.008
% Anions: Cl-: 35.45, Br-: 79.90, I-: 126.90, OH-: 17.01, NO3-: 62.00

salts = {
    % Salt name, MW, Cation, Cat_Charge, Cat_Radius(pm), Cat_MW, An_Charge, An_Radius(pm), An_MW, Calc_func, RH_min, RH_max, Color
    'LiCl',    42.4,    'Li^+',     1,  76,  6.94,   1,  181,  35.45,  @calculate_mf_LiCl,   0.12, 0.9, [0, 0.5, 0];
    % 'LiOH',    24,      'Li^+',     1,  76,  6.94,   1,  140,  17.01,  @calculate_mf_LiOH,   0.85, 0.9, [1, 0, 1];
    'LiBr',    86.85,   'Li^+',     1,  76,  6.94,   1,  196,  79.90,  @calculate_mf_LiBr,   0.07, 0.9, [0.3010, 0.7450, 0.9330];
    'LiI',     133.85,  'Li^+',     1,  76,  6.94,   1,  220,  126.90, @calculate_mf_LiI,    0.18, 0.9, [0.4940, 0.1840, 0.5560];
    'NaOH',    40,      'Na^+',     1,  102, 22.99,  1,  140,  17.01,  @calculate_mf_NaOH,   0.23, 0.9, [0.75, 0, 0.75];
    'CaCl_2',  111,     'Ca^{2+}',  2,  100, 40.08,  1,  181,  35.45,  @calculate_mf_CaCl,   0.31, 0.9, [0.9290, 0.6940, 0.1250];
    'MgCl_2',  95.2,    'Mg^{2+}',  2,  72,  24.31,  1,  181,  35.45,  @calculate_mf_MgCl,   0.33, 0.9, [0.8500, 0.3250, 0.0980];
    'Mg(NO_3)_2', 148.3, 'Mg^{2+}', 2,  72,  24.31,  1,  200,  62.00,  @calculate_mf_MgNO3,  0.55, 0.9, [0, 1, 1];
    'ZnCl_2',  136.3,   'Zn^{2+}',  2,  74,  65.38,  1,  181,  35.45,  @calculate_mf_ZnCl,   0.07, 0.8, [0.6350, 0.0780, 0.1840];
    'ZnBr_2',  225.2,   'Zn^{2+}',  2,  74,  65.38,  1,  196,  79.90,  @calculate_mf_ZnBr,   0.08, 0.85, [0.4660, 0.6740, 0.1880];
    'ZnI_2',   319.18,  'Zn^{2+}',  2,  74,  65.38,  1,  220,  126.90, @calculate_mf_ZnI,    0.25, 0.9, [0.9290, 0.6940, 0.1250];
    'HCl',     36.5,    'H^+',      1,  140, 1.008,  1,  181,  35.45,  @calculate_mf_HCl,    0.17, 0.9, [0, 0.4470, 0.7410];
};

n_salts = size(salts, 1);

%% Calculate charge density, charge per molar mass, and activity coefficients for each salt
charge_density = zeros(n_salts, 1);
charge_per_MW = zeros(n_salts, 1);
gamma_w = zeros(n_salts, 1);
cation_names = cell(n_salts, 1);
salt_names = cell(n_salts, 1);
colors = zeros(n_salts, 3);

for i = 1:n_salts
    salt_name = salts{i, 1};
    MW_salt = salts{i, 2};
    cation = salts{i, 3};
    cat_charge = salts{i, 4};
    cat_radius_pm = salts{i, 5};
    cat_MW = salts{i, 6};
    an_charge = salts{i, 7};
    an_radius_pm = salts{i, 8};
    an_MW = salts{i, 9};
    calc_func = salts{i, 10};
    RH_min = salts{i, 11};
    RH_max = salts{i, 12};
    color = salts{i, 13};
    
    % Calculate charge density for cation and anion (charge per unit volume)
    % z / r^3, where r is in pm, result in units of e/pm^3
    cat_charge_density = cat_charge / (cat_radius_pm^3);
    an_charge_density = an_charge / (an_radius_pm^3);
    
    % Calculate mean charge density
    charge_density(i) = (cat_charge_density + an_charge_density) / 2;
    
    % Calculate charge per molar mass for cation and anion (e/(g/mol))
    cat_charge_per_MW = cat_charge / cat_MW;
    an_charge_per_MW = an_charge / an_MW;
    
    % Calculate mean charge per molar mass
    charge_per_MW(i) = (cat_charge_per_MW + an_charge_per_MW) / 2;
    
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
xlabel('Mean Charge Density (e/pm^3) = (z_{cat}/r_{cat}^3 + z_{an}/r_{an}^3)/2', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('Water Activity Coefficient vs Mean Charge Density (RH approx %.0f%%)', RH_eval*100), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 10)
xlim([0, max(charge_density)*1.1])
ylim([0.7, 0.8])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save figure
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print('figures/Activity_Coefficient_vs_Charge_Density', '-dpng', '-r600')
savefig('figures/Activity_Coefficient_vs_Charge_Density.fig')

%% Create second plot: Water Activity Coefficient vs Mean Charge per Molar Mass
figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% Plot each salt with its specific color and label
for i = 1:n_salts
    scatter(charge_per_MW(i), gamma_w(i), 150, colors(i, :), 'filled', ...
            'DisplayName', sprintf('%s (%s)', salt_names{i}, cation_names{i}), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% Add ideal line (gamma_w = 1)
plot([0, max(charge_per_MW)*1.1], [1, 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');

% Labels and formatting
xlabel('Mean Charge per Molar Mass (e/(g/mol)) = (z_{cat}/MW_{cat} + z_{an}/MW_{an})/2', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Water Activity Coefficient (\gamma_w)', 'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('Water Activity Coefficient vs Mean Charge per Molar Mass (RH approx %.0f%%)', RH_eval*100), ...
      'FontSize', 16, 'FontWeight', 'bold')
legend('Location', 'best', 'FontSize', 10)
xlim([0, max(charge_per_MW)*1.1])
ylim([0.7, 0.8])
set(gca, 'FontSize', 12)
set(gcf, 'color', 'w');

% Save second figure
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 7]; 
print('figures/Activity_Coefficient_vs_Charge_per_MW', '-dpng', '-r600')
savefig('figures/Activity_Coefficient_vs_Charge_per_MW.fig')

%% Display results table
fprintf('\n=== Water Activity Coefficient vs Mean Charge Density ===\n');
fprintf('%-15s %-10s %-8s %-12s %-8s %-12s %-15s\n', 'Salt', 'Cation', 'Cat_z', 'Cat_r (pm)', 'An_z', 'An_r (pm)', 'Mean z/r^3');
fprintf('%-15s %-10s %-8s %-12s %-8s %-12s %-15s\n', repmat('-', 1, 15), repmat('-', 1, 10), ...
        repmat('-', 1, 8), repmat('-', 1, 12), repmat('-', 1, 8), repmat('-', 1, 12), repmat('-', 1, 15));

% Sort by charge density for display
[charge_density_sorted, sort_idx] = sort(charge_density);
for i = 1:n_salts
    idx = sort_idx(i);
    fprintf('%-15s %-10s %-8d %-12.1f %-8d %-12.1f %-15.6e\n', ...
            salt_names{idx}, strrep(cation_names{idx}, '^', ''), ...
            salts{idx, 4}, salts{idx, 5}, salts{idx, 7}, salts{idx, 8}, charge_density_sorted(i));
end

fprintf('\n=== Water Activity Coefficient vs Mean Charge per Molar Mass ===\n');
fprintf('%-15s %-10s %-8s %-10s %-8s %-10s %-15s\n', 'Salt', 'Cation', 'Cat_z', 'Cat_MW', 'An_z', 'An_MW', 'Mean z/MW');
fprintf('%-15s %-10s %-8s %-10s %-8s %-10s %-15s\n', repmat('-', 1, 15), repmat('-', 1, 10), ...
        repmat('-', 1, 8), repmat('-', 1, 10), repmat('-', 1, 8), repmat('-', 1, 10), repmat('-', 1, 15));

% Sort by charge per MW for display
[charge_per_MW_sorted, sort_idx2] = sort(charge_per_MW);
for i = 1:n_salts
    idx = sort_idx2(i);
    fprintf('%-15s %-10s %-8d %-10.2f %-8d %-10.2f %-15.6e\n', ...
            salt_names{idx}, strrep(cation_names{idx}, '^', ''), ...
            salts{idx, 4}, salts{idx, 6}, salts{idx, 7}, salts{idx, 9}, charge_per_MW_sorted(i));
end

fprintf('\n%-15s %-20s\n', 'Salt', 'gamma_w');
fprintf('%-15s %-20s\n', repmat('-', 1, 15), repmat('-', 1, 20));
for i = 1:n_salts
    idx = sort_idx(i);
    fprintf('%-15s %-20.4f\n', salt_names{idx}, gamma_w(idx));
end

fprintf('\nPlots generated successfully!\n');
fprintf('  - figures/Activity_Coefficient_vs_Charge_Density.png\n');
fprintf('  - figures/Activity_Coefficient_vs_Charge_Density.fig\n');
fprintf('  - figures/Activity_Coefficient_vs_Charge_per_MW.png\n');
fprintf('  - figures/Activity_Coefficient_vs_Charge_per_MW.fig\n');