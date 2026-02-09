function plot_water_activity_all_salts()
    close all 
    clear
    clc 
    
    % --- 1. SETUP & PATHS ---
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'calculate_mf'));
    addpath(fullfile(filepath, '..', 'util'));
    
    fig_out_dir = fullfile(filepath, '..', 'figures', 'activity_coefficient');
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    fig_out_dir2 = fullfile(filepath, '..', 'figures', 'pitzer_comparisons');
    if ~exist(fig_out_dir2, 'dir')
        mkdir(fig_out_dir2);
    end
    
    T = 25; 
    MWw = 18.015; % g/mol for water
    MWw_kg = 0.018015; % kg/mol for Pitzer calculations
    
    % --- 2. SALT DEFINITIONS ---
    % {salt_name, MW, RH_min, RH_max, function_name, function_args, num_cations, num_anions}
    salt_data = {
        % Endothermic
        {'NaCl', 58.443, 0.765, 0.99, 'calculate_mf_NaCl', 0, 1, 1};
        {'KCl', 74.551, 0.855, 0.99, 'calculate_mf_KCl', 0, 1, 1};
        {'NH4Cl', 53.491, 0.815, 0.99, 'calculate_mf_NH4Cl', 0, 1, 1};
        {'CsCl', 168.363, 0.82, 0.99, 'calculate_mf_CsCl', 0, 1, 1};
        {'NaNO3', 85.00, 0.971, 0.995, 'calculate_mf_NaNO3', 0, 1, 1};
        {'AgNO3', 169.87, 0.865, 0.985, 'calculate_mf_AgNO3', 0, 1, 1};
        {'KI', 165.998, 0.97, 0.995, 'calculate_mf_KI', 0, 1, 1};
        {'LiNO3', 68.95, 0.736, 0.99, 'calculate_mf_LiNO3', 0, 1, 1};
        {'KNO3', 101.10, 0.932, 0.995, 'calculate_mf_KNO3', 0, 1, 1};
        {'NaClO4', 122.44, 0.778, 0.99, 'calculate_mf_NaClO4', 0, 1, 1};
        {'KClO3', 122.55, 0.981, 0.9926, 'calculate_mf_KClO3', 0, 1, 1};
        {'NaBr', 102.89, 0.614, 0.9280, 'calculate_mf_NaBr', 0, 1, 1};
        {'NaI', 149.89, 0.581, 0.9659, 'calculate_mf_NaI', 0, 1, 1};
        {'KBr', 119.00, 0.833, 0.9518, 'calculate_mf_KBr', 0, 1, 1};
        {'RbCl', 120.92, 0.743, 0.9517, 'calculate_mf_RbCl', 0, 1, 1};
        {'CsBr', 212.81, 0.848, 0.9472, 'calculate_mf_CsBr', 0, 1, 1};
        {'CsI', 259.81, 0.913, 0.9614, 'calculate_mf_CsI', 0, 1, 1};
        
        % Exothermic
        {'LiCl', 42.4, 0.12, 0.97, 'calculate_mf_LiCl', 1, 1, 1};
        {'LiOH', 24, 0.85, 0.97, 'calculate_mf_LiOH', 0, 1, 1};
        {'NaOH', 40, 0.23, 0.97, 'calculate_mf_NaOH', 0, 1, 1};
        {'HCl', 36.5, 0.17, 0.97, 'calculate_mf_HCl', 0, 1, 1};
        {'CaCl2', 111, 0.31, 0.97, 'calculate_mf_CaCl', 1, 1, 2};
        {'MgCl2', 95.2, 0.33, 0.97, 'calculate_mf_MgCl', 0, 1, 2};
        {'MgNO3', 148.3, 0.55, 0.9, 'calculate_mf_MgNO3', 0, 1, 2};
        {'LiBr', 86.85, 0.07, 0.97, 'calculate_mf_LiBr', 0, 1, 1};
        {'ZnCl2', 136.3, 0.07, 0.97, 'calculate_mf_ZnCl', 0, 1, 2};
        {'ZnI2', 319.18, 0.25, 0.97, 'calculate_mf_ZnI', 0, 1, 2};
        {'ZnBr2', 225.2, 0.08, 0.85, 'calculate_mf_ZnBr', 0, 1, 2};
        {'LiI', 133.85, 0.18, 0.97, 'calculate_mf_LiI', 0, 1, 1};
        
        % Sulfates
        {'Na2SO4', 142.04, 0.9000, 0.9947, 'calculate_mf_Na2SO4', 0, 2, 1};
        {'K2SO4', 174.26, 0.9730, 0.9948, 'calculate_mf_K2SO4', 0, 2, 1};
        {'NH42SO4', 132.14, 0.8320, 0.9949, 'calculate_mf_NH42SO4', 0, 2, 1};
        {'MgSO4', 120.37, 0.9060, 0.9950, 'calculate_mf_MgSO4', 0, 1, 1};
        {'MnSO4', 151.00, 0.9200, 0.9951, 'calculate_mf_MnSO4', 0, 1, 1};
        {'Li2SO4', 109.94, 0.8540, 0.9946, 'calculate_mf_Li2SO4', 0, 2, 1};
        {'NiSO4', 154.75, 0.9720, 0.9952, 'calculate_mf_NiSO4', 0, 1, 1};
        {'CuSO4', 159.61, 0.9760, 0.9953, 'calculate_mf_CuSO4', 0, 1, 1};
        {'ZnSO4', 161.44, 0.9390, 0.9952, 'calculate_mf_ZnSO4', 0, 1, 1};
        
        % Nitrates/Halides/Chlorates
        {'BaNO3', 261.34, 0.9869, 0.9948, 'calculate_mf_BaNO32', 0, 1, 2};
        {'CaNO3', 164.09, 0.6474, 0.9945, 'calculate_mf_CaNO32', 0, 1, 2};
        {'CaBr2', 199.89, 0.5114, 0.9318, 'calculate_mf_CaBr2', 0, 1, 2};
        {'CaI2', 293.89, 0.7590, 0.9294, 'calculate_mf_CaI2', 0, 1, 2};
        {'SrCl2', 158.53, 0.7235, 0.9669, 'calculate_mf_SrCl2', 0, 1, 2};
        {'SrBr2', 247.43, 0.6857, 0.9364, 'calculate_mf_SrBr2', 0, 1, 2};
        {'SrI2', 341.43, 0.5589, 0.9361, 'calculate_mf_SrI2', 0, 1, 2};
        {'BaCl2', 208.23, 0.9078, 0.9600, 'calculate_mf_BaCl2', 0, 1, 2};
        {'BaBr2', 297.14, 0.7454, 0.9387, 'calculate_mf_BaBr2', 0, 1, 2};
        {'NH4NO3', 80.043, 0.62, 0.95, 'calculate_mf_NH4NO3', 0, 1, 1};
        {'LiClO4', 106.39, 0.7785, 0.9869, 'calculate_mf_LiClO4', 0, 1, 1};
    };
    
    % --- 3. DATA PROCESSING LOOP ---
    num_points = 100;
    all_salts_data = struct();
    
    % Arrays to aggregate ALL data for the master fitting
    master_RH = [];         % 0-100 scale for fitting
    master_GammaIon = [];   % Ionic basis gamma
    master_Molality = [];   % Molality (m)
    master_Nu = [];         % Stoichiometric coefficient (nu)
    
    fprintf('Processing %d salts...\n', length(salt_data));
    
    for s = 1:length(salt_data)
        salt_name = salt_data{s}{1};
        MW = salt_data{s}{2};
        RH_min = salt_data{s}{3};
        RH_max = salt_data{s}{4};
        func_name = salt_data{s}{5};
        func_args = salt_data{s}{6};
        
        n_cat = salt_data{s}{7};
        n_an  = salt_data{s}{8};
        nu = n_cat + n_an; 
        
        RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
        mf_salt = zeros(size(RH_vec));
        mf_water = zeros(size(RH_vec));
        x_water_mol = zeros(size(RH_vec));
        x_water_ion = zeros(size(RH_vec));
        
        success = true;
        for i = 1:length(RH_vec)
            try
                if func_args == 1
                    mf_salt(i) = feval(func_name, RH_vec(i), T);
                else
                    mf_salt(i) = feval(func_name, RH_vec(i));
                end
                mf_water(i) = 1 - mf_salt(i);
                
                n_w = mf_water(i) / MWw;
                n_s = mf_salt(i) / MW;
                
                x_water_mol(i) = n_w / (n_w + n_s);
                x_water_ion(i) = n_w / (n_w + (nu * n_s));
                
            catch
                success = false;
                break;
            end
        end
        
        if success
            % Calc Gamma
            gamma_w_mol = RH_vec ./ x_water_mol;
            gamma_w_ion = RH_vec ./ x_water_ion;
            
            % Calc Molality (for Pitzer)
            % m = moles_salt / kg_water
            % m = (n_s) / (n_w * MWw_kg) 
            % From x_water_ion: m = (1/x_ion - 1) / (nu * MWw_kg)
            molality = (1./x_water_ion - 1) ./ (nu * MWw_kg);
            
            % Store Individual Data
            all_salts_data.(salt_name) = struct(...
                'RH', RH_vec, ...
                'x_water_mol', x_water_mol, 'gamma_w_mol', gamma_w_mol, ...
                'x_water_ion', x_water_ion, 'gamma_w_ion', gamma_w_ion, ...
                'molality', molality, ...
                'display_name', salt_name);
                
            % Aggregate Data for Master Fits
            % Filter out NaNs or Infs
            valid_mask = ~isnan(gamma_w_ion) & ~isinf(gamma_w_ion) & ~isnan(molality);
            
            master_RH = [master_RH, RH_vec(valid_mask) * 100]; % Scale 0-100
            master_GammaIon = [master_GammaIon, gamma_w_ion(valid_mask)];
            master_Molality = [master_Molality, molality(valid_mask)];
            master_Nu = [master_Nu, repmat(nu, 1, sum(valid_mask))];
        end
    end
    
    salt_names = fieldnames(all_salts_data);
    num_salts = length(salt_names);
    
    % --- Color Generation ---
    colors = generate_colors(num_salts);
    
    % --- 4. ORIGINAL PLOTTING (Code 1 Functionality) ---
    plot_original_curves(all_salts_data, salt_names, colors, fig_out_dir);
    
    % --- 5. NEW: 4 POLYNOMIAL FITS (Global) ---
    fprintf('Calculating Global Polynomial Fits...\n');
    perform_polynomial_fits(master_RH, master_GammaIon, fig_out_dir);
    
    % --- 6. NEW: PITZER FITTING (Global) ---
    fprintf('Calculating Global Pitzer Fit...\n');
    perform_pitzer_fit(master_RH, master_GammaIon, master_Molality, master_Nu, T, fig_out_dir2);
    
    disp('Done. Check figures folder.');
end

% =========================================================================
%                           SUB-FUNCTIONS
% =========================================================================

function perform_polynomial_fits(x_data, y_data, outdir)
    % x_data is RH (0-100), y_data is Gamma (Ionic)
    
    % Sort for plotting lines cleanly
    [x_sorted, sort_idx] = sort(x_data);
    
    % Definitions
    % Fit 1: Unconstrained Quadratic (y = p1*x^2 + p2*x + p3)
    p1 = polyfit(x_data, y_data, 2);
    y1 = polyval(p1, x_sorted);
    rmse1 = sqrt(mean((polyval(p1, x_data) - y_data).^2));
    eq1 = sprintf('y = %.2e x^2 + %.2e x + %.4f', p1(1), p1(2), p1(3));
    
    % Fit 2: Unconstrained Linear (y = p1*x + p2)
    p2 = polyfit(x_data, y_data, 1);
    y2 = polyval(p2, x_sorted);
    rmse2 = sqrt(mean((polyval(p2, x_data) - y_data).^2));
    eq2 = sprintf('y = %.4f x + %.4f', p2(1), p2(2));
    
    % Fit 3: Constrained Quadratic through (100, 1)
    % Model: y = a(x-100)^2 + b(x-100) + 1
    % Shift vars: Y = y - 1, X1 = (x-100)^2, X2 = (x-100)
    Y_shift = y_data(:) - 1;
    X_shift = x_data(:) - 100;
    X_mat = [X_shift.^2, X_shift];
    coeff3 = X_mat \ Y_shift; % [a; b]
    y3 = coeff3(1)*(x_sorted-100).^2 + coeff3(2)*(x_sorted-100) + 1;
    rmse3 = sqrt(mean(( (coeff3(1)*(x_data-100).^2 + coeff3(2)*(x_data-100) + 1) - y_data ).^2));
    eq3 = sprintf('y = %.2e(x-100)^2 + %.2e(x-100) + 1', coeff3(1), coeff3(2));
    
    % Fit 4: Constrained Linear through (100, 1)
    % Model: y = m(x-100) + 1
    X_lin = x_data(:) - 100;
    coeff4 = X_lin \ Y_shift; % [m]
    y4 = coeff4(1)*(x_sorted-100) + 1;
    rmse4 = sqrt(mean(( (coeff4(1)*(x_data-100) + 1) - y_data ).^2));
    eq4 = sprintf('y = %.4f(x-100) + 1', coeff4(1));

    % Fit 5: Constrained sqrt (fixed)
    y5 = 0.1*sqrt(x_sorted);
    rmse5 = sqrt(mean((y5 - y_data).^2));
    eq5 = sprintf('y = (0.1*sqrt(x))');

    % % Fit 6: Unconstrained sqrt with intercept
    % % Model: y = a*sqrt(x) + b
    % X_sqrt = [sqrt(x_data(:)), ones(size(x_data(:)))];
    % coeff6 = X_sqrt \ y_data(:);   % [a; b]
    % y6 = coeff6(1)*sqrt(x_sorted) + coeff6(2);
    % rmse6 = sqrt(mean((coeff6(1)*sqrt(x_data) + coeff6(2) - y_data).^2));
    % eq6 = sprintf('y = %.4f*sqrt(x) + %.4f', coeff6(1), coeff6(2));
        % Fit 6: Unconstrained sqrt + x*sqrt(x) with intercept
        % Model: y = a*sqrt(x) + b*x*sqrt(x) + c
        X_sqrt2 = [ ...
            sqrt(x_data(:)), ...
            x_data(:).*sqrt(x_data(:)), ...
            ones(size(x_data(:))) ...
        ];
        coeff6 = X_sqrt2 \ y_data(:);   % [a; b; c]
    
        y6 = coeff6(1)*sqrt(x_sorted) + ...
             coeff6(2)*x_sorted.*sqrt(x_sorted) + ...
             coeff6(3);
    
        rmse6 = sqrt(mean( ...
            (coeff6(1)*sqrt(x_data) + ...
             coeff6(2)*x_data.*sqrt(x_data) + ...
             coeff6(3) - y_data).^2 ));
    
        eq6 = sprintf('y = %.4f*sqrt(x) + %.4f*x*sqrt(x) + %.4f', ...
                      coeff6(1), coeff6(2), coeff6(3));


    
    % Plotting
    figure('Position', [100, 100, 1000, 800]);
    scatter(x_data, y_data, 10, [0.7 0.7 0.7], 'filled', 'DisplayName', 'Aggregated Data');
    hold on; grid on; box on;
    
    plot(x_sorted, y1, 'r--', 'LineWidth', 2, 'DisplayName', ['Fit 1 (Quad): ' eq1]);
    plot(x_sorted, y2, 'c--', 'LineWidth', 2, 'DisplayName', ['Fit 2 (Lin): ' eq2]);
    plot(x_sorted, y6, 'k--', 'LineWidth', 3, 'DisplayName', ['Fit 3 (Sqrt): ' eq6]);
    plot(x_sorted, y3, 'b-', 'LineWidth', 2, 'DisplayName', ['Fit 4 (Constr Quad): ' eq3]);
    plot(x_sorted, y4, 'm-', 'LineWidth', 2, 'DisplayName', ['Fit 5 (Constr Lin): ' eq4]);
    plot(x_sorted, y5, 'g-', 'LineWidth', 3, 'DisplayName', ['Fit 6 (Constr Sqrt): ' eq5]);


    
    plot(100, 1, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'DisplayName', 'Constraint (100, 1)');
    
    xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Activity Coeff \gamma_w (Ionic Basis)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Global Polynomial Fits to Aggregated Data', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'southwest', 'FontSize', 9);
    % ylim([0.8 1.4]);
    xlim([0 100]);
    
    % Add Text Box for RMSE
    stats_str = {
        'RMSE Values:';
        sprintf('1. Quad: %.4f', rmse1);
        sprintf('2. Lin: %.4f', rmse2);
        sprintf('3. Sqrt: %.4f', rmse6);
        sprintf('4. Constr Quad: %.4f', rmse3);
        sprintf('5. Constr Lin: %.4f', rmse4);
        sprintf('6. Constr Sqrt: %.4f', rmse5);

    };
    annotation('textbox', [0.15 0.6 0.2 0.2], 'String', stats_str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    
    print(fullfile(outdir, 'Activity_Coefficient_Ionic_vs_RH_Polynomial_Fits'), '-dpng', '-r300');
end

function perform_pitzer_fit(rh_vec, gamma_obs_vec, m_vec, nu_vec, T, outdir)
    % --- DEFINE SUBSETS ---
    idx_11 = (nu_vec == 2);      % 1:1 salts (e.g., NaCl)
    idx_non11 = (nu_vec ~= 2);   % Non 1:1 salts (e.g., MgCl2, Na2SO4)

    % --- HELPER FOR OPTIMIZATION (UNWEIGHTED) ---
    function [b0, b1, Cp, rmse, nu_mean] = run_fit(indices)
        rh_sub = rh_vec(indices);
        g_sub = gamma_obs_vec(indices);
        m_sub = m_vec(indices);
        nu_sub = nu_vec(indices);
        
        nu_mean = mean(nu_sub);
        
        % Use standard objective (no weights)
        cost_func = @(p) pitzer_objective(p, rh_sub, g_sub, m_sub, nu_sub, T);
        
        opts = optimset('Display', 'off', 'TolFun', 1e-7, 'MaxFunEvals', 3000, 'MaxIter', 1000);
        p_opt = fminsearch(cost_func, [0.1, 0.1, 0.001], opts);
        p_opt = fminsearch(cost_func, p_opt, opts); % Restart
        
        b0 = p_opt(1); b1 = p_opt(2); Cp = p_opt(3);
        
        % Standard RMSE
        g_fit = zeros(size(rh_sub));
        for k=1:length(rh_sub)
            [g_fit(k), ~] = calculate_pitzer_point(m_sub(k), nu_sub(k), T, b0, b1, Cp);
        end
        rmse = sqrt(mean((g_fit - g_sub).^2));
    end

    % --- RUN FITS ---
    fprintf('  > Fitting Global (Unweighted)...\n');
    [b0_all, b1_all, Cp_all, rmse_all, nu_all] = run_fit(1:length(rh_vec));
    
    fprintf('  > Fitting 1:1 Salts (Unweighted)...\n');
    [b0_11, b1_11, Cp_11, rmse_11, nu_11] = run_fit(idx_11);
    
    fprintf('  > Fitting Non-1:1 Salts (Unweighted)...\n');
    [b0_non, b1_non, Cp_non, rmse_non, nu_non] = run_fit(idx_non11);

    % --- PLOTTING ---
    figure('Position', [100, 50, 800, 1200], 'Color', 'w');
    t = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Helper to plot curve
    plot_model = @(b0, b1, cp, nu_val, col, style) ...
        plot_pitzer_curve(b0, b1, cp, nu_val, max(m_vec), T, col, style);

    % --- SUBPLOT 1: GLOBAL ---
    nexttile;
    hold on; grid on; box on;
    scatter(rh_vec, gamma_obs_vec, 15, [0.6 0.6 0.6], 'filled'); % All Grey
    plot_model(b0_all, b1_all, Cp_all, nu_all, 'r', '-');
    plot(100, 1, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    ylabel('\gamma_w (Ionic)');
    title(sprintf('Global Fit (All Salts, Mean \\nu=%.2f)', nu_all), 'FontWeight', 'bold');
    text(5, 1.05, sprintf('RMSE: %.4f\n\\beta_0: %.4f\n\\beta_1: %.4f\nC_{\\phi}: %.4f', ...
        rmse_all, b0_all, b1_all, Cp_all), 'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    xlim([0 100]);

    % --- SUBPLOT 2: 1:1 ONLY ---
    nexttile;
    hold on; grid on; box on;
    scatter(rh_vec(~idx_11), gamma_obs_vec(~idx_11), 10, [0.8 0.8 0.8], 'filled'); 
    scatter(rh_vec(idx_11), gamma_obs_vec(idx_11), 15, 'b', 'filled'); 
    plot_model(b0_11, b1_11, Cp_11, nu_11, 'b', '-');
    plot(100, 1, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    ylabel('\gamma_w (Ionic)');
    title('1:1 Salts Only (e.g., NaCl, KCl)', 'FontWeight', 'bold', 'Color', 'b');
    text(5, 1.05, sprintf('RMSE: %.4f\n\\beta_0: %.4f\n\\beta_1: %.4f\nC_{\\phi}: %.4f', ...
        rmse_11, b0_11, b1_11, Cp_11), 'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'b');
    xlim([0 100]);

    % --- SUBPLOT 3: NON-1:1 ONLY ---
    nexttile;
    hold on; grid on; box on;
    scatter(rh_vec(idx_11), gamma_obs_vec(idx_11), 10, [0.8 0.8 0.8], 'filled'); 
    scatter(rh_vec(idx_non11), gamma_obs_vec(idx_non11), 15, [0 0.5 0], 'filled'); 
    plot_model(b0_non, b1_non, Cp_non, nu_non, [0 0.5 0], '-');
    plot(100, 1, 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    xlabel('Relative Humidity (%)');
    ylabel('\gamma_w (Ionic)');
    title(sprintf('Non-1:1 Salts Only (e.g., MgCl2, Na2SO4, Mean \\nu=%.2f)', nu_non), ...
        'FontWeight', 'bold', 'Color', [0 0.5 0]);
    text(5, 1.05, sprintf('RMSE: %.4f\n\\beta_0: %.4f\n\\beta_1: %.4f\nC_{\\phi}: %.4f', ...
        rmse_non, b0_non, b1_non, Cp_non), 'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', [0 0.5 0]);
    xlim([0 100]);

    exportgraphics(gcf, fullfile(outdir, 'Pitzer_Subplots_Analysis.png'), 'Resolution', 300);
end

function sse = pitzer_objective(params, rh_vec, gamma_obs, m_vec, nu_vec, T)
    b0 = params(1); b1 = params(2); Cp = params(3);
    sse = 0;
    count = 0;
    
    for i = 1:length(m_vec)
        m = m_vec(i);
        nu = nu_vec(i);
        
        [gamma_calc, ~] = calculate_pitzer_point(m, nu, T, b0, b1, Cp);
        
        if ~isnan(gamma_calc)
            sse = sse + (gamma_calc - gamma_obs(i))^2;
            count = count + 1;
        end
    end
    if count == 0, sse = 1e9; end
end

function [gamma_calc, aw_calc] = calculate_pitzer_point(m, nu, T, b0, b1, Cp)
    % Calculates Gamma and Aw for a single point
    % Constants
    MWw_kg = 0.018015;
    b = 1.2;
    alpha1 = 2.0; alpha2 = 0; beta2 = 0;
    
    % Ionic Strength Estimate (Approximate based on nu)
    % For 1:1, I=m. For 2:1, I=3m.
    if nu == 2, I = m; zM=1; zX=1; vM=1; vX=1;
    elseif nu == 3, I = 3*m; zM=2; zX=1; vM=1; vX=2; % Assume M X2 type
    else, I = 0.5*nu*m; zM=1; zX=1; vM=1; vX=1; % Fallback
    end
    
    sqrt_I = sqrt(I);
    
    % Debye-Huckel
    A_phi = 0.3915 * (1 + 0.0018 * (T + 273.15 - 298.15));
    term_DH = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    
    % Pitzer Terms
    B_phi = b0 + b1 * exp(-alpha1 * sqrt_I);
    
    % Stoich factors
    z_factor = abs(zM * zX);
    factor_B = (2 * vM * vX) / nu;
    factor_C = (2 * (vM * vX)^1.5) / nu;
    
    % Osmotic Coeff Phi
    phi = 1 + z_factor * term_DH + m * factor_B * B_phi + m^2 * factor_C * Cp;
    
    % Water Activity
    ln_aw = -phi * nu * m * MWw_kg;
    aw_calc = exp(ln_aw);
    
    % Gamma Calculation
    % gamma = aw / x_w_ion
    % x_w_ion = 1 / (1 + nu * m * MWw_kg)
    x_w_ion = 1 / (1 + nu * m * MWw_kg);
    gamma_calc = aw_calc / x_w_ion;
end

function plot_original_curves(all_salts_data, salt_names, colors, fig_out_dir)
    num_salts = length(salt_names);
    
    % 1. Gamma vs Mol Fraction (Ionic)
    figure('Position', [100, 100, 1200, 800]); hold on; grid on; box on;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        plot(data.x_water_ion, data.gamma_w_ion, 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        text(data.x_water_ion(1), data.gamma_w_ion(1), ['  ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    plot([0 1], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal');
    xlabel('Mole Fraction of Water (Ionic)', 'FontSize', 14); ylabel('Activity Coefficient \gamma_w', 'FontSize', 14);
    title('Water Activity Coefficient (Ionic) vs Mole Fraction', 'FontSize', 16);
    xlim([0.0 1.0]);
    print(fullfile(fig_out_dir, 'Activity_Coefficient_Ionic_vs_MoleFraction'), '-dpng', '-r300');
    
    % 2. Gamma vs RH
    figure('Position', [150, 150, 1200, 800]); hold on; grid on; box on;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        plot(data.RH*100, data.gamma_w_ion, 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        text(data.RH(1)*100, data.gamma_w_ion(1), ['  ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal');
    xlabel('Relative Humidity (%)', 'FontSize', 14); ylabel('Activity Coefficient \gamma_w', 'FontSize', 14);
    title('Water Activity Coefficient (Ionic) vs RH', 'FontSize', 16);
    xlim([0 100]);
    print(fullfile(fig_out_dir, 'Activity_Coefficient_Ionic_vs_RH'), '-dpng', '-r300');

    % 2b. Gamma vs RH and x_w vs RH
    figure('Position', [150, 150, 1200, 800]); hold on; grid on; box on;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        plot(data.RH*100, data.gamma_w_ion, 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        plot(data.RH*100, data.x_water_ion, 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        text(data.RH(1)*100, data.gamma_w_ion(1), ['  ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    plot([0 100], [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal');
    xlabel('Relative Humidity (%)', 'FontSize', 14); ylabel('Activity Coefficient $\gamma_w$ or Mole Fraction $x_w$','FontSize', 14, 'Interpreter', 'latex');
    title('Water Activity Coefficient (Ionic) and Mole Fraction vs RH', 'FontSize', 16);
    xlim([0 100]);
    print(fullfile(fig_out_dir, 'Activity_Coefficient_and_Mol_Ionic_vs_RH'), '-dpng', '-r300');

    % 2c. Gamma vs RH and x_w vs RH ln
    figure('Position', [150, 150, 1200, 800]); hold on; grid on; box on;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        plot(data.RH*100, log(data.gamma_w_ion), 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        plot(data.RH*100, log(data.x_water_ion), 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        text((data.RH(1)*100), log(data.gamma_w_ion(1)), ['  ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    plot([0 100], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal');
    xlabel('Relative Humidity (%)', 'FontSize', 14); ylabel('ln(Activity Coefficient $\gamma_w$) or ln(Mole Fraction $x_w$)','FontSize', 14, 'Interpreter', 'latex');
    title('ln(Water Activity Coefficient (Ionic)) and ln(Water Mole Fraction (Ionic))  vs RH', 'FontSize', 16);
    % xlim([0 100]);
    print(fullfile(fig_out_dir, 'ln_Activity_Coefficient_and_ln_Mol_Ionic_vs_RH'), '-dpng', '-r300');

    % 2d. Gamma vs RH and x_w vs RH ln
    figure('Position', [150, 150, 1200, 800]); hold on; grid on; box on;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        plot(log(data.RH), log(data.gamma_w_ion), 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        plot(log(data.RH), log(data.x_water_ion), 'LineWidth', 2.5, 'color', colors(s,:), 'DisplayName', data.display_name);
        text(log(data.RH(1)), log(data.gamma_w_ion(1)), ['  ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
    end
    plot([-2.5 0], [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal');
    xlabel('ln (Relative Humidity (Fractional)', 'FontSize', 14); ylabel('ln(Activity Coefficient $\gamma_w$) or ln(Mole Fraction $x_w$)','FontSize', 14, 'Interpreter', 'latex');
    title('ln(Water Activity Coefficient (Ionic)) and ln(Water Mole Fraction (Ionic))  vs ln(RH)', 'FontSize', 16);
    % xlim([0 100]);
    print(fullfile(fig_out_dir, 'ln_Activity_Coefficient_and_ln_Mol_Ionic_vs_ln_aw'), '-dpng', '-r300');
    
    % 3. Scatter at 90% RH
    figure('Position', [200, 200, 900, 700]); hold on; grid on; box on;
    target_RH = 0.90;
    for s = 1:num_salts
        data = all_salts_data.(salt_names{s});
        rh_lower = min(data.RH);
        gamma_at_90 = interp1(data.RH, data.gamma_w_ion, target_RH, 'linear', NaN);
        if ~isnan(gamma_at_90)
            scatter(rh_lower * 100, gamma_at_90, 80, colors(s,:), 'filled', 'MarkerEdgeColor', 'k');
            text(rh_lower * 100, gamma_at_90+0.002, [' ' data.display_name], 'Color', colors(s,:), 'FontSize', 9, 'FontWeight', 'bold');
        end
    end
    yline(1.0, 'k--', 'Ideal');
    xlabel('Lower Bound RH (%)'); ylabel(['\gamma_w at ' num2str(target_RH*100) '% RH']);
    title('Screening: Non-Ideality vs Deliquescence', 'FontSize', 14);
    print(fullfile(fig_out_dir, 'Gamma90_vs_LowerBoundRH'), '-dpng', '-r300');
end

function colors = generate_colors(num_salts)
    colors = zeros(num_salts, 3);
    golden_angle = 0.38196601125; 
    for i = 1:num_salts
        hue = mod((i-1) * golden_angle, 1.0);
        sat = 0.85 + 0.15 * sin(i * 0.5); 
        val = 0.85 + 0.15 * cos(i * 0.4);
        colors(i,:) = hsv2rgb([hue, sat, val]);
    end
end

function plot_pitzer_curve(b0, b1, cp, nu_val, max_m, T, col, style)
    m_plot = linspace(0.1, max_m, 100);
    g_curve = zeros(size(m_plot));
    rh_curve = zeros(size(m_plot));
    for i=1:length(m_plot)
        [g_curve(i), aw_tmp] = calculate_pitzer_point(m_plot(i), nu_val, T, b0, b1, cp);
        rh_curve(i) = aw_tmp * 100;
    end
    plot(rh_curve, g_curve, 'Color', col, 'LineStyle', style, 'LineWidth', 2.5);
end