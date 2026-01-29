function fit_pitzer_from_exothermic()
    % Main execution function
    clc; clear; close all;

    % ================= CONFIGURATION =================
    % 1. Define the input polynomial (Activity Coeff (y) vs Activity (x))
    %    y = -0.000029 x^2 + 0.012902 x + 0.000849
    poly_func = @(x) -0.000029.*x.^2 + 0.012902.*x + 0.000849;
    
    % 2. System Constants
    nu = 2;         % 1:1 Electrolyte (e.g., NaCl) -> nu = 2
    T  = 25;        % Temperature in Celsius
    Mw = 0.018015;  % Molar mass of water (kg/mol)
    
    % 3. Optimization Settings
    %    Initial guesses for [beta0, beta1, Cphi]
    initial_params = [0.1, 0.1, 0.001]; 
    
    % ================= DATA GENERATION =================
    fprintf('Generating data from polynomial...\n');
    
    % We assume x is Percent RH (0-100) based on the coefficients provided.
    % If x is strictly 0-1, change this range to linspace(0.1, 0.99, 50).
    x_gen = linspace(10, 95, 50); 
    
    y_gen = poly_func(x_gen); % Gamma (Activity Coefficient)
    
    % Convert to physical units (0-1 scale)
    % Assumption: x is %RH, so aw = x / 100.
    aw_data = x_gen ./ 100; 
    gamma_data = y_gen;
    
    % Calculate Mole Fraction of Water (xw)
    % aw = xw * gamma  =>  xw = aw / gamma
    xw_data = aw_data ./ gamma_data;
    
    % Filter invalid physics (xw must be <= 1)
    valid_indices = xw_data < 1 & xw_data > 0;
    xw_data = xw_data(valid_indices);
    aw_data = aw_data(valid_indices);
    gamma_data = gamma_data(valid_indices);
    
    if isempty(xw_data)
        error('Data generation failed: Calculated mole fractions are > 1. Check if input x is %RH (0-100) or activity (0-1).');
    end

    % Calculate Molality (m)
    % xw = nw / (nw + nu*m_salt)
    % ... algebra ...
    % m = (1/xw - 1) * (1 / (nu * Mw))
    m_data = (1./xw_data - 1) .* (1 / (nu * Mw));

    fprintf('Valid data points generated: %d\n', length(m_data));
    fprintf('Max Molality in range: %.2f mol/kg\n', max(m_data));

    % ================= OPTIMIZATION =================
    % Define the cost function (Sum of Squared Errors)
    cost_func = @(p) objective_function(p, m_data, aw_data, nu, T);

    % Run Optimization (Nelder-Mead Simplex)
    options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxFunEvals', 1000);
    fprintf('\nStarting Optimization...\n');
    p_opt = fminsearch(cost_func, initial_params, options);

    % Unpack results
    beta0_opt = p_opt(1);
    beta1_opt = p_opt(2);
    Cphi_opt  = p_opt(3);

    % ================= REPORTING =================
    fprintf('\n--- Optimization Results (1:1 Salt) ---\n');
    fprintf('Beta0 : %.5f\n', beta0_opt);
    fprintf('Beta1 : %.5f\n', beta1_opt);
    fprintf('Cphi  : %.5f\n', Cphi_opt);
    
    % ================= PLOTTING =================
    % Calculate model curve using optimized parameters
    aw_model = zeros(size(m_data));
    for i = 1:length(m_data)
        % Fixed alpha1=2.0, alpha2=0, beta2=0 for standard 1:1 Pitzer
        aw_model(i) = pitzer_water_activity(m_data(i), nu, beta0_opt, beta1_opt, 0, Cphi_opt, 2.0, 0, T);
    end

    figure('Color', 'w');
    
    subplot(2,1,1);
    plot(m_data, aw_data, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Input Data (from Poly)');
    hold on;
    plot(m_data, aw_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Pitzer Model');
    xlabel('Molality (m)');
    ylabel('Water Activity (a_w)');
    title('Pitzer Fit: Water Activity vs Molality');
    legend('Location', 'Best');
    grid on;

    subplot(2,1,2);
    % Back-calculate Gamma to compare with original polynomial y
    % gamma_model = aw_model / xw
    gamma_model = aw_model ./ xw_data; 
    plot(aw_data*100, gamma_data, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Original Poly (y)');
    hold on;
    plot(aw_model*100, gamma_model, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Pitzer Back-calc');
    xlabel('RH');
    ylabel('Activity Coefficient (\gamma_w)');
    title('Verification: Activity Coefficient Match');
    legend('Location', 'Best');
    eqn_str = sprintf([ ...
        '$\\beta_0 = %.6f$\\\\ ' ...
        '$\\beta_1 = %.6f$\\\\ ' ...
        '$C_{\\phi} = %.6f$'], ...
        beta0_opt, beta1_opt, Cphi_opt);
    
    annotation('textbox', [0.62 0.78 0.30 0.05], ...
        'String', eqn_str, ...
        'Interpreter', 'latex', ...
        'FontSize', 12, ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1.1);
    
    grid on

    % ================= SAVE FIGURE =================
    outdir = fullfile('..', 'figures', 'pitzer_comparisons');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    outfile = fullfile(outdir, 'pitzer_fit_comparison.png');
    exportgraphics(gcf, outfile, 'Resolution', 300);

end

% ================= HELPER FUNCTIONS =================

function sse = objective_function(params, m_data, aw_target, nu, T)
    % Unpack parameters
    b0 = params(1);
    b1 = params(2);
    Cp = params(3);
    
    % Standard Pitzer constants for 1:1
    alpha1 = 2.0;
    alpha2 = 0.0;
    beta2  = 0.0; 
    
    sse = 0;
    for i = 1:length(m_data)
        aw_calc = pitzer_water_activity(m_data(i), nu, b0, b1, beta2, Cp, alpha1, alpha2, T);
        
        % Error accumulation (Sum of Squared Errors)
        sse = sse + (aw_calc - aw_target(i))^2;
    end
end

% --- YOUR PROVIDED PITZER FUNCTION ---
function aw = pitzer_water_activity(m, nu, beta0, beta1, beta2, Cphi, alpha1, alpha2, T)
% PITZER_WATER_ACTIVITY Calculate water activity using corrected Pitzer model
% Reference: Pitzer (1973), Eq. 13 for Osmotic Coefficient
    % if nargin < 9, error('All 9 input arguments are required'); end
    % --- 1. Constants & Pre-calculations ---
    Mw = 0.018015; % kg/mol
    b  = 1.2;      % Standard Pitzer parameter
    
    % Ionic Strength (I)
    if nu == 2      % 1:1
        I = m;
        vM = 1; vX = 1;
        zM = 1; zX = 1;
    elseif nu == 3  % 2:1 or 1:2
        I = 3 * m;
        vM = 1; vX = 2; 
        zM = 2; zX = 1;
    else            % Fallback (approx)
        I = 0.5 * nu * m; 
        vM = 1; vX = 1; zM = 1; zX = 1; 
    end
    sqrt_I = sqrt(I);
    A_phi = debye_huckel_parameter(T);
    
    % --- 2. Debye-Huckel Term (Osmotic) ---
    term_DH = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    
    % --- 3. Second Virial Coefficient (B_phi) ---
    exp1 = exp(-alpha1 * sqrt_I);
    B_phi = beta0 + beta1 * exp1;
    
    if abs(beta2) > 1e-9 || alpha2 > 1e-9
        exp2 = exp(-alpha2 * sqrt_I);
        B_phi = B_phi + beta2 * exp2;
    end
    
    % --- 4. Stoichiometric Multipliers ---
    z_factor = abs(zM * zX); 
    factor_B = (2 * vM * vX) / nu;
    factor_C = (2 * (vM * vX)^1.5) / nu;
    
    % --- 5. Calculate Osmotic Coefficient (Phi) ---
    phi = 1 + z_factor * term_DH + m * factor_B * B_phi + m^2 * factor_C * Cphi;
    
    % --- 6. Calculate Water Activity ---
    ln_aw = -phi * nu * m * Mw;
    aw = exp(ln_aw);
    
    % Physical Bounds
    aw = max(0, min(1, aw));
end

function A_phi = debye_huckel_parameter(T)
    % A_phi approx at 25C is 0.392 (Pitzer 1973 used 0.392)
    T_K = T + 273.15;
    T_ref = 298.15;
    A_phi_ref = 0.3915; 
    A_phi = A_phi_ref * (1 + 0.0018 * (T_K - T_ref)); 
end