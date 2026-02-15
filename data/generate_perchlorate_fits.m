% Generate polynomial fits for Ca(ClO4)2, Sr(ClO4)2, and Ba(ClO4)2
% using Pitzer model with parameters from the thermodynamic database
%
% This script:
% 1. Uses Pitzer parameters to calculate water activity
% 2. Generates mass fraction vs RH data points
% 3. Fits 4th degree polynomials
% 4. Outputs coefficients for calculate_mf functions

clear; clc; close all;

% Physical constants
R = 8.314462;  % J/(mol·K)
T = 298.15;    % K (25°C)
MWw = 18.015;  % g/mol (water)

% Debye-Hückel parameters at 25°C
A_phi = 0.392;  % Debye-Hückel parameter for osmotic coefficient

%% Salt definitions
salts = {
    % Name, MW (g/mol), z_cat, z_an, nu, Pitzer parameters [beta0, beta1, beta2, Cphi]
    'Ca(ClO4)2', 238.98, 2, -1, 3, [0.4511, 1.757, 0.0, -0.005];
    'Ba(ClO4)2', 336.23, 2, -1, 3, [0.3614, 1.5758, 0.0, -0.03126];
    'Sr(ClO4)2', 286.52, 2, -1, 3, [0.427, 1.7, 0.0, -0.02];  % Estimated from Ca and Ba
};

%% Generate fits for each salt
for s = 1:size(salts, 1)
    salt_name = salts{s, 1};
    MW = salts{s, 2};
    z_cat = salts{s, 3};
    z_an = salts{s, 4};
    nu = salts{s, 5};
    pitzer = salts{s, 6};
    
    beta0 = pitzer(1);
    beta1 = pitzer(2);
    beta2 = pitzer(3);
    Cphi = pitzer(4);
    
    alpha1 = 2.0;   % Standard value for 2-1 electrolytes
    alpha2 = 0.0;   % Not used when beta2 = 0
    
    fprintf('\n==================================================\n');
    fprintf('Processing: %s\n', salt_name);
    fprintf('==================================================\n');
    fprintf('Pitzer parameters: beta0=%.4f, beta1=%.4f, Cphi=%.4f\n', beta0, beta1, Cphi);
    
    %% Generate molality range
    % For 2-1 electrolytes, typical range is 0.1 to saturation (~ 3-6 mol/kg)
    % Let's generate points from dilute to concentrated solutions
    m_vec = linspace(0.01, 5.0, 50)';  % molality in mol/kg
    
    %% Calculate ionic strength
    % For M(ClO4)2: I = 1/2 * (2^2 * m + 1^2 * 2*m) = 3m
    I_vec = 3 * m_vec;
    
    %% Calculate water activity using Pitzer model
    aw_vec = zeros(size(m_vec));
    
    for i = 1:length(m_vec)
        m = m_vec(i);
        I = I_vec(i);
        
        % Debye-Hückel term
        f_gamma = -A_phi * (sqrt(I) / (1 + 1.2*sqrt(I)));
        
        % g functions for Pitzer model
        x = alpha1 * sqrt(I);
        if x < 1e-4
            g_beta1 = 0;
        else
            g_beta1 = 2 * (1 - (1 + x) * exp(-x)) / (x^2);
        end
        
        % BMX term (for 2-1 electrolyte)
        BMX = beta0 + beta1 * g_beta1 + beta2 * 0;  % beta2 = 0 for these salts
        
        % Osmotic coefficient
        nu_M = 1;  % number of cations
        nu_X = 2;  % number of anions
        nu_total = nu_M + nu_X;
        z_M = abs(z_cat);
        z_X = abs(z_an);
        
        phi = 1 + f_gamma + ...
              m * (2 * nu_M * nu_X / nu_total) * (BMX + m * z_M * z_X * Cphi);
        
        % Water activity
        % aw = exp(-phi * nu * m * MWw / 1000)
        aw = exp(-phi * nu_total * m * MWw / 1000);
        aw_vec(i) = aw;
    end
    
    %% Convert to mass fraction
    % m = mol salt / kg water
    % mass of salt = m * MW (in g)
    % mass fraction = mass_salt / (mass_salt + mass_water)
    %               = (m * MW) / (m * MW + 1000)
    mf_vec = (m_vec * MW) ./ (m_vec * MW + 1000);
    
    %% Filter valid range (RH typically 0.3 to 0.99)
    valid_idx = (aw_vec >= 0.3) & (aw_vec <= 0.99);
    mf_valid = mf_vec(valid_idx);
    aw_valid = aw_vec(valid_idx);
    
    fprintf('Valid RH range: %.4f to %.4f\n', min(aw_valid), max(aw_valid));
    fprintf('Valid mf range: %.4f to %.4f\n', min(mf_valid), max(mf_valid));
    fprintf('Number of data points: %d\n', length(mf_valid));
    
    if length(mf_valid) < 10
        warning('Too few data points for %s, adjusting molality range...', salt_name);
        continue;
    end
    
    %% Fit 4th degree polynomial: RH = a0 + a1*mf + a2*mf^2 + a3*mf^3 + a4*mf^4
    % This matches the form used in calculate_mf functions
    p = polyfit(mf_valid, aw_valid, 4);
    
    % Extract coefficients
    A_4 = p(1);
    A_3 = p(2);
    A_2 = p(3);
    A_1 = p(4);
    A_0 = p(5);
    
    % Calculate R^2
    aw_fit = polyval(p, mf_valid);
    SS_res = sum((aw_valid - aw_fit).^2);
    SS_tot = sum((aw_valid - mean(aw_valid)).^2);
    R_squared = 1 - SS_res / SS_tot;
    
    % Calculate RMSE
    RMSE = sqrt(mean((aw_valid - aw_fit).^2));
    
    fprintf('\nPolynomial coefficients (RH = A_0 + A_1*mf + ... + A_4*mf^4):\n');
    fprintf('A_4 = %.4f\n', A_4);
    fprintf('A_3 = %.4f\n', A_3);
    fprintf('A_2 = %.4f\n', A_2);
    fprintf('A_1 = %.4f\n', A_1);
    fprintf('A_0 = %.4f\n', A_0);
    fprintf('\nFit quality:\n');
    fprintf('R^2 = %.6f\n', R_squared);
    fprintf('RMSE = %.6f\n', RMSE);
    
    %% Skip plotting in non-GUI mode
    % Figure generation commented out for command-line execution
    % Uncomment below if running in GUI mode
    % figure('Position', [100, 100, 800, 600]);
    % subplot(2,1,1);
    % plot(mf_valid, aw_valid, 'bo', 'MarkerSize', 6, 'DisplayName', 'Pitzer Model');
    % hold on;
    % mf_plot = linspace(min(mf_valid), max(mf_valid), 200);
    % aw_plot = polyval(p, mf_plot);
    % plot(mf_plot, aw_plot, 'r-', 'LineWidth', 2, 'DisplayName', '4th Order Fit');
    % xlabel('Mass Fraction');
    % ylabel('Relative Humidity (Water Activity)');
    % title(sprintf('%s: RH vs Mass Fraction', strrep(salt_name, '_', '\_')));
    % legend('Location', 'best');
    % grid on;
    % 
    % subplot(2,1,2);
    % residuals = aw_valid - aw_fit;
    % plot(mf_valid, residuals, 'ko', 'MarkerSize', 4);
    % xlabel('Mass Fraction');
    % ylabel('Residuals');
    % title(sprintf('Fit Residuals (RMSE = %.6f)', RMSE));
    % grid on;
    % yline(0, 'r--');
    % 
    % % Save figure
    % saveas(gcf, sprintf('%s_polynomial_fit.png', strrep(salt_name, '(', ''), strrep(salt_name, ')', '')));
    
    %% Generate MATLAB code snippet for calculate_mf function
    fprintf('\n--- Code snippet for calculate_mf_%s.m ---\n', ...
        strrep(strrep(salt_name, '(', ''), ')', '2'));
    fprintf('if RH < %.4f || RH > %.4f \n', min(aw_valid), max(aw_valid));
    fprintf('    error("below deliquescence relative humidity or above range") \n');
    fprintf('end\n');
    fprintf('A_4 = %.4f;\n', A_4);
    fprintf('A_3 = %.4f;\n', A_3);
    fprintf('A_2 = %.4f;\n', A_2);
    fprintf('A_1 = %.4f;\n', A_1);
    fprintf('A_0 = %.4f;\n', A_0);
    fprintf('f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;\n');
    fprintf('mf = robust_fzero(f, %.4f, %.4f, %.4f);\n', ...
        min(mf_valid), max(mf_valid), mean(mf_valid));
    fprintf('--- End code snippet ---\n\n');
    
    %% Store data for export
    results(s).name = salt_name;
    results(s).MW = MW;
    results(s).RH_min = min(aw_valid);
    results(s).RH_max = max(aw_valid);
    results(s).mf_min = min(mf_valid);
    results(s).mf_max = max(mf_valid);
    results(s).A = [A_0, A_1, A_2, A_3, A_4];
    results(s).R_squared = R_squared;
    results(s).RMSE = RMSE;
    results(s).data_mf = mf_valid;
    results(s).data_aw = aw_valid;
end

%% Summary
fprintf('\n\n');
fprintf('=========================================================\n');
fprintf('SUMMARY OF POLYNOMIAL FITS\n');
fprintf('=========================================================\n');
fprintf('%-15s %10s %10s %10s %8s\n', 'Salt', 'RH_min', 'RH_max', 'mf_max', 'R^2');
fprintf('---------------------------------------------------------\n');
for s = 1:length(results)
    fprintf('%-15s %10.4f %10.4f %10.4f %8.6f\n', ...
        results(s).name, results(s).RH_min, results(s).RH_max, ...
        results(s).mf_max, results(s).R_squared);
end
fprintf('=========================================================\n');

fprintf('\nPolynomial fits complete! Update the calculate_mf files with the coefficients shown above.\n');
