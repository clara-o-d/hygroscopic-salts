function aw = pitzer_water_activity(m, nu, beta0, beta1, beta2, Cphi, alpha1, alpha2, T)
% PITZER_WATER_ACTIVITY Calculate water activity using corrected Pitzer model
% Reference: Pitzer (1973), Eq. 13 for Osmotic Coefficient

    if nargin < 9, error('All 9 input arguments are required'); end

    % --- 1. Constants & Pre-calculations ---
    Mw = 0.018015; % kg/mol
    b  = 1.2;      % Standard Pitzer parameter
    
    % Ionic Strength (I)
    % Assuming simple single salt M(vM)X(vX)
    % 1:1 salt (NaCl) -> I = m
    % 2:1 salt (MgCl2)-> I = 3m
    % We approximate I using your input 'nu':
    if nu == 2      % 1:1
        I = m;
        vM = 1; vX = 1;
        zM = 1; zX = 1;
    elseif nu == 3  % 2:1 or 1:2
        I = 3 * m;
        vM = 1; vX = 2; % Assuming M X2 type generally
        zM = 2; zX = 1;
    else            % Fallback (approx)
        I = 0.5 * nu * m; 
        vM = 1; vX = 1; zM = 1; zX = 1; 
    end

    sqrt_I = sqrt(I);
    A_phi = debye_huckel_parameter(T);

    % --- 2. Debye-Huckel Term (Osmotic) ---
    % CORRECTED: This is f_phi, NOT f_gamma. No log term.
    term_DH = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    
    % --- 3. Second Virial Coefficient (B_phi) ---
    % Formula: Beta0 + Beta1 * exp(-alpha1 * sqrt(I)) ...
    exp1 = exp(-alpha1 * sqrt_I);
    B_phi = beta0 + beta1 * exp1;
    
    if abs(beta2) > 1e-9 || alpha2 > 1e-9
        exp2 = exp(-alpha2 * sqrt_I);
        B_phi = B_phi + beta2 * exp2;
    end

    % --- 4. Stoichiometric Multipliers ---
    % The Pitzer equation for Phi is:
    % Phi - 1 = |zM*zX| * f_phi 
    %           + m * (2*vM*vX/nu) * B_phi 
    %           + m^2 * (2*(vM*vX)^1.5/nu) * C_phi
    
    % Charge factor for DH term
    z_factor = abs(zM * zX); 
    
    % Factor for B term: (2 * vM * vX) / nu
    % 1:1 -> 2*1*1 / 2 = 1
    % 2:1 -> 2*1*2 / 3 = 1.333
    factor_B = (2 * vM * vX) / nu;
    
    % Factor for C term: (2 * (vM*vX)^1.5) / nu
    % 1:1 -> 2 * 1 / 2 = 1
    % 2:1 -> 2 * (2)^1.5 / 3 = 1.885
    factor_C = (2 * (vM * vX)^1.5) / nu;

    % --- 5. Calculate Osmotic Coefficient (Phi) ---
    phi = 1 + z_factor * term_DH + m * factor_B * B_phi + m^2 * factor_C * Cphi;
    
    % --- 6. Calculate Water Activity ---
    % ln(aw) = -phi * nu * m * Mw
    ln_aw = -phi * nu * m * Mw;
    aw = exp(ln_aw);
    
    % Physical Bounds
    aw = max(0, min(1, aw));
end

function A_phi = debye_huckel_parameter(T)
    % A_phi approx at 25C is 0.392 (Pitzer 1973 used 0.392)
    % Temperature dependence
    T_K = T + 273.15;
    T_ref = 298.15;
    A_phi_ref = 0.3915; % Slightly more precise ref
    A_phi = A_phi_ref * (1 + 0.0018 * (T_K - T_ref)); 
end


% function aw = pitzer_water_activity(m, nu, beta0, beta1, beta2, Cphi, alpha1, alpha2, T)
% % PITZER_WATER_ACTIVITY Calculate water activity using Pitzer model
% %
% % Inputs:
% %   m       - molality (mol salt/kg water)
% %   nu      - number of ions (2 for 1:1, 3 for 2:1, etc.)
% %   beta0   - Pitzer parameter beta0
% %   beta1   - Pitzer parameter beta1
% %   beta2   - Pitzer parameter beta2 (use 0 for 1:1 electrolytes)
% %   Cphi    - Pitzer parameter Cphi
% %   alpha1  - Pitzer parameter alpha1 (typically 2.0)
% %   alpha2  - Pitzer parameter alpha2 (typically 12.0 for 2:1, 0 for 1:1)
% %   T       - Temperature in Celsius
% %
% % Output:
% %   aw      - water activity
% %
% % Reference: Pitzer, K.S. (1973). "Thermodynamics of electrolytes. I. 
% %            Theoretical basis and general equations." 
% %            J. Phys. Chem., 77(2), 268-277.
% %
% % The Pitzer model relates water activity to the osmotic coefficient:
% %   ln(aw) = -phi * nu * m * Mw
% % where Mw = 0.018015 kg/mol (molecular weight of water)
% 
% if nargin < 9
%     error('All 9 input arguments are required')
% end
% 
% % Constants
% Mw = 0.018015; % kg/mol - molecular weight of water
% T_K = T + 273.15; % Convert to Kelvin
% 
% % Debye-Hückel limiting slope (Aphi)
% % Temperature-dependent formula (valid for water at different temperatures)
% % At 25°C, Aphi ≈ 0.392
% A_phi = debye_huckel_parameter(T);
% 
% % Ionic strength
% I = 0.5 * nu * m; % For simple salts like MX or MX2
% 
% % Calculate ionic charge product
% % For MX (1:1): z_M = 1, z_X = 1, |z_M*z_X| = 1
% % For MX2 (2:1): z_M = 2, z_X = 1, |z_M*z_X| = 2
% if nu == 2
%     z_product = 1; % 1:1 electrolyte
% elseif nu == 3
%     z_product = 2; % 2:1 or 1:2 electrolyte
% else
%     z_product = 1; % Default
% end
% 
% % Calculate b parameter (typically 1.2 kg^0.5/mol^0.5)
% b = 1.2;
% 
% % Calculate f(I) - Debye-Hückel term
% sqrt_I = sqrt(I);
% f_gamma = -A_phi * (sqrt_I / (1 + b * sqrt_I) + (2/b) * log(1 + b * sqrt_I));
% 
% % Calculate B_phi and B'_phi
% % For 1:1 electrolytes (alpha2 = 0, beta2 = 0)
% % For 2:1 electrolytes (use both terms)
% 
% if alpha2 == 0 % 1:1 electrolyte
%     g_alpha1 = 2 * (1 - (1 + alpha1 * sqrt_I) * exp(-alpha1 * sqrt_I)) / (alpha1^2 * I);
%     B_phi = beta0 + beta1 * g_alpha1;
% else % 2:1 electrolyte
%     g_alpha1 = 2 * (1 - (1 + alpha1 * sqrt_I) * exp(-alpha1 * sqrt_I)) / (alpha1^2 * I);
%     g_alpha2 = 2 * (1 - (1 + alpha2 * sqrt_I) * exp(-alpha2 * sqrt_I)) / (alpha2^2 * I);
%     B_phi = beta0 + beta1 * g_alpha1 + beta2 * g_alpha2;
% end
% 
% % Calculate osmotic coefficient
% % phi = 1 + f_gamma + m * B_phi + m^2 * C_phi
% % For mixed electrolytes, additional terms would be needed
% 
% phi = 1 + f_gamma + m * B_phi + m^2 * Cphi;
% 
% % Handle negative or extremely small osmotic coefficients
% if phi < 0.01
%     phi = 0.01;
%     warning('Osmotic coefficient calculated as negative or very small. Setting to 0.01')
% end
% 
% % Calculate water activity
% % ln(aw) = -phi * nu * m * Mw
% ln_aw = -phi * nu * m * Mw;
% aw = exp(ln_aw);
% 
% % Ensure water activity is physical (between 0 and 1)
% if aw > 1
%     aw = 1;
% elseif aw < 0
%     aw = 0;
% end
% 
% end
% 
% 
% function A_phi = debye_huckel_parameter(T)
% % Calculate Debye-Hückel parameter A_phi as a function of temperature
% % for water as solvent
% %
% % Input:
% %   T - Temperature in Celsius
% %
% % Output:
% %   A_phi - Debye-Hückel limiting slope
% %
% % Reference: Archer & Wang (1990) and Pitzer (1991)
% 
% T_K = T + 273.15; % Convert to Kelvin
% 
% % Simplified temperature-dependent formula
% % At 25°C (298.15 K): A_phi ≈ 0.392
% % At 0°C (273.15 K): A_phi ≈ 0.377
% % At 100°C (373.15 K): A_phi ≈ 0.465
% 
% % Linear approximation (good for 0-100°C range)
% % More accurate formula would use density and dielectric constant of water
% T_ref = 298.15; % Reference temperature (25°C)
% A_phi_ref = 0.392; % Value at 25°C
% 
% % Temperature dependence (simplified)
% A_phi = A_phi_ref * (1 + 0.0018 * (T_K - T_ref));
% 
% end
