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
