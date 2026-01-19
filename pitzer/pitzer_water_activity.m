function aw = pitzer_water_activity(m, nu, beta0, beta1, beta2, Cphi, alpha1, alpha2, T)
% PITZER_WATER_ACTIVITY Calculate water activity using Pitzer model
%
% Inputs:
%   m       - molality (mol salt/kg water)
%   nu      - number of ions (2 for 1:1, 3 for 2:1, etc.)
%   beta0   - Pitzer parameter beta0
%   beta1   - Pitzer parameter beta1
%   beta2   - Pitzer parameter beta2 (use 0 for 1:1 electrolytes)
%   Cphi    - Pitzer parameter Cphi
%   alpha1  - Pitzer parameter alpha1 (typically 2.0)
%   alpha2  - Pitzer parameter alpha2 (typically 12.0 for 2:1, 0 for 1:1)
%   T       - Temperature in Celsius
%
% Output:
%   aw      - water activity
%
% Reference: Pitzer, K.S. (1973). "Thermodynamics of electrolytes. I. 
%            Theoretical basis and general equations." 
%            J. Phys. Chem., 77(2), 268-277.
%
% The Pitzer model relates water activity to the osmotic coefficient:
%   ln(aw) = -phi * nu * m * Mw
% where Mw = 0.018015 kg/mol (molecular weight of water)

if nargin < 9
    error('All 9 input arguments are required')
end

% Constants
Mw = 0.018015; % kg/mol - molecular weight of water
T_K = T + 273.15; % Convert to Kelvin

% Debye-Hückel limiting slope (Aphi)
% Temperature-dependent formula (valid for water at different temperatures)
% At 25°C, Aphi ≈ 0.392
A_phi = debye_huckel_parameter(T);

% Ionic strength
I = 0.5 * nu * m; % For simple salts like MX or MX2

% Calculate ionic charge product
% For MX (1:1): z_M = 1, z_X = 1, |z_M*z_X| = 1
% For MX2 (2:1): z_M = 2, z_X = 1, |z_M*z_X| = 2
if nu == 2
    z_product = 1; % 1:1 electrolyte
elseif nu == 3
    z_product = 2; % 2:1 or 1:2 electrolyte
else
    z_product = 1; % Default
end

% Calculate b parameter (typically 1.2 kg^0.5/mol^0.5)
b = 1.2;

% Calculate f(I) - Debye-Hückel term
sqrt_I = sqrt(I);
f_gamma = -A_phi * (sqrt_I / (1 + b * sqrt_I) + (2/b) * log(1 + b * sqrt_I));

% Calculate B_phi and B'_phi
% For 1:1 electrolytes (alpha2 = 0, beta2 = 0)
% For 2:1 electrolytes (use both terms)

if alpha2 == 0 % 1:1 electrolyte
    g_alpha1 = 2 * (1 - (1 + alpha1 * sqrt_I) * exp(-alpha1 * sqrt_I)) / (alpha1^2 * I);
    B_phi = beta0 + beta1 * g_alpha1;
else % 2:1 electrolyte
    g_alpha1 = 2 * (1 - (1 + alpha1 * sqrt_I) * exp(-alpha1 * sqrt_I)) / (alpha1^2 * I);
    g_alpha2 = 2 * (1 - (1 + alpha2 * sqrt_I) * exp(-alpha2 * sqrt_I)) / (alpha2^2 * I);
    B_phi = beta0 + beta1 * g_alpha1 + beta2 * g_alpha2;
end

% Calculate osmotic coefficient
% phi = 1 + f_gamma + m * B_phi + m^2 * C_phi
% For mixed electrolytes, additional terms would be needed

phi = 1 + f_gamma + m * B_phi + m^2 * Cphi;

% Handle negative or extremely small osmotic coefficients
if phi < 0.01
    phi = 0.01;
    warning('Osmotic coefficient calculated as negative or very small. Setting to 0.01')
end

% Calculate water activity
% ln(aw) = -phi * nu * m * Mw
ln_aw = -phi * nu * m * Mw;
aw = exp(ln_aw);

% Ensure water activity is physical (between 0 and 1)
if aw > 1
    aw = 1;
elseif aw < 0
    aw = 0;
end

end


function A_phi = debye_huckel_parameter(T)
% Calculate Debye-Hückel parameter A_phi as a function of temperature
% for water as solvent
%
% Input:
%   T - Temperature in Celsius
%
% Output:
%   A_phi - Debye-Hückel limiting slope
%
% Reference: Archer & Wang (1990) and Pitzer (1991)

T_K = T + 273.15; % Convert to Kelvin

% Simplified temperature-dependent formula
% At 25°C (298.15 K): A_phi ≈ 0.392
% At 0°C (273.15 K): A_phi ≈ 0.377
% At 100°C (373.15 K): A_phi ≈ 0.465

% Linear approximation (good for 0-100°C range)
% More accurate formula would use density and dielectric constant of water
T_ref = 298.15; % Reference temperature (25°C)
A_phi_ref = 0.392; % Value at 25°C

% Temperature dependence (simplified)
A_phi = A_phi_ref * (1 + 0.0018 * (T_K - T_ref));

end
