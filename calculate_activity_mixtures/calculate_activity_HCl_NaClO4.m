function aw = calculate_activity_HCl_NaClO4(mf1, mf2)
% Calculate water activity for HCl + NaClO4
% Component 1: HCl
% Component 2: NaClO4
%
% Inputs:
%   mf1 - mass fraction of HCl
%   mf2 - mass fraction of NaClO4
% Output:
%   aw - water activity
%
% Data range:
%   HCl: 0.0182 to 0.0182
%   NaClO4: 0.0590 to 0.6272
%   aw: 0.5250 to 0.9640

% Validate inputs
if mf1 < 0.018227 || mf1 > 0.018227
    warning('HCl mass fraction outside calibrated range');
end
if mf2 < 0.058959 || mf2 > 0.627189
    warning('NaClO4 mass fraction outside calibrated range');
end

% Polynomial fit using NaClO4 mass fraction
% Degree: 4, RMSE: 0.002429
A_4 = 1.1234801987e+01;
A_3 = -1.3998936307e+01;
A_2 = 4.6762351093e+00;
A_1 = -9.7024572631e-01;
A_0 = 1.0088993542e+00;

% Calculate water activity
mf = mf2;  % Using NaClO4
aw = A_0 + A_1*mf^1 + A_2*mf^2 + A_3*mf^3 + A_4*mf^4;

end
