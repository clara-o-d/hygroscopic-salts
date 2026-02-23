function aw = calculate_activity_HCl_Na2SO4(mf1, mf2)
% Calculate water activity for HCl + Na2SO4
% Component 1: HCl
% Component 2: Na2SO4
%
% Inputs:
%   mf1 - mass fraction of HCl
%   mf2 - mass fraction of Na2SO4
% Output:
%   aw - water activity
%
% Data range:
%   HCl: 0.0036 to 0.0344
%   Na2SO4: 0.0141 to 0.2073
%   aw: 0.9030 to 0.9940

% Validate inputs
if mf1 < 0.003582 || mf1 > 0.034369
    warning('HCl mass fraction outside calibrated range');
end
if mf2 < 0.014074 || mf2 > 0.207290
    warning('Na2SO4 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000218

% Calculate water activity
aw = 0 + 1.0010319232e+00 + -9.0373273068e-01*mf1 + -2.6942708615e-01*mf2 + -4.0929238877e+00*mf1.^2 + 1.5212581198e+00*mf1.*mf2 + -2.8144856520e-01*mf2.^2 + -2.2408140475e-01*mf1.^3 + 8.2019042730e+00*mf1.^2.*mf2 + -6.6660341974e+00*mf1.*mf2.^2 + 3.3926417276e-01*mf2.^3;

end
