function aw = calculate_activity_HCl_NaCl(mf1, mf2)
% Calculate water activity for HCl + NaCl
% Component 1: HCl
% Component 2: NaCl
%
% Inputs:
%   mf1 - mass fraction of HCl
%   mf2 - mass fraction of NaCl
% Output:
%   aw - water activity
%
% Data range:
%   HCl: 0.0172 to 0.0649
%   NaCl: 0.0210 to 0.1347
%   aw: 0.7910 to 0.9650

% Validate inputs
if mf1 < 0.017175 || mf1 > 0.064882
    warning('HCl mass fraction outside calibrated range');
end
if mf2 < 0.021025 || mf2 > 0.134668
    warning('NaCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000000

% Calculate water activity
aw = 0 + 1.3097712580e+00 + -3.4671003106e+01*mf1 + 7.8094147486e-01*mf2 + 9.3515567463e+02*mf1.^2 + 2.6893621615e+01*mf1.*mf2 + -2.5334324019e+01*mf2.^2 + -7.2748263120e+03*mf1.^3 + -7.9620996064e+02*mf1.^2.*mf2 + 8.9813275573e+01*mf1.*mf2.^2 + 9.1022668532e+01*mf2.^3;

end
