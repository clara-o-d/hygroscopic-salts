function aw = calculate_activity_NH4Cl_LiCl(mf1, mf2)
% Calculate water activity for NH4Cl + LiCl
% Component 1: NH4Cl
% Component 2: LiCl
%
% Inputs:
%   mf1 - mass fraction of NH4Cl
%   mf2 - mass fraction of LiCl
% Output:
%   aw - water activity
%
% Data range:
%   NH4Cl: 0.0053 to 0.1762
%   LiCl: 0.0042 to 0.1450
%   aw: 0.7220 to 0.9900

% Validate inputs
if mf1 < 0.005321 || mf1 > 0.176250
    warning('NH4Cl mass fraction outside calibrated range');
end
if mf2 < 0.004222 || mf2 > 0.144989
    warning('LiCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.002840

% Calculate water activity
aw = 0 + 1.0038033732e+00 + -8.1163629248e-01*mf1 + -1.0230918052e+00*mf2 + 3.7441990355e+00*mf1.^2 + 1.7338650791e+00*mf1.*mf2 + -2.0112009404e+00*mf2.^2 + -3.0346356867e+02*mf1.^3 + 1.2035511291e+03*mf1.^2.*mf2 + -1.5618278318e+03*mf1.*mf2.^2 + 5.7247928902e+02*mf2.^3;

end
