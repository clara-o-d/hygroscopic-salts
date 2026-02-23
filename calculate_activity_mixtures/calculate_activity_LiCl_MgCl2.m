function aw = calculate_activity_LiCl_MgCl2(mf1, mf2)
% Calculate water activity for LiCl + MgCl2
% Component 1: LiCl
% Component 2: MgCl2
%
% Inputs:
%   mf1 - mass fraction of LiCl
%   mf2 - mass fraction of MgCl2
% Output:
%   aw - water activity
%
% Data range:
%   LiCl: 0.0080 to 0.2762
%   MgCl2: 0.0047 to 0.1860
%   aw: 0.4610 to 0.9810

% Validate inputs
if mf1 < 0.007990 || mf1 > 0.276173
    warning('LiCl mass fraction outside calibrated range');
end
if mf2 < 0.004738 || mf2 > 0.186003
    warning('MgCl2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000697

% Calculate water activity
aw = 0 + 9.9535119708e-01 + -5.2737261468e-01*mf1 + -2.8590585386e-01*mf2 + -5.7568901473e+00*mf1.^2 + -8.6476989021e+00*mf1.*mf2 + -3.7307170295e+00*mf2.^2 + -8.3366300653e+00*mf1.^3 + 1.0403223784e+02*mf1.^2.*mf2 + -9.4819177082e+01*mf1.*mf2.^2 + 2.7533767977e+01*mf2.^3;

end
