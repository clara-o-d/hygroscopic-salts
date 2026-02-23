function aw = calculate_activity_LiCl_BaCl2(mf1, mf2)
% Calculate water activity for LiCl + BaCl2
% Component 1: LiCl
% Component 2: BaCl2
%
% Inputs:
%   mf1 - mass fraction of LiCl
%   mf2 - mass fraction of BaCl2
% Output:
%   aw - water activity
%
% Data range:
%   LiCl: 0.0051 to 0.1620
%   BaCl2: 0.0103 to 0.2104
%   aw: 0.7084 to 0.9886

% Validate inputs
if mf1 < 0.005062 || mf1 > 0.161999
    warning('LiCl mass fraction outside calibrated range');
end
if mf2 < 0.010304 || mf2 > 0.210444
    warning('BaCl2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.002126

% Calculate water activity
aw = 0 + 9.9901122730e-01 + -6.4036738667e-01*mf1 + -4.2466691153e-01*mf2 + 1.4767275225e+01*mf1.^2 + -5.4386362138e+01*mf1.*mf2 + 1.2588020702e+01*mf2.^2 + 3.6632995001e+02*mf1.^3 + -1.1676221896e+03*mf1.^2.*mf2 + 9.1614291505e+02*mf1.*mf2.^2 + -1.4616160980e+02*mf2.^3;

end
