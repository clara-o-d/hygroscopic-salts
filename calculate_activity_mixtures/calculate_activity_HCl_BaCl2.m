function aw = calculate_activity_HCl_BaCl2(mf1, mf2)
% Calculate water activity for HCl + BaCl2
% Component 1: HCl
% Component 2: BaCl2
%
% Inputs:
%   mf1 - mass fraction of HCl
%   mf2 - mass fraction of BaCl2
% Output:
%   aw - water activity
%
% Data range:
%   HCl: 0.0129 to 0.0650
%   BaCl2: 0.0189 to 0.2071
%   aw: 0.9150 to 0.9750

% Validate inputs
if mf1 < 0.012902 || mf1 > 0.065009
    warning('HCl mass fraction outside calibrated range');
end
if mf2 < 0.018917 || mf2 > 0.207054
    warning('BaCl2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000640

% Calculate water activity
aw = 0 + 1.1886741624e+01 + -7.4865917666e+02*mf1 + -6.8933244862e+01*mf2 + 7.0379864648e+03*mf1.^2 + 3.8406482014e+03*mf1.*mf2 + -1.0638191031e+02*mf2.^2 + 1.4910715301e+04*mf1.^3 + 4.4312258142e+03*mf1.^2.*mf2 + 5.8911206655e+03*mf1.*mf2.^2 + 6.7380798242e+00*mf2.^3;

end
