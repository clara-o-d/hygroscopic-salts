function aw = calculate_activity_HCl_BaClO42(mf1, mf2)
% Calculate water activity for HCl + Ba(ClO4)2
% Component 1: HCl
% Component 2: Ba(ClO4)2
%
% Inputs:
%   mf1 - mass fraction of HCl
%   mf2 - mass fraction of Ba(ClO4)2
% Output:
%   aw - water activity
%
% Data range:
%   HCl: 0.0172 to 0.0506
%   Ba(ClO4)2: 0.0658 to 0.5049
%   aw: 0.7200 to 0.9380

% Validate inputs
if mf1 < 0.017175 || mf1 > 0.050640
    warning('HCl mass fraction outside calibrated range');
end
if mf2 < 0.065832 || mf2 > 0.504891
    warning('Ba(ClO4)2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.003412

% Calculate water activity
aw = 0 + 9.2979978017e-01 + -1.1102396060e+00*mf1 + 3.5638872330e-01*mf2 + 1.6592607708e+01*mf1.^2 + 8.2326969724e+00*mf1.*mf2 + -1.6943987743e+00*mf2.^2 + 1.7182100008e+00*mf1.^3 + -2.3396055950e+02*mf1.^2.*mf2 + 2.2349282021e+00*mf1.*mf2.^2 + 3.5266573660e-01*mf2.^3;

end
