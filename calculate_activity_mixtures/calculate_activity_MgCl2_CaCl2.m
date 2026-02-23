function aw = calculate_activity_MgCl2_CaCl2(mf1, mf2)
% Calculate water activity for MgCl2 + CaCl2
% Component 1: MgCl2
% Component 2: CaCl2
%
% Inputs:
%   mf1 - mass fraction of MgCl2
%   mf2 - mass fraction of CaCl2
% Output:
%   aw - water activity
%
% Data range:
%   MgCl2: 0.0094 to 0.2222
%   CaCl2: 0.0110 to 0.2498
%   aw: 0.6030 to 0.9900

% Validate inputs
if mf1 < 0.009431 || mf1 > 0.222171
    warning('MgCl2 mass fraction outside calibrated range');
end
if mf2 < 0.010976 || mf2 > 0.249779
    warning('CaCl2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.015789

% Calculate water activity
aw = 0 + 9.6939297855e-01 + 1.6601387887e+00*mf1 + -2.6508680893e-01*mf2 + -3.9827646186e+01*mf1.^2 + 2.1113353607e+00*mf1.*mf2 + -3.9583281122e+00*mf2.^2 + 1.1346889451e+03*mf1.^3 + -2.9798932706e+03*mf1.^2.*mf2 + 2.6233729200e+03*mf1.*mf2.^2 + -6.8723113060e+02*mf2.^3;

end
