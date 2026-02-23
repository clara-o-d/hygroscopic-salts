function aw = calculate_activity_HClO4_LiClO4(mf1, mf2)
% Calculate water activity for HClO4 + LiClO4
% Component 1: HClO4
% Component 2: LiClO4
%
% Inputs:
%   mf1 - mass fraction of HClO4
%   mf2 - mass fraction of LiClO4
% Output:
%   aw - water activity
%
% Data range:
%   HClO4: 0.0405 to 0.2497
%   LiClO4: 0.0430 to 0.2614
%   aw: 0.7610 to 0.9317

% Validate inputs
if mf1 < 0.040508 || mf1 > 0.249703
    warning('HClO4 mass fraction outside calibrated range');
end
if mf2 < 0.043003 || mf2 > 0.261420
    warning('LiClO4 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000056

% Calculate water activity
aw = 0 + 9.9817237603e-01 + -3.1398468751e-01*mf1 + -2.7783933351e-01*mf2 + -4.7783569608e-01*mf1.^2 + -1.3512219571e+00*mf1.*mf2 + -8.0456442174e-01*mf2.^2 + -5.1205168883e+00*mf1.^3 + 6.6297953016e+00*mf1.^2.*mf2 + -6.7476444452e+00*mf1.*mf2.^2 + -1.4669331214e-01*mf2.^3;

end
