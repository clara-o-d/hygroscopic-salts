function aw = calculate_activity_LiCl_CaCl2(mf1, mf2)
% Calculate water activity for LiCl + CaCl2
% Component 1: LiCl
% Component 2: CaCl2
%
% Inputs:
%   mf1 - mass fraction of LiCl
%   mf2 - mass fraction of CaCl2
% Output:
%   aw - water activity
%
% Data range:
%   LiCl: 0.0038 to 0.3851
%   CaCl2: 0.0011 to 0.3879
%   aw: 0.2490 to 0.9930

% Validate inputs
if mf1 < 0.003801 || mf1 > 0.385054
    warning('LiCl mass fraction outside calibrated range');
end
if mf2 < 0.001109 || mf2 > 0.387891
    warning('CaCl2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.004293

% Calculate water activity
aw = 0 + 9.9600758904e-01 + -4.3798308761e-01*mf1 + -2.7531589608e-01*mf2 + -8.0824528644e+00*mf1.^2 + -5.2936555490e+00*mf1.*mf2 + -2.6770297887e+00*mf2.^2 + 2.1836004826e+01*mf1.^3 + -3.0283407983e+01*mf1.^2.*mf2 + 5.6582998348e+01*mf1.*mf2.^2 + -1.3297586288e+01*mf2.^3;

end
