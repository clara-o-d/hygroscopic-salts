function aw = calculate_activity_NaCl_LiCl(mf1, mf2)
% Calculate water activity for NaCl + LiCl
% Component 1: NaCl
% Component 2: LiCl
%
% Inputs:
%   mf1 - mass fraction of NaCl
%   mf2 - mass fraction of LiCl
% Output:
%   aw - water activity
%
% Data range:
%   NaCl: 0.0057 to 0.2366
%   LiCl: 0.0047 to 0.1549
%   aw: 0.7604 to 0.9826

% Validate inputs
if mf1 < 0.005745 || mf1 > 0.236635
    warning('NaCl mass fraction outside calibrated range');
end
if mf2 < 0.004704 || mf2 > 0.154883
    warning('LiCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000135

% Calculate water activity
aw = 0 + 9.9923363828e-01 + -5.2943359056e-01*mf1 + -7.2843309057e-01*mf2 + -8.4028561842e-01*mf1.^2 + -2.6453389541e+00*mf1.*mf2 + -2.9350244352e+00*mf2.^2 + -2.5554851709e+00*mf1.^3 + -6.9702282620e-01*mf1.^2.*mf2 + -1.0632207011e-01*mf1.*mf2.^2 + -6.4191257893e+00*mf2.^3;

end
