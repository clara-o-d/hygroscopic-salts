function aw = calculate_activity_LiCl_KCl(mf1, mf2)
% Calculate water activity for LiCl + KCl
% Component 1: LiCl
% Component 2: KCl
%
% Inputs:
%   mf1 - mass fraction of LiCl
%   mf2 - mass fraction of KCl
% Output:
%   aw - water activity
%
% Data range:
%   LiCl: 0.0042 to 0.1450
%   KCl: 0.0074 to 0.2297
%   aw: 0.7310 to 0.9900

% Validate inputs
if mf1 < 0.004222 || mf1 > 0.144989
    warning('LiCl mass fraction outside calibrated range');
end
if mf2 < 0.007400 || mf2 > 0.229703
    warning('KCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.002768

% Calculate water activity
aw = 0 + 1.0008255643e+00 + -7.5139780685e-01*mf1 + -4.5145263689e-01*mf2 + 5.0067807561e+00*mf1.^2 + -1.2715386616e+01*mf1.*mf2 + 3.6767238683e+00*mf2.^2 + 3.1596358949e+02*mf1.^3 + -7.1557935049e+02*mf1.^2.*mf2 + 4.7894724318e+02*mf1.*mf2.^2 + -9.7498345050e+01*mf2.^3;

end
