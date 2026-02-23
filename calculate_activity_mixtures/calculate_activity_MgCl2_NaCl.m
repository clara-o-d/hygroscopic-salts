function aw = calculate_activity_MgCl2_NaCl(mf1, mf2)
% Calculate water activity for MgCl2 + NaCl
% Component 1: MgCl2
% Component 2: NaCl
%
% Inputs:
%   mf1 - mass fraction of MgCl2
%   mf2 - mass fraction of NaCl
% Output:
%   aw - water activity
%
% Data range:
%   MgCl2: 0.0094 to 0.1033
%   NaCl: 0.0087 to 0.1746
%   aw: 0.7993 to 0.9904

% Validate inputs
if mf1 < 0.009431 || mf1 > 0.103303
    warning('MgCl2 mass fraction outside calibrated range');
end
if mf2 < 0.008690 || mf2 > 0.174613
    warning('NaCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.002177

% Calculate water activity
aw = 0 + 1.0001142213e+00 + -4.5555102614e-01*mf1 + -6.1611360611e-01*mf2 + 7.1025962416e+00*mf1.^2 + -1.3544046734e+01*mf1.*mf2 + 3.7751820019e+00*mf2.^2 + 1.2988371386e+03*mf1.^3 + -2.6699910320e+03*mf1.^2.*mf2 + 1.5937208772e+03*mf1.*mf2.^2 + -2.8870812705e+02*mf2.^3;

end
