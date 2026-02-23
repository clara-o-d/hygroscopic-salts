function aw = calculate_activity_LiCl_CsCl(mf1, mf2)
% Calculate water activity for LiCl + CsCl
% Component 1: LiCl
% Component 2: CsCl
%
% Inputs:
%   mf1 - mass fraction of LiCl
%   mf2 - mass fraction of CsCl
% Output:
%   aw - water activity
%
% Data range:
%   LiCl: 0.0042 to 0.2028
%   CsCl: 0.0326 to 0.5025
%   aw: 0.5560 to 0.9900

% Validate inputs
if mf1 < 0.004222 || mf1 > 0.202783
    warning('LiCl mass fraction outside calibrated range');
end
if mf2 < 0.032575 || mf2 > 0.502527
    warning('CsCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.009326

% Calculate water activity
aw = 0 + 9.8240979993e-01 + -1.2603386150e+00*mf1 + 4.5540053348e-01*mf2 + -4.4633445675e+01*mf1.^2 + 3.2921770808e+01*mf1.*mf2 + -8.8023674885e+00*mf2.^2 + -4.3695914134e+02*mf1.^3 + 6.2030147589e+02*mf1.^2.*mf2 + -2.6527515707e+02*mf1.*mf2.^2 + 3.7070216565e+01*mf2.^3;

end
