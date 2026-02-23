function aw = calculate_activity_NaCl_CsCl(mf1, mf2)
% Calculate water activity for NaCl + CsCl
% Component 1: NaCl
% Component 2: CsCl
%
% Inputs:
%   mf1 - mass fraction of NaCl
%   mf2 - mass fraction of CsCl
% Output:
%   aw - water activity
%
% Data range:
%   NaCl: 0.0058 to 0.2375
%   CsCl: 0.0326 to 0.4730
%   aw: 0.7020 to 0.9900

% Validate inputs
if mf1 < 0.005810 || mf1 > 0.237506
    warning('NaCl mass fraction outside calibrated range');
end
if mf2 < 0.032575 || mf2 > 0.472952
    warning('CsCl mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.006851

% Calculate water activity
aw = 0 + 9.8819774962e-01 + -1.0190968844e+00*mf1 + 3.7639054396e-01*mf2 + -2.5027723913e+01*mf1.^2 + 2.7313754515e+01*mf1.*mf2 + -9.0982413347e+00*mf2.^2 + -2.9357514974e+02*mf1.^3 + 5.0507781541e+02*mf1.^2.*mf2 + -2.7139986611e+02*mf1.*mf2.^2 + 4.7273886882e+01*mf2.^3;

end
