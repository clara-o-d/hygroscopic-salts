function aw = calculate_activity_HClO4_NaClO4(mf1, mf2)
% Calculate water activity for HClO4 + NaClO4
% Component 1: HClO4
% Component 2: NaClO4
%
% Inputs:
%   mf1 - mass fraction of HClO4
%   mf2 - mass fraction of NaClO4
% Output:
%   aw - water activity
%
% Data range:
%   HClO4: 0.0686 to 0.4191
%   NaClO4: 0.0718 to 0.5598
%   aw: 0.4444 to 0.8942

% Validate inputs
if mf1 < 0.068551 || mf1 > 0.419129
    warning('HClO4 mass fraction outside calibrated range');
end
if mf2 < 0.071760 || mf2 > 0.559805
    warning('NaClO4 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000342

% Calculate water activity
aw = 0 + 9.6003201775e-01 + 7.9528817501e-02*mf1 + 1.3441209926e-02*mf2 + -1.8093388939e+00*mf1.^2 + -2.5195442622e+00*mf1.*mf2 + -7.9957697904e-01*mf2.^2 + -2.5840074359e+00*mf1.^3 + 6.0797513464e+00*mf1.^2.*mf2 + -4.3703411340e-01*mf1.*mf2.^2 + 2.5633380117e-01*mf2.^3;

end
