function aw = calculate_activity_NaNO3_LiNO3(mf1, mf2)
% Calculate water activity for NaNO3 + LiNO3
% Component 1: NaNO3
% Component 2: LiNO3
%
% Inputs:
%   mf1 - mass fraction of NaNO3
%   mf2 - mass fraction of LiNO3
% Output:
%   aw - water activity
%
% Data range:
%   NaNO3: 0.0176 to 0.3599
%   LiNO3: 0.0178 to 0.2471
%   aw: 0.7627 to 0.9841

% Validate inputs
if mf1 < 0.017575 || mf1 > 0.359893
    warning('NaNO3 mass fraction outside calibrated range');
end
if mf2 < 0.017757 || mf2 > 0.247139
    warning('LiNO3 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000863

% Calculate water activity
aw = 0 + 1.0009023788e+00 + -4.1915240838e-01*mf1 + -4.7424287880e-01*mf2 + 1.8367426574e-01*mf1.^2 + -1.1198738472e-01*mf1.*mf2 + -1.1118592163e+00*mf2.^2 + -1.0758124055e+00*mf1.^3 + 3.6527023272e-01*mf1.^2.*mf2 + 2.8536041670e-01*mf1.*mf2.^2 + -1.0663149694e+00*mf2.^3;

end
