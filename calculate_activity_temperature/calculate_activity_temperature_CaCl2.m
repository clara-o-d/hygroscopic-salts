function aw = calculate_activity_temperature_CaCl2(mf, T)
% Calculate water activity for CaCl2 as a function of mass fraction and temperature
%
% Inputs:
%   mf - mass fraction of CaCl2
%   T  - temperature (°C)
% Output:
%   aw - water activity
%
% Data range:
%   Mass fraction: 0.0017 to 0.5505
%   Temperature: 25.0°C to 90.0°C
%   Water activity: 0.2868 to 0.9995
%
% Fit quality:
%   Polynomial degree: 4
%   RMSE: 0.003498
%   Number of data points: 74

% Validate inputs
if mf < 0.001725 || mf > 0.550458
    warning('CaCl2 mass fraction outside calibrated range (%.4f to %.4f)', 0.001725, 0.550458);
end
if T < 25.00 || T > 90.00
    warning('Temperature outside calibrated range (%.1f°C to %.1f°C)', 25.00, 90.00);
end

% Bivariate polynomial fit: aw = f(mf, T)
% Polynomial degree: 4, RMSE: 0.003498

% Calculate water activity
aw = 0 + 1.282814544321e-05 + -1.682900107898e-03*mf + 2.540532176048e-04*T + 1.136498423826e+00*mf.^2 + -2.365655086021e-02*mf.*T + 3.570648737787e-03*T.^2 + -1.610859893139e+01*mf.^3 + 1.461050527949e-02*mf.^2.*T + 4.476474757904e-04*mf.*T.^2 + -9.498565186915e-05*T.^3 + 1.911375099507e+01*mf.^4 + -7.860188957046e-03*mf.^3.*T + -1.664861731100e-05*mf.^2.*T.^2 + -2.680867780484e-06*mf.*T.^3 + 6.295105819828e-07*T.^4;

% Ensure water activity is between 0 and 1
aw = max(0, min(1, aw));

end
