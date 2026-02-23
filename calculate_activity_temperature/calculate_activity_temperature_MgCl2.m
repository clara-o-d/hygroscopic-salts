function aw = calculate_activity_temperature_MgCl2(mf, T)
% Calculate water activity for MgCl2 as a function of mass fraction and temperature
%
% Inputs:
%   mf - mass fraction of MgCl2
%   T  - temperature (°C)
% Output:
%   aw - water activity
%
% Data range:
%   Mass fraction: 0.0460 to 0.2268
%   Temperature: 25.0°C to 139.9°C
%   Water activity: 0.8335 to 0.9828
%
% Fit quality:
%   Polynomial degree: 4
%   RMSE: 0.000112
%   Number of data points: 49

% Validate inputs
if mf < 0.045957 || mf > 0.226786
    warning('MgCl2 mass fraction outside calibrated range (%.4f to %.4f)', 0.045957, 0.226786);
end
if T < 25.00 || T > 139.85
    warning('Temperature outside calibrated range (%.1f°C to %.1f°C)', 25.00, 139.85);
end

% Bivariate polynomial fit: aw = f(mf, T)
% Polynomial degree: 4, RMSE: 0.000112

% Calculate water activity
aw = 0 + 6.001722862689e-06 + -8.535045834631e-04*mf + 1.377179682261e-04*T + -1.641030040384e+00*mf.^2 + -1.509901969766e-02*mf.*T + 2.443976359762e-03*T.^2 + -5.562706110888e+00*mf.^3 + 9.473362452980e-03*mf.^2.*T + 1.877304191869e-04*mf.*T.^2 + -3.774222717514e-05*T.^3 + -1.274650791665e+00*mf.^4 + 9.852253500359e-03*mf.^3.*T + -2.825263586951e-05*mf.^2.*T.^2 + -6.600797824553e-07*mf.*T.^3 + 1.474760782590e-07*T.^4;

% Ensure water activity is between 0 and 1
aw = max(0, min(1, aw));

end
