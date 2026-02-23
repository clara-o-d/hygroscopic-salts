function aw = calculate_activity_temperature_LiCl(mf, T)
% Calculate water activity for LiCl as a function of mass fraction and temperature
%
% Inputs:
%   mf - mass fraction of LiCl
%   T  - temperature (°C)
% Output:
%   aw - water activity
%
% Data range:
%   Mass fraction: 0.0415 to 0.4407
%   Temperature: 25.0°C to 100.0°C
%   Water activity: 0.1300 to 0.9647
%
% Fit quality:
%   Polynomial degree: 4
%   RMSE: 0.001100
%   Number of data points: 54

% Validate inputs
if mf < 0.041481 || mf > 0.440682
    warning('LiCl mass fraction outside calibrated range (%.4f to %.4f)', 0.041481, 0.440682);
end
if T < 25.00 || T > 100.00
    warning('Temperature outside calibrated range (%.1f°C to %.1f°C)', 25.00, 100.00);
end

% Bivariate polynomial fit: aw = f(mf, T)
% Polynomial degree: 4, RMSE: 0.001100

% Calculate water activity
aw = 0 + 9.883383998075e-01 + -3.734438377826e-01*mf + -2.951877030381e-07*T + -5.081274298923e+00*mf.^2 + -4.549756219824e-03*mf.*T + 3.412248725098e-06*T.^2 + -1.038920698914e+01*mf.^3 + 4.367333741891e-02*mf.^2.*T + -6.277839354456e-06*mf.*T.^2 + -2.199667256385e-08*T.^3 + 3.073699674209e+01*mf.^4 + -6.578372514926e-02*mf.^3.*T + 1.128053857359e-05*mf.^2.*T.^2 + 2.416158947514e-08*mf.*T.^3 + 2.800426152386e-11*T.^4;

% Ensure water activity is between 0 and 1
aw = max(0, min(1, aw));

end
