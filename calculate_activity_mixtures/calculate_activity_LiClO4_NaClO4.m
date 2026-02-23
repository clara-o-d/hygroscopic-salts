function aw = calculate_activity_LiClO4_NaClO4(mf1, mf2)
% Calculate water activity for LiClO4 + NaClO4
% Component 1: LiClO4
% Component 2: NaClO4
%
% Inputs:
%   mf1 - mass fraction of LiClO4
%   mf2 - mass fraction of NaClO4
% Output:
%   aw - water activity
%
% Data range:
%   LiClO4: 0.0209 to 0.2723
%   NaClO4: 0.0232 to 0.3502
%   aw: 0.7715 to 0.9723

% Validate inputs
if mf1 < 0.020940 || mf1 > 0.272325
    warning('LiClO4 mass fraction outside calibrated range');
end
if mf2 < 0.023213 || mf2 > 0.350219
    warning('NaClO4 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000056

% Calculate water activity
aw = 0 + 9.9946365151e-01 + -3.0660901590e-01*mf1 + -2.5653248859e-01*mf2 + -6.2146831342e-01*mf1.^2 + -5.8623849723e-01*mf1.*mf2 + -2.2433077868e-01*mf2.^2 + -3.2101083063e+00*mf1.^3 + 3.8682383419e+00*mf1.^2.*mf2 + -3.1650876008e+00*mf1.*mf2.^2 + 4.8725831456e-02*mf2.^3;

end
