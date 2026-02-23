function aw = calculate_activity_NaC2H3O2_LiC2H3O2(mf1, mf2)
% Calculate water activity for NaC2H3O2 + LiC2H3O2
% Component 1: NaC2H3O2
% Component 2: LiC2H3O2
%
% Inputs:
%   mf1 - mass fraction of NaC2H3O2
%   mf2 - mass fraction of LiC2H3O2
% Output:
%   aw - water activity
%
% Data range:
%   NaC2H3O2: 0.0184 to 0.2437
%   LiC2H3O2: 0.0148 to 0.2244
%   aw: 0.7552 to 0.9846

% Validate inputs
if mf1 < 0.018404 || mf1 > 0.243676
    warning('NaC2H3O2 mass fraction outside calibrated range');
end
if mf2 < 0.014785 || mf2 > 0.224370
    warning('LiC2H3O2 mass fraction outside calibrated range');
end

% Bivariate polynomial fit
% Degree: 3, RMSE: 0.000189

% Calculate water activity
aw = 0 + 9.9973523386e-01 + -3.9645448659e-01*mf1 + -4.6775065503e-01*mf2 + -7.3637291786e-01*mf1.^2 + -1.6116330289e+00*mf1.*mf2 + -1.0579693465e+00*mf2.^2 + -1.4068193996e+00*mf1.^3 + 2.6994181989e+00*mf1.^2.*mf2 + 3.4744216935e+00*mf1.*mf2.^2 + -7.4140780633e-01*mf2.^3;

end
