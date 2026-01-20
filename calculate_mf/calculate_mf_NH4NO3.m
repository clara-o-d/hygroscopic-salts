function mf = calculate_mf_NH4NO3(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Ammonium Nitrate (NH4NO3) as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from ANPredictiveModel.m
% Valid RH range: 0.118 to 0.732 (11.8% to 73.2%)
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
if RH < 0.118
    error("below minimum relative humidity (11.8%)")
end
if RH > 0.732
    error("above maximum relative humidity (73.2%)")
end

% 4th order polynomial fit coefficients (RH as function of mass fraction)
% Fitted from data in fit_polynomial_endothermic_and_sulfates.m
A_4 = -3.2882786275;
A_3 = 9.5395623029;
A_2 = -11.3332115718;
A_1 = 5.1665983820;
A_0 = 0.0321027637;

% Invert: find mf such that RH = A_0 + A_1*mf + A_2*mf^2 + A_3*mf^3 + A_4*mf^4
f = @(xi) RH - (A_0 + A_1.*xi + A_2.*xi.^2 + A_3.*xi.^3 + A_4.*xi.^4);
mf = robust_fzero(f, 0.5456, 1.0, 0.7);
end
