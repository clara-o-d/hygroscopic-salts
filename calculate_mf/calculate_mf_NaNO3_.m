function mf = calculate_mf_NaNO3_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Sodium Nitrate as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = 0.009245; 
A_3 = -0.1583; 
A_2 = -0.02788;
A_1 = -0.08334; 
A_0 = 0.9999;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0041, 0.2866, 0.14);
if RH < 0.71    % mf > 0.2866
    error("below deliquescence relative humidity")
end 
end