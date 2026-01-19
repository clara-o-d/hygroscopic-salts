function mf = calculate_mf_KCl_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Potassium Chloride as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = -1.306; 
A_3 = 0.9188; 
A_2 = -0.4134;
A_1 = -0.0655; 
A_0 = 0.9984;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0581, 0.5813, 0.32);
if RH < 0.83
    error("below deliquescence relative humidity")
end 
end