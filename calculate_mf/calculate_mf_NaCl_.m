function mf = calculate_mf_NaCl_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Sodium Chloride as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = 5.863; 
A_3 = -5.545; 
A_2 = -0.332;
A_1 = -0.5597; 
A_0 = 0.9998;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0116, 0.2596, 0.13);
if RH < 0.73     % if mf > 0.2596
    error("below deliquescence relative humidity")
end 
end