function mf = calculate_mf_NH4Cl_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Ammonium Chloride as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = 0.4475; 
A_3 = -1.067; 
A_2 = 0.08897;
A_1 = -0.2293; 
A_0 = 1;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0308, 0.4880, 0.25);
if RH < 0.76   % mf > 0.4880
    error("below deliquescence relative humidity")
end 
end