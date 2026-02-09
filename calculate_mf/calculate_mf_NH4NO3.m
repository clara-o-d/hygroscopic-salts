function mf = calculate_mf_NH4NO3(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of NH4NO3 as a
% function of the Relative Humidity at a temperature of 25C
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
if RH < 0.62 || RH > 0.95 
    error("below deliquescence relative humidity or above range") 
end  
A_4 = -4.648; 
A_3 = 6.018; 
A_2 = -3.084;
A_1 = 0.2392; 
A_0 = 0.9565;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1268, 0.6589, 0.4);

end
% ---------------------------------------------------------
