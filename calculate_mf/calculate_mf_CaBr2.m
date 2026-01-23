function mf = calculate_mf_CaBr2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of CaBr2 as a
% function of the Relative Humidity at a temperature of 25C
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
if RH < 0.6395 || RH > 0.954 
    error("below deliquescence relative humidity or above range") 
end  
A_4 = 7.902; 
A_3 = -12.19; 
A_2 = 4.102;
A_1 = -0.8593; 
A_0 = 1.034;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1674, 0.4788, 0.3231);

end
% ---------------------------------------------------------
