function mf = calculate_mf_SrI2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of SrI2 as a
% function of the Relative Humidity at a temperature of 25C
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
if RH < 0.6785 || RH > 0.9569 
    error("below deliquescence relative humidity or above range") 
end  
A_4 = -17.11; 
A_3 = 23.0; 
A_2 = -12.67;
A_1 = 2.873; 
A_0 = 0.7376;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.2488, 0.5866, 0.4177);

end
% ---------------------------------------------------------
