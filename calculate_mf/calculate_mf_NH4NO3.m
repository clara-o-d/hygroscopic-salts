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
if RH < 0.118 || RH > 0.732 
    error("below deliquescence relative humidity or above range") 
end  
A_4 = -3.288; 
A_3 = 9.54; 
A_2 = -11.33;
A_1 = 5.167; 
A_0 = 0.0321;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.5456, 1.0, 0.7728);

end
% ---------------------------------------------------------
