function mf = calculate_mf_SrI2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.8626, 0.9641] -> RH [0.6785, 0.9569]
if RH < 0.6785
    error("Input RH (%.4f) is below the lower fit limit for SrI2 (0.6785)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -11202.092662;
A_3 = 40360.069556;
A_2 = -54527.467244;
A_1 = 32738.541412;
A_0 = -7369.305883;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.8626, 0.9641, 0.91);
end
