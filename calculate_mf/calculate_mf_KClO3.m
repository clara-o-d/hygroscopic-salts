function mf = calculate_mf_KClO3(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.1429, 0.3685] -> RH [0.9800, 0.9936]
if RH < 0.9800
    error("Input RH (%.4f) is below the lower fit limit for KClO3 (0.9800)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = 0.009347;
A_3 = -0.050959;
A_2 = -0.018595;
A_1 = -0.041055;
A_0 = 1.000016;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1429, 0.3685, 0.26);
end
