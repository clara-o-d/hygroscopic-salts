function mf = calculate_mf_Li2SO4_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.0629, 0.6681] -> RH [0.8530, 0.9956]
if RH < 0.8530
    error("Input RH (%.4f) is below the lower fit limit for Li2SO4 (0.8530)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -1.861;
A_3 = 1.79;
A_2 = -0.7378;
A_1 = 0.03636;
A_0 = 0.9955;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.35);
end