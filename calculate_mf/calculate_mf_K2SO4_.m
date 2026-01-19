function mf = calculate_mf_K2SO4_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.1442, 0.5742] -> RH [0.9720, 0.9958]
if RH < 0.9720
    error("Input RH (%.4f) is below the lower fit limit for K2SO4 (0.9720)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -0.03918;
A_3 = -0.05269;
A_2 = 0.01542;
A_1 = -0.03377;
A_0 = 1.001;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.35);
end