function mf = calculate_mf_NiSO4_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.2100, 0.7687] -> RH [0.9390, 0.9962]
if RH < 0.9390
    error("Input RH (%.4f) is below the lower fit limit for NiSO4 (0.9390)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -2.827;
A_3 = 4.673;
A_2 = -2.835;
A_1 = 0.71;
A_0 = 0.934;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.5);
end