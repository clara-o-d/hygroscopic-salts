function mf = calculate_mf_ZnSO4_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.2244, 0.8127] -> RH [0.9130, 0.9962]
if RH < 0.9130
    error("Input RH (%.4f) is below the lower fit limit for ZnSO4 (0.9130)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -4.09;
A_3 = 7.274;
A_2 = -4.723;
A_1 = 1.281;
A_0 = 0.874;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.5);
end