function mf = calculate_mf_CaBr2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.6905, 0.9107] -> RH [0.6395, 0.9540]
if RH < 0.6395
    error("Input RH (%.4f) is below the lower fit limit for CaBr2 (0.6395)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -254.497771;
A_3 = 760.127187;
A_2 = -852.774526;
A_1 = 425.454914;
A_0 = -78.625399;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.6905, 0.9107, 0.80);
end
