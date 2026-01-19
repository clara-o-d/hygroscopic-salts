function mf = calculate_mf_BaCl2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.5481, 0.7696] -> RH [0.9375, 0.9731]
if RH < 0.9375
    error("Input RH (%.4f) is below the lower fit limit for BaCl2 (0.9375)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = 51.737901;
A_3 = -141.918321;
A_2 = 144.292411;
A_1 = -64.600718;
A_0 = 11.731931;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.5481, 0.7696, 0.66);
end
