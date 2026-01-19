function mf = calculate_mf_CaI2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.8451, 0.9332] -> RH [0.8321, 0.9524]
if RH < 0.8321
    error("Input RH (%.4f) is below the lower fit limit for CaI2 (0.8321)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -14038.400517;
A_3 = 49710.013994;
A_2 = -65994.841633;
A_1 = 38930.542523;
A_0 = -8608.760658;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.8451, 0.9332, 0.89);
end
