function mf = calculate_mf_LiClO4(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.1427, 0.8440] -> RH [0.7775, 0.9931]
if RH < 0.7775
    error("Input RH (%.4f) is below the lower fit limit for LiClO4 (0.7775)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -3.381979;
A_3 = 5.227194;
A_2 = -2.982794;
A_1 = 0.631381;
A_0 = 0.947961;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1427, 0.8440, 0.49);
end
