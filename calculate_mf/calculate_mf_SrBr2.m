function mf = calculate_mf_SrBr2(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.7311, 0.9190] -> RH [0.7776, 0.9571]
if RH < 0.7776
    error("Input RH (%.4f) is below the lower fit limit for SrBr2 (0.7776)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -371.066326;
A_3 = 1170.024726;
A_2 = -1385.238827;
A_1 = 729.569711;
A_0 = -143.218582;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.7311, 0.9190, 0.83);
end
