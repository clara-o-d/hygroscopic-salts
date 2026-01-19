function mf = calculate_mf_RbCl(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.5479, 0.8494] -> RH [0.7423, 0.9527]
if RH < 0.7423
    error("Input RH (%.4f) is below the lower fit limit for RbCl (0.7423)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -235.724643;
A_3 = 660.252219;
A_2 = -690.848892;
A_1 = 319.307884;
A_0 = -53.960292;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.5479, 0.8494, 0.70);
end
