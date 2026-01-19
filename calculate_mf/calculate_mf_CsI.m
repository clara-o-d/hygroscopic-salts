function mf = calculate_mf_CsI(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.8344, 0.9067] -> RH [0.9124, 0.9624]
if RH < 0.9124
    error("Input RH (%.4f) is below the lower fit limit for CsI (0.9124)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -21853.278178;
A_3 = 75825.598171;
A_2 = -98634.949546;
A_1 = 57009.062207;
A_0 = -12351.841453;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.8344, 0.9067, 0.87);
end
