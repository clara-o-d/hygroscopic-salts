function mf = calculate_mf_BaNO3(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.2749, 0.6026] -> RH [0.9859, 0.9958]
if RH < 0.9859
    error("Input RH (%.4f) is below the lower fit limit for Ba(NO3)2 (0.9859)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -0.239984;
A_3 = 0.363612;
A_2 = -0.241044;
A_1 = 0.053942;
A_0 = 0.993004;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.2749, 0.6026, 0.44);
end
