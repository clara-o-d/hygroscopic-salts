function mf = calculate_mf_NaI(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% Fit on data range: MassFrac [0.5560, 0.9128] -> RH [0.5801, 0.9669]
if RH < 0.5801
    error("Input RH (%.4f) is below the lower fit limit for NaI (0.5801)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -111.341752;
A_3 = 309.657672;
A_2 = -322.717222;
A_1 = 148.853016;
A_0 = -24.615377;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.5560, 0.9128, 0.73);
end
