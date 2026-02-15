function mf = calculate_mf_NaBrO3(RH)
% calculate_mf_NaBrO3
% Calculates mass fraction of NaBrO3 based on Relative Humidity.
% Range: 0.8985 < RH < 0.9951
% Data source: Paper 25 in sources_list_.xlsx - Sodium bromate data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.8985 || RH > 0.9951
    error("RH %.4f is outside valid range [%.4f, %.4f] for NaBrO3", RH, 0.8985, 0.9951);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = -0.661799371345;
A_3 = -0.331374618964;
A_2 = -0.027904331786;
A_1 = -0.324058058181;
A_0 = 0.999888949337;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0149, 0.2739, 0.1444);

end
