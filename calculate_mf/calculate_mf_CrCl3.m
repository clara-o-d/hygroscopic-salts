function mf = calculate_mf_CrCl3(RH)
% calculate_mf_CrCl3
% Calculates mass fraction of CrCl3 based on Relative Humidity.
% Range: 0.8890 < RH < 0.9956
% Data source: sources_list_.xlsx - Chromium chloride data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.8890 || RH > 0.9956
    error("RH %.4f is outside valid range [%.4f, %.4f] for CrCl3", RH, 0.8890, 0.9956);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = 44.097300085036338;
A_3 = -18.289170404594756;
A_2 = 0.054411660546782;
A_1 = -0.283216723233178;
A_0 = 1.000088607099477;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0156, 0.1815, 0.0985);

end
