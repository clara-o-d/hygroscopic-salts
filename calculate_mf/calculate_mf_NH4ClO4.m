function mf = calculate_mf_NH4ClO4(RH)
% calculate_mf_NH4ClO4
% Calculates mass fraction of NH4ClO4 based on Relative Humidity.
% Range: 0.9444 < RH < 0.9968

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.9444 || RH > 0.9968
    error("RH %.4f is outside valid range [%.4f, %.4f] for NH4ClO4", RH, 0.9444, 0.9968);
end

% Coefficients
A_4 = 7.116;
A_3 = -2.743;
A_2 = 0.2342;
A_1 = -0.2754;
A_0 = 0.9999;

% Polynomial function: RH = A0 + A1*x + ... + A4*x^4
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0116, 0.1979, 0.1047);

end
