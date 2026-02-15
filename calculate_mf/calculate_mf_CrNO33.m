function mf = calculate_mf_CrNO33(RH)
% calculate_mf_CrNO33
% Calculates mass fraction of Cr(NO3)3 based on Relative Humidity.
% Range: 0.8838 < RH < 0.9957
% Data source: sources_list_.xlsx - Chromium nitrate data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.8838 || RH > 0.9957
    error("RH %.4f is outside valid range [%.4f, %.4f] for Cr(NO3)3", RH, 0.8838, 0.9957);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = 1.321591844547575;
A_3 = -3.279370134297320;
A_2 = -0.224316552945786;
A_1 = -0.180573115077002;
A_0 = 1.000121180415139;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0232, 0.2631, 0.1432);

end
