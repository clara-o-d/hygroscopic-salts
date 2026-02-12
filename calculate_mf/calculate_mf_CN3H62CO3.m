function mf = calculate_mf_CN3H62CO3(RH)
% calculate_mf_CN3H62CO3
% Calculates mass fraction of (CN3H6)2CO3 based on Relative Humidity.
% Range: 0.9439 < RH < 0.9999

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.9439 || RH > 0.9999
    error("RH %.4f is outside valid range [%.4f, %.4f] for (CN3H6)2CO3", RH, 0.9439, 0.9999);
end

% Coefficients
A_4 = 2.515;
A_3 = -2.481;
A_2 = 0.7115;
A_1 = -0.231;
A_0 = 0.9999;

% Polynomial function: RH = A0 + A1*x + ... + A4*x^4
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0002, 0.3201, 0.1601);

end
