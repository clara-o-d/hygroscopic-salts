function mf = calculate_mf_NaClO3(RH)
% calculate_mf_NaClO3
% Calculates mass fraction of NaClO3 based on Relative Humidity.
% Range: 0.8457 < RH < 0.9950
% Data source: Paper 25 in sources_list_.xlsx - Sodium chlorate data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.8457 || RH > 0.9950
    error("RH %.4f is outside valid range [%.4f, %.4f] for NaClO3", RH, 0.8457, 0.9950);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = 0.877804743405;
A_3 = -1.351933751306;
A_2 = -0.094123911525;
A_1 = -0.460194467223;
A_0 = 0.999841300531;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0105, 0.2714, 0.1410);

end
