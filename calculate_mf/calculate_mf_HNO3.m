function mf = calculate_mf_HNO3(RH)
% calculate_mf_HNO3
% Calculates mass fraction of HNO3 based on Relative Humidity.
% Range: 0.8294 < RH < 0.9949
% Data source: Paper 25 in sources_list_.xlsx - Nitric acid data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.8294 || RH > 0.9949
    error("RH %.4f is outside valid range [%.4f, %.4f] for HNO3", RH, 0.8294, 0.9949);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = 1.031175777640;
A_3 = -4.068871311689;
A_2 = -1.219295149493;
A_1 = -0.780026716912;
A_0 = 0.999863154244;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0063, 0.1590, 0.0826);

end
