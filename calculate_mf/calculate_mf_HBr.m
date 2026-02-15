function mf = calculate_mf_HBr(RH)
% calculate_mf_HBr
% Calculates mass fraction of HBr based on Relative Humidity.
% Range: 0.9437 < RH < 0.9949
% Data source: Paper 25 in sources_list_.xlsx - Hydrobromic acid data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.9437 || RH > 0.9949
    error("RH %.4f is outside valid range [%.4f, %.4f] for HBr", RH, 0.9437, 0.9949);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = -94.539636995614;
A_3 = 10.319767966967;
A_2 = -2.262092131945;
A_1 = -0.598417949706;
A_0 = 0.999826895839;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0080, 0.0749, 0.0414);

end
