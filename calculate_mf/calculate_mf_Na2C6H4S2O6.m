function mf = calculate_mf_Na2C6H4S2O6(RH)
% calculate_mf_Na2C6H4S2O6
% Calculates mass fraction of Na2C6H4S2O6 based on Relative Humidity.
% Range: 0.828 < RH < 1.0

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.828 || RH > 1.0
    error("RH %.4f is outside valid range [%.4f, %.4f] for Na2C6H4S2O6", RH, 0.828, 1.0);
end

% Coefficients
A_4 = -0.5676;
A_3 = -0.3648;
A_2 = -0.1752;
A_1 = -0.1635;
A_0 = 1.0;

% Polynomial function: RH = A0 + A1*x + ... + A4*x^4
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0003, 0.4585, 0.2294);

end
