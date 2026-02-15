function mf = calculate_mf_HClO4(RH)
% calculate_mf_HClO4
% Calculates mass fraction of HClO4 based on Relative Humidity.
% Range: 0.5051 < RH < 0.9949
% Data source: Paper 25 in sources_list_.xlsx - Perchloric acid data at 25°

% Add util folder to path if robust_fzero is not found
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end

% Check RH bounds
if RH > 1
    error("RH should be 0 < RH < 1");
end
if RH < 0.5051 || RH > 0.9949
    error("RH %.4f is outside valid range [%.4f, %.4f] for HClO4", RH, 0.5051, 0.9949);
end

% Coefficients from polynomial fit: RH = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4
A_4 = 3.344100308206;
A_3 = -5.944831210097;
A_2 = -0.343637298039;
A_1 = -0.526099514446;
A_0 = 1.000400113623;

% Polynomial function
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;

% Solve for mass fraction using robust_fzero
mf = robust_fzero(f, 0.0099, 0.3761, 0.1930);

end
