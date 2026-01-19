function mf = calculate_mf_CsCl_(RH)
% Add util folder to path if needed
if ~exist('robust_fzero', 'file')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'util'));
end
% This function calculates the mass fraction of Cesium Chloride as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = -6.013; 
A_3 = 11.96; 
A_2 = -8.614;
A_1 = 2.566; 
A_0 = 0.7265;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.2394, 0.9042, 0.57);
% if mf > 0.9042
%     error("below deliquescence relative humidity")
% end 
end