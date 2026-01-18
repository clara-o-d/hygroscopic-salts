function mf = calculate_mf_AgNO3_(RH)
% This function calculates the mass fraction of Silver Nitrate as a
% function of the Relative Humidity at a temperature of 25C
% Fit on water activity data at 25C from: 
if RH > 1 
    error("RH should be 0 < RH < 1")
end 
A_4 = -14.18; 
A_3 = 37.82; 
A_2 = -37.46;
A_1 = 16.16; 
A_0 = -1.565;
f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.4442, 0.9464, 0.70);
if mf > 0.9464
    error("below deliquescence relative humidity")
end 
end