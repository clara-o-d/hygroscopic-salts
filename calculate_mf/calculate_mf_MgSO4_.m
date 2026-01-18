function mf = calculate_mf_MgSO4_(RH)
% Fit on data range: MassFrac [0.1386, 0.7070] -> RH [0.9050, 0.9960]
if RH < 0.9050
    error("Input RH (%.4f) is below the lower fit limit for MgSO4 (0.9050)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -3.071;
A_3 = 4.028;
A_2 = -1.963;
A_1 = 0.3657;
A_0 = 0.9729;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.4);
end