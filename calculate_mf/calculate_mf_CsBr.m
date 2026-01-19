function mf = calculate_mf_CsBr(RH)
% Fit on data range: MassFrac [0.7979, 0.9367] -> RH [0.8475, 0.9482]
if RH < 0.8475
    error("Input RH (%.4f) is below the lower fit limit for CsBr (0.8475)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = 14.805575;
A_3 = -107.715189;
A_2 = 207.749596;
A_1 = -156.378018;
A_0 = 42.175648;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.7979, 0.9367, 0.87);
end
