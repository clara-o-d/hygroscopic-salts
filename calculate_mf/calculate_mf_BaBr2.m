function mf = calculate_mf_BaBr2(RH)
% Fit on data range: MassFrac [0.8307, 0.9434] -> RH [0.8221, 0.9587]
if RH < 0.8221
    error("Input RH (%.4f) is below the lower fit limit for BaBr2 (0.8221)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -2298.823313;
A_3 = 8034.982435;
A_2 = -10533.863744;
A_1 = 6138.288520;
A_0 = -1340.357793;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.8307, 0.9434, 0.89);
end
