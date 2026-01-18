function mf = calculate_mf_CaNO3(RH)
% Fit on data range: MassFrac [0.1300, 0.8997] -> RH [0.6464, 0.9955]
if RH < 0.6464
    error("Input RH (%.4f) is below the lower fit limit for Ca(NO3)2 (0.6464)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -8.408009;
A_3 = 14.647404;
A_2 = -8.912854;
A_1 = 2.088723;
A_0 = 0.839698;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1300, 0.8997, 0.51);
end
