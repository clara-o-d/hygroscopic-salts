function mf = calculate_mf_SrCl2(RH)
% Fit on data range: MassFrac [0.4801, 0.8171] -> RH [0.8059, 0.9778]
if RH < 0.8059
    error("Input RH (%.4f) is below the lower fit limit for SrCl2 (0.8059)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -71.923743;
A_3 = 179.896788;
A_2 = -168.213554;
A_1 = 69.355567;
A_0 = -9.634025;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.4801, 0.8171, 0.65);
end
