function mf = calculate_mf_LiNO3(RH)
% Fit on data range: MassFrac [0.0257, 0.6130] -> RH [0.7353, 0.9967]
if RH < 0.7353
    error("Input RH (%.4f) is below the lower fit limit for LiNO3 (0.7353)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -6.794033;
A_3 = 7.459413;
A_2 = -3.060065;
A_1 = 0.250592;
A_0 = 0.987894;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0257, 0.6130, 0.32);
end
