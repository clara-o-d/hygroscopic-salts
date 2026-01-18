function mf = calculate_mf_CuSO4_(RH)
% Fit on data range: MassFrac [0.2205, 0.6644] -> RH [0.9750, 0.9963]
if RH < 0.9750
    error("Input RH (%.4f) is below the lower fit limit for CuSO4 (0.9750)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -0.7273;
A_3 = 1.019;
A_2 = -0.5633;
A_1 = 0.1177;
A_0 = 0.9885;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.4);
end