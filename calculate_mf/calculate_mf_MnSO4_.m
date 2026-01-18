function mf = calculate_mf_MnSO4_(RH)
% Fit on data range: MassFrac [0.2020, 0.8351] -> RH [0.8620, 0.9961]
if RH < 0.8620
    error("Input RH (%.4f) is below the lower fit limit for MnSO4 (0.8620)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -5.731;
A_3 = 10.12;
A_2 = -6.464;
A_1 = 1.71;
A_0 = 0.8388;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.5);
end