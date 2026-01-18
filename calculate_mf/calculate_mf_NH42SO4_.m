function mf = calculate_mf_NH42SO4_(RH)
% Fit on data range: MassFrac [0.0884, 0.8289] -> RH [0.8310, 0.9959]
if RH < 0.8310
    error("Input RH (%.4f) is below the lower fit limit for (NH4)2SO4 (0.8310)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -2.179;
A_3 = 3.06;
A_2 = -1.563;
A_1 = 0.2564;
A_0 = 0.9821;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0, 1, 0.45);
end