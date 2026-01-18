function mf = calculate_mf_KNO3(RH)
% Fit on data range: MassFrac [0.0537, 0.6651] -> RH [0.9315, 0.9967]
if RH < 0.9315
    error("Input RH (%.4f) is below the lower fit limit for KNO3 (0.9315)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = 0.012512;
A_3 = -0.083721;
A_2 = -0.017550;
A_1 = -0.057961;
A_0 = 0.999924;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.0537, 0.6651, 0.36);
end
