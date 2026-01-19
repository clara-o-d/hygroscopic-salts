function mf = calculate_mf_NaClO4(RH)
% Fit on data range: MassFrac [0.1586, 0.8033] -> RH [0.7775, 0.9869]
if RH < 0.7775
    error("Input RH (%.4f) is below the lower fit limit for NaClO4 (0.7775)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -3.344665;
A_3 = 5.066175;
A_2 = -2.986222;
A_1 = 0.670195;
A_0 = 0.935590;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.1586, 0.8033, 0.48);
end
