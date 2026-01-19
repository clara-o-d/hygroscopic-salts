function mf = calculate_mf_KBr(RH)
% Fit on data range: MassFrac [0.4966, 0.7737] -> RH [0.8325, 0.9528]
if RH < 0.8325
    error("Input RH (%.4f) is below the lower fit limit for KBr (0.8325)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -116.205197;
A_3 = 285.586109;
A_2 = -262.274244;
A_1 = 106.397529;
A_0 = -15.112296;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.4966, 0.7737, 0.64);
end
