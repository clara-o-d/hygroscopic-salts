function mf = calculate_mf_NaBr(RH)
% Fit on data range: MassFrac [0.5407, 0.8243] -> RH [0.6133, 0.9290]
if RH < 0.6133
    error("Input RH (%.4f) is below the lower fit limit for NaBr (0.6133)", RH);
end
if RH > 1 
    error("RH should be 0 < RH < 1");
end 

A_4 = -130.776222;
A_3 = 334.244087;
A_2 = -320.221115;
A_1 = 135.668459;
A_0 = -20.466786;

f = @(xi) RH - A_0 - A_1.*xi - A_2.*xi.^2 - A_3.*xi.^3 - A_4.*xi.^4;
mf = robust_fzero(f, 0.5407, 0.8243, 0.68);
end
