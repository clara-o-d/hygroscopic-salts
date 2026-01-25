function yi = interp1_with_nearest(x, y, xi)
    % Interpolate with pchip, then use nearest neighbor for extrapolation
    % to avoid NaN values
    
    yi = interp1(x, y, xi, 'pchip');
    
    % For any NaN values (extrapolation), use nearest neighbor
    nan_idx = isnan(yi);
    if any(nan_idx)
        yi(nan_idx) = interp1(x, y, xi(nan_idx), 'nearest', 'extrap');
    end
end
