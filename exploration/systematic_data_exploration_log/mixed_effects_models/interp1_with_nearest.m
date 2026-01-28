function yi = interp1_with_nearest(x, y, xi)
    % Interpolate with pchip or linear, then use nearest neighbor for extrapolation
    % to avoid NaN values
    
    % Check if we have enough points
    if length(x) < 2
        % Not enough points, return constant value
        yi = repmat(median(y), size(xi));
        return;
    elseif length(x) >= 4
        % Use pchip for smooth interpolation (requires at least 4 points)
        yi = interp1(x, y, xi, 'pchip');
    else
        % Use linear interpolation if we have 2-3 points
        yi = interp1(x, y, xi, 'linear');
    end
    
    % For any NaN values (extrapolation), use nearest neighbor
    nan_idx = isnan(yi);
    if any(nan_idx)
        yi(nan_idx) = interp1(x, y, xi(nan_idx), 'nearest', 'extrap');
    end
    
    % Final fallback: if still NaN, use median
    if any(isnan(yi))
        yi(isnan(yi)) = median(y);
    end
end
