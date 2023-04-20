function X = naninterp(X)
% Interpolate over NaNs, extrapolate with mean value.
% See INTERP1 for more info
M = mean(X,'omitnan');
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'linear',M);
end
