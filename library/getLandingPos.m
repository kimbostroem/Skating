function landingPos = getLandingPos(COP, footPos)

minGapLength = 0.3;
nAvg = 20; % take median over these many samples after landing to get a more stable estimate of the landing position
nSamples = size(COP, 2);
isGap = any(isnan(COP));
% find the beginnings of each gap
gapStarts = [-1 find(isGap)];
dGapStarts = [1 diff(gapStarts)];
gapStarts(dGapStarts == 1) = [];
% find the endings of the gap
gapStops = [find(isGap) nSamples+2];
dGapStops = [diff(gapStops) 1];
gapStops(dGapStops == 1) = [];
% remove leading and trailing gaps
idx = (gapStarts <= 1 | gapStops+nAvg >= nSamples);
gapStarts(idx) = [];
gapStops(idx) = [];
% remove gaps that are too short
gapSizes = vecnorm(COP(:, gapStops+1) - COP(:, gapStarts-1));
idx = (gapSizes < minGapLength);
gapStops(idx) = [];

if ~isempty(gapStops)
    idxAvg = gapStops(1)+1:gapStops(1)+1+nAvg;
    landingPos = median(footPos(:, idxAvg), 2);
else
    landingPos = nan(3, 1);
end

end