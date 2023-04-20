function [OutSize,IsEstID] = calcSegmentSize(Index,BasePos,FollowerPos,DefBody,SubjectHeight,Window,DoErrorPlot)
%% Calculate segment lengths (normalized to SubjectHeight)
% end point can be in fifth dimension of the BasePos segment
% BasePos, FollowerPos can be 3xNSamp or 4x3xNSamp or 5x3xNSamp
%
% SYNTAX
% [OutSize,IsEstID] = calcSegmentSize(Id,BasePos,FollowerPos,DefBody,SubjectHeight,Window,DoErrorPlot)
% [OutSize,IsEstID] = calcSegmentSize(Id,BasePos,FollowerPos,DefBody,SubjectHeight,Window)
%
% INPUT
%     Id          (double) index of current segment
%     BasePos     (double) base segment : 3xNSamp position or 4x3xNSamp or 5x3xNSamp segment
%     FollowerPos (double) follower segment : 3xNSamp position or 4x3xNSamp or 5x3xNSamp segment
%     DefBody     (struct) MNData structure
%     SubjectHeight (double) estimated subject height in m
%     Window      (double) 1-s Window of samples to use to estimate the segment length (selected by minimal movment)
%     DoErrorPlot (logical) flag if an error plot is to be made
%           
% OUTPUT
%     OutSize     (double) the length of the segment
%     IsEstID     (logical) has the segment length been estimated from the data
% 
% Local functions: 
%
% See also: estimateSegLengths, extractBodyFromSegments
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230213 (MdL & LK) use selected Window rather than omitting nan values; improved output

narginchk(6,7);
if nargin<7, DoErrorPlot = false; end
IsEstID = false;
OutSize = nan;

% if not a marker but a segment is provided, get the position of the segment's base
BaseIsSegment     = length(size(BasePos)) == 3;
FollowerIsSegment = length(size(FollowerPos)) == 3;
FollowerIsPosition= ~isempty(FollowerPos) && length(size(FollowerPos)) == 2;
BaseIsNaN         = all(isnan(BasePos),'all');
FollowerIsNaN     = all(isnan(FollowerPos),'all');
% End segments do not have a follower
IsDistalSegment   = isempty(FollowerPos) && BaseIsSegment && size(BasePos,1)>4;
NoFollowerSegment = FollowerIsNaN        && BaseIsSegment && size(BasePos,1)>4;

% define positions from segments, if necessary
if BaseIsNaN
    % if there is no base position, the length cannot be calculated
    return;
elseif FollowerIsSegment
    FollowerPos = squeeze(FollowerPos(4,:,:));
elseif FollowerIsPosition
    % no action needed
elseif IsDistalSegment
    FollowerPos = squeeze(BasePos(5,:,:));
elseif NoFollowerSegment
    FollowerPos = squeeze(BasePos(5,:,:));
else
    % if there is no follower position, the length cannot be calculated
    return;
end

% convert base segment to position, if necessary
if BaseIsSegment
    BasePos = squeeze(BasePos(4,:,:));
end

% select the time window for which to estimate the length
BasePos     = BasePos(:,Window);
FollowerPos = FollowerPos(:,Window);

% get the mean and standard deviation of the length
Distance_Window = vecnorm(FollowerPos - BasePos, 2);
Distance = mean(Distance_Window,'omitnan');
StDev    = std( Distance_Window,'omitnan');

%% checks
ScaleToDefault = Distance / (DefBody.segments.size(Index) * DefBody.general.height);
HasLength = Distance>0.001; % length is more than 1 mm, i.e. base and follower aren't the same position
% Variable length : 
% - standard deviation > 5 cm, or
% - standard deviation > 10 %
HasVariableLength = (StDev>0.05 || StDev/Distance > 0.1);
if ~HasLength
    fprintf('Size of %s could not be determined\n',DefBody.segments.name{Index});
    return;
elseif HasVariableLength % check on variable length
    fprintf('Highly variable size of %s: %.1f +/- %.1fcm (%.0f%%) Scale %.1f\n',DefBody.segments.name{Index},Distance*100,StDev*100,100*StDev/Distance,ScaleToDefault);
    if DoErrorPlot
        figure; plot(sqrt(sum((FollowerPos-BasePos).^2))); 
        title(sprintf('calcSegmentSize: variable length of %s',strrep(DefBody.segments.name{Index},'_','\_')))
    end
else
    if  contains(DefBody.segments.name{Index},'pelvis')
        fprintf('Estimated width of  %22s: %4.1f +/- %.1fcm (%.0f%%) Scale %.2f\n',DefBody.segments.name{Index},Distance*100,StDev*100,100*StDev/Distance,ScaleToDefault);
    else
        fprintf('Estimated length of %22s: %4.1f +/- %.1fcm (%.0f%%) Scale %.2f\n',DefBody.segments.name{Index},Distance*100,StDev*100,100*StDev/Distance,ScaleToDefault);
    end
end
IsEstID = true;
OutSize = Distance / SubjectHeight;
end
