function StablePeriod = getStablePeriod(Markers,Present,Info,DoPlot)
%% find the period of 1 second when the movements are slowest
%
% SYNTAX
% StablePeriod = getStablePeriod(Markers,Present,Info,DoPlot)
%
% INPUT
%     Markers      (struct of 3 x N) time series for each marker that is present
%     Present      (struct of logical) for each marker: is it present
%     Info         (struct) timing information
%     DoPlot       (logical) 
%
% OUTPUT
%     StablePeriod (double array) indices of the selected period
%
% Local functions: 
%
% See also: extractBodyFromSegments, calcMainCS
% 
% (c) 2023 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230208 (MdL) corrected border condition

% Init
narginchk(3,4);
if nargin < 4, DoPlot = false; end
IsPresent  = cell2mat(struct2cell(Present));
MarkerData = struct2cell(Markers);
MarkerData = MarkerData(IsPresent);
NMarkers   = length(MarkerData);
NSamples   = size(MarkerData{1},2);
Speed      = zeros(NMarkers,NSamples);
Time       = Info.QTMTimeSRint;
Frequency  = Info.QTMFreqSRint;
FilterSpan = 0.1; % s
OneSecond  = Frequency;
HalfSecond = floor(Frequency/2);

% Get the velocity and speed of each of the segments
for Marker=1:NMarkers
    Pos = squeeze(MarkerData{Marker});
    [~,Velocity] = filterData(Time, Pos, FilterSpan, 'NoFill', 2, 2);
    Velocity(isnan(Velocity)) = 0; % ignore samples with gaps
    Speed(Marker,:) = sqrt(sum(Velocity.^2)); % speed is the tangential of the velocity
end
% ignore periods where all segments are empty
Speed(all(Speed==0)) = max(Speed,[],'all');

% find the 1-s window during which the maximum speed (across all segments) is lowest
[~,MinimumInMaximumSpeed] = min(movmean(max(Speed),OneSecond));
if MinimumInMaximumSpeed < 1 + HalfSecond
    Start = 1; 
    % if the measurement is shorter that a second, use the entire measurement
    End   = min(OneSecond, NSamples); 
elseif MinimumInMaximumSpeed >= NSamples - HalfSecond
    % if the measurement is shorter that a second, use the entire measurement    
    Start = max(1, NSamples-OneSecond);
    End   = NSamples;
else
    Start = MinimumInMaximumSpeed - HalfSecond + 1;
    End   = MinimumInMaximumSpeed + HalfSecond;
end
StablePeriod = Start:End;

% check with a plot
if DoPlot
    figure; hold on; title('extracBodyFromSegments: selection of stable phase')
    xlabel('time (s)')
    ylabel('speed (m/s)')
    for Marker=1:NMarkers
        plot(Time,Speed(Marker,:))
    end
    xline(Time(Start));
    xline(Time(End));
end
end
