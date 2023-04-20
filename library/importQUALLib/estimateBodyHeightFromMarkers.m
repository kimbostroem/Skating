function Height = estimateBodyHeightFromMarkers(Segments,DefaultHeight,Window)
%% Estimate the body height from the head markers
%
% SYNTAX
% Height = estimateBodyHeightFromMarkers(Segments,DefaultHeight,Window)
%
% INPUT
%     Segments      (struct) Coordinate systems of the segments
%     DefaultHeight (double) default height from general.def_body.xlsx
%     Window        (double) the samples during which a stable pose was detected
%
% OUTPUT
%     Height        (double) estimated actual height
%
% Local functions:
%
% See also: extractBodyFromSegments, calcMainCS
% 
% (c) 2023 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230215 (MdL) fix in warning for unplausible height; workaround for higher ground level
% version 230221 (MdL) outsourced from extractBodyFromSegments
% version 230222 (MdL) use lowest position across markers rather than toes segment

% stop if no markers are present
Height = DefaultHeight;
MinPlausibleHeight = 1;
MaxPlausibleHeight = 2.5;
MaxPlausibleFloorOffset = 0.1;

% Check if the necessary information is present
if size(Segments.skull,1) < 5 || all(isnan(Segments.skull(5,:,:)),'all')
    warning('Not enough head markers found: using default value of %.2f',Height);
    return;
end

% Get Top
HeadTop = squeeze(Segments.skull(5,:,:));
% Get Floor
if isfield(Segments,'bottom_position')
    LowestMarker = squeeze(Segments.bottom_position(4,3,:));
    Floor = mean(LowestMarker(Window),'omitnan') - 0.06; % subtract 6 cm for marker size and foot thickness
    % The might be a remaining error due to shoe wear or missing toe markers etc.
    if abs(Floor) < MaxPlausibleFloorOffset
        Floor=0;
    end
else
    Floor = 0;
end


% get the height from the vertical component of the top marker
HeightMean = mean(HeadTop(3,Window),'omitnan');
HeightStd = std(HeadTop(3,Window),'omitnan');
if Floor > MaxPlausibleFloorOffset
    warning('compensating body height for higher floor level (%.2fm)',Floor);
    HeightMean = HeightMean - Floor;
end

IsPlausible = ...
    HeightMean < MaxPlausibleHeight && ...
    HeightMean > MinPlausibleHeight && ...
    HeightStd  < 0.05;
if IsPlausible 
    Height  = round(HeightMean,2); % round to cm
    fprintf('Estimated height %.0f +/- %.1fcm\n',100*Height,100*HeightStd);
else
    warning('The estimated height %.2f +/- %.2fm is not plausible -> using default value of %.2f\n',HeightMean,HeightStd,Height);
end
end
