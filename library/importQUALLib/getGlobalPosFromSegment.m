function PosGlobal = getGlobalPosFromSegment(Segment,LocalPosition,Which)
%% get the global position of a location on a segment, as defined by a CS
%
% SYNTAX
%   PosLocal = getGlobalPosFromSegment(Segment,LocalPosition)
%
% INPUT
%     Segment       (4x3[xNS] double) Coordinate system (rotation matrix + position); may be a 
%                   time series of NS samples
%     LocalPosition (1x3 double) local position on the segment
%     Which         (char) 'base' (default) or 'end'
%
% OUTPUT
%    PosGlobal      (3xNS double) global position of a point with LocalPosition on the segment 
%
% See also: calcMainCS, createCS
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Version 320213 (MdL) allow position with respect to end of segment

% Detect the kind of input variable: if the size has two positions, it is a 1-sample CS (i.e.
% Coordinate System of rotation matrix + position vector); if size has three positions, it is 
% an array of Cs with length NSamples
IsTimeSeries = length(size(Segment)) == 3;

% Init
X=1;Y=2;Z=3;B=4;E=5;
if nargin == 3 && strcmp(Which,'end') 
    Ref = E;
else % default: with respect to segment base
    Ref = B;
end

% do the matrix multiplications and translation
if IsTimeSeries
    Reference = squeeze(Segment(Ref,:,:));
    Rx   = squeeze(Segment(X,:,:));
    Ry   = squeeze(Segment(Y,:,:));
    Rz   = squeeze(Segment(Z,:,:));
    PosGlobal = Reference + Rx * LocalPosition(X) + Ry * LocalPosition(Y) + Rz * LocalPosition(Z);
else
    Reference = Segment(Ref,:);
    R = Segment(X:Z,:);
    PosGlobal = Reference + R * LocalPosition;
end
end
