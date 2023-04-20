function LocalPos = getLocalPosOnSegment(Segment,GlobalPosition,DoInvert)
%% get the global position of a location on a segment, as defined by a CS
%
% SYNTAX
%   PosLocal = getLocalPosOnSegment(Segment,GlobalPosition)
%
% INPUT
%     Segment        (4x3[xNS] double) Coordinate system (rotation matrix + position); may be a 
%                    time series of NS samples
%     GlobalPosition (1x3[xNS] double) local position on the segment
%     DoInvert       (char array) if 'invert', then invert the rotation matrix
%
% OUTPUT
%    LocalPos  (3xNS double) Local position of point in segment coords
%
% See also: calcMainCS, getGlobalPosFromSegment
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Version 230206 (MdL) optional parameter DoInvert

% Init
narginchk(2,3);
if nargin < 3
    DoInvert = false;
else
    DoInvert = matches(DoInvert,'invert');
end

IsTimeSeries = length(size(Segment)) == 3;
X=1;Y=2;Z=3;B=4;
if IsTimeSeries
    Base = squeeze(Segment(B,:,:));
    if DoInvert
        Rx   = squeeze(Segment(X,:,:));
        Ry   = squeeze(Segment(Y,:,:));
        Rz   = squeeze(Segment(Z,:,:));
    else
        Rx   = squeeze(Segment(X:Z, X, :));
        Ry   = squeeze(Segment(X:Z, Y, :));
        Rz   = squeeze(Segment(X:Z, Z, :));
    end
    PosOffset = GlobalPosition - Base;
    LocalPos = Rx .* PosOffset(X,:) + Ry .* PosOffset(Y,:) + Rz .* PosOffset(Z,:);
else
    Base = Segment(B,:);
    if DoInvert
        R = Segment(X:Z,:)';
    else
        R = Segment(X:Z,:);
    end
    LocalPos = R * (GlobalPosition - Base)';
end
end
