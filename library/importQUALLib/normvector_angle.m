function [nvec,angle] = normvector_angle(p1,p2,p3)
%% Returns norm vector to the plane defined by points p1,p2,p3, and the angle p1-p2-p3 in degrees
% p1,p2,p3 may be of size [3,ns], [1,3] or [3,1]
% if p3 is empty or missing, p1 and p2 are treated as vectors v1 and v2 for the cross product
% 
% INPUT
%     p1 (double) [3,ns], [1,3] or [3,1] 3D point or vector (or time series thereof)
%     p2 (double) [3,ns], [1,3] or [3,1] 3D point or vector (or time series thereof)
%     p3 (double) [3,ns], [1,3] or [3,1] 3D point (or time series thereof)
%
% OUTPUT
%     nvec  (double) 3D vector (or time series thereof): normal vector to plane p1-3
%     angle (double) angle (time series) of p1-p2-p3
%
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 20200529 (MdL)
% version 221029 (MdL) allow to vectors as input; added header

% Init
narginchk(2,3);
if nargin<3
    p3 = [];
end
% check dimensionality of the input:
if isempty(p1) || isempty(p2)
    nvec = [];
    angle= [];
    return;
end
sz1 = size(p1);
sz2 = size(p2);
if sum(sz1-sz2) || (sz1(1)~=3 && sz1(2)~=3)
    error('unexpected dimensionality');
end

% the two vectors defined by the three points
if ~isempty(p3)
    v1   = normLength(p1-p2);
    v2   = normLength(p3-p2);
else
    % if p3 was not given, then take p1 and p2 as vectors
    v1 = p1;
    v2 = p2;
end
nvec = normLength(cross(v1,v2));

if nargout>1
    % Geometrically, the dot product is the cosine of the angle between two vectors
    angle= acosd(dot(v1,v2));
end
end
