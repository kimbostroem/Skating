function CS_out = rotateCS(CS_in,Angle,CS_about,RotateAbout)
%% Rotate coordinate system (time series) CS_in about axis x, y, or z of CS_About by Angle
% This function is about TEN TIMES FASTER than a loop of matrix multiplications
% 
% SYNTAX
% CS_out = rotateCS(CS_in,Angle,CS_about,RotateAbout);
%
% INPUT
%     CS_in       ([3+]x[3]x[nSamples] double) time series of nSamples of a 3x3 rotation matrix; there
%                 may be a fourth (and fifth) column (4,:,:), the origin of the cs
%     Angle       (double) angle [deg]
%     CS_about    (3x3 or [3+]x[3]x[nSamples] double) the refrence CS, which may also be a stationary
%                 rotation matrix
%     RotateAbout (char) 'x', 'y', or 'z'
%           
% OUTPUT
%     CS_out      (3+,3,nSamples double -> dimension of CS_in) 
%
% EXAMPLE 
% 
% See also: calcJointStream, calcMainCS
% 
% (c) 2020 by Predimo
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% 230129 (MdL & LK) added header and comments


% Check Dimensions
DimsIn = size(CS_in);
DimsAbout = size(CS_about);
% the first D must be 3<=D<=5; the second D==3 and the third (NSamples) must exist
if DimsIn(1)<3 || DimsIn(1)>5 || DimsIn(2) ~=3 || length(DimsIn) ~= 3
    error('ill-shaped CS %dx%dx%d',DimsIn(1),DimsIn(2),length(DimsIn));
end
% require that NSamples >3 so that we are not confusing the structure of the input variable
NSamples = DimsIn(3);
if NSamples<4
    error('ill-shaped CS (%d =too few samples)',DimsIn(3));
end

% copy the coordinates of the origin
CS_out = CS_in;

% if the CS_about = stationary, repeat it for the length of the time series
if length(DimsAbout) == 2
    CS_about = repmat(CS_about,[1 1 NSamples]);
end

% axis is the axis of rotation
switch RotateAbout
    case 'x',  axis = squeeze(CS_about(1,1:3,:));
    case 'y',  axis = squeeze(CS_about(2,1:3,:));
    case 'z',  axis = squeeze(CS_about(3,1:3,:));
    % case 'x',  axis = squeeze(CS_about(1:3,1,:));
    % case 'y',  axis = squeeze(CS_about(1:3,2,:));
    % case 'z',  axis = squeeze(CS_about(1:3,3,:));
    otherwise, error('non-defined axis "%s"',RotateAbout);
end

%rotate each of the co-ordinate axes about the axis
for Ax=1:3
    CS_out(Ax,:,:) = rotateVectorDeg(squeeze(CS_in(Ax,1:3,:)),axis,Angle);
    % CS_out(:,Ax,:) = rotateVectorDeg(squeeze(CS_in(1:3,Ax,:)),axis,Angle);
end
end

%% =================================================================================================

function out = rotateVectorDeg(v,axis,angdeg)
% rotate vector v about the vector "axis" by an angle angdeg (degrees)
% generalized for time series of vectors
% v and axis may be of size [3,ns], [1,3] or [3,1]
% 
% (c) 2020 by Predimo
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% 200529 (MdL) renamed function

% check dimensionality
szv=size(v); sza=size(axis);
if any(sum(szv-sza)) || (sza(1)~=3 && sza(2)~=3)
    error('unexpected dimensionality');
end

% get the axes and rotations
axis = normLength(axis);
% dot product gives the cos of the angle between the vectors
vn = dot(axis,v).*axis;
vx = v - vn;
vy = cross(axis,vx);
out = vn + vx*cosd(angdeg) + vy*sind(angdeg);
end

%% =================================================================================================






%% zo zag de functie er eerst uit
% % this version calculates inverted rot mat
% CS_in    = permute(CS_in(   1:3,:,:), [2,1,3]);
% CS_about = permute(CS_about(1:3,:,:), [2,1,3]);
% 
% % axis is the axis of rotation
% switch RotateAbout
%     case 'x',  axis = squeeze(CS_about(1,:,:));
%     case 'y',  axis = squeeze(CS_about(2,:,:));
%     case 'z',  axis = squeeze(CS_about(3,:,:));
%     otherwise, error('non-defined axis "%s"',RotateAbout);
% end
% 
% %rotate each of the co-ordinate axes about the axis
% for Ax=1:3
%     CS_out(Ax,:,:) = rotVecd(squeeze(CS_in(Ax,:,:)),axis,Angle);
% end
% 
% % bring to the original form:
% CS_out(1:3,:,:) = permute(CS_out(1:3,:,:), [2,1,3]);



%% test sequence for checking the time
% CS_in = rand([3 3 1000000]);
% CS_about = rand([3 3 1000000]);
% Angle = 11;
% xyz = 'x';
% tic
% CS_out = rotateCS(CS_in,Angle,CS_about,xyz);
% toc
% 
% CS_about = rand([3 3]);
% tic
% CS_out = rotateCS(CS_in,Angle,CS_about,xyz);
% toc

% for comparison: 
%% a version with a conventional loop is a factor 10 slower:

% % Init
% CS_out = CS_in;
% 
% % RotateAbout is the axis of rotation
% switch RotateAbout
%     case 'x', R_axis = [1 0 0;    0 cosd(Angle) -sind(Angle);    0 sind(Angle) cosd(Angle)];
%     case 'y', R_axis = [cosd(Angle) 0 sind(Angle);    0 1 0;    -sind(Angle) 0 cosd(Angle)];
%     case 'z', R_axis = [cosd(Angle) -sind(Angle) 0;    sind(Angle) cosd(Angle) 0;    0 0 1];
%     otherwise, error('Wrong value "%s" for RotateAbout: expected "x", "y" or "z"',RotateAbout);
% end
% 
% % loop the samples
% for i=1:length(CS_in)
%     R_in_i = squeeze(CS_in(1:3,1:3,i));
%     R_about_i = squeeze(CS_about(1:3,1:3,i)); 
%     % get R_in relative to parent CS
%     R_relative_i = R_about_i' * R_in_i;
%     % rotate about axis by angle
%     CS_r = R_relative_i * R_axis;
%     % rotate back from parent CS
%     CS_out(1:3,:,i) = CS_r * R_about_i;
% end
% end




