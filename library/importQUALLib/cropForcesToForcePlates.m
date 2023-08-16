function [COPcropped,ForceCropped] = cropForcesToForcePlates(COPs,Forces,Locations,Fill)
%% Set COP values outside the force plate edges to the Fill level
% Fill == 0 is at the force plate center
% 
% SYNTAX
% [COPcropped,ForceCropped] = cropForcesToForcePlates(COP,Forces,Locations,Fill)
% 
%   COPs        (Plates x 3d x NSamples double) non-cropped COP
%   Forces      (double Plates x 3d x NSamples) optional argument
%   Locations   (double Plates x 4-corners x 3d) locations of the corners of the force plates [m, meter]
%   Fill        (double) [default NaN] COP value for non-loaded phases (0=force plate center)
%            
% OUTPUT
%   COPcropped  (Plates x 3d x NSamples double) cropped COP
%   Forces      (double Plates x 3d x NSamples) optional
%
% EXAMPLES
% - Crop COP and force:
% [COPcropped,ForceCropped] = cropForcesToForcePlates(COP,Forces,Locations,0)
% - Only crop COP:
% COPcropped = cropForcesToForcePlates(COP,[],Locations)
%
% Local functions: 
% 
% See also: getCOPfromAnalog
% 
% (c) 2023 by Predimo GmbH
% Author: Marc de Lussanet
% 
% version 230803 (MdL) first version

% Arguments
narginchk(3,4);
if nargin<4 || isempty(Fill),   Fill =  0;   end

% Init
X=1;Y=2;Z=3;
NForcePlates = size(COPs,1);
COPcropped = COPs;
ForceCropped = Forces;

% Loop the force plates
for FPNo = 1:NForcePlates
    Location = squeeze(Locations(FPNo,:,:));
    MinEdges = min(Location,[],1);
    MaxEdges = max(Location,[],1);
    LcX      = mean(Location(:,X)); % plate center in X
    LcY      = mean(Location(:,Y)); % plate centre in Y
    % Check if the COP is outside the edges of the force plate
    IsOutside = ...
        COPs(FPNo,X,:) < MinEdges(X) | COPs(FPNo,Y,:) < MinEdges(Y) | ...
        COPs(FPNo,X,:) > MaxEdges(X) | COPs(FPNo,Y,:) > MaxEdges(Y);
    % If the ground level is not at zero, the COP can be outside the edges of the plate
    IsAtGroundBase = abs(COPs(FPNo,Z,:)) < 1e-10;
    DoOmit = IsOutside & IsAtGroundBase;
    % Omit by setting the COP to the desired position
    COPcropped(FPNo,X,DoOmit) = Fill + LcX;
    COPcropped(FPNo,Y,DoOmit) = Fill + LcY;
    % Crop the forces by setting them to zero
    if ~isempty(ForceCropped)
        ForceCropped(FPNo,X,DoOmit) = 0;
        ForceCropped(FPNo,Y,DoOmit) = 0;
    end
end
end