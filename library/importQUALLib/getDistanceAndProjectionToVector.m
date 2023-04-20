function [Distance,Projection] = getDistanceAndProjectionToVector(A,B,C)
%% Get the distances of 
%% – point C to line AB and 
%% - the height line to A
%
% SYNTAX
%     [Distance,Projection] = getDistanceAndProjectionToVector(A,B,C)
%
% INPUT
%     A   (3xNS double) time series of Origin, point A
%     B   (3xNS double) time series of second point on line, B
%     C   (3xNS double) time series of point C
%
% OUTPUT
%    Distance   (double) median distance SC
%    Projection (double) median distance AS
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Version 221117 first version

%% Triangle ABC, get distance of line AB to C
AC = C - A;
AB = B - A;

%% get point S (SC being the height line of ABC)
% where P is the relative location of S on AB
P     = sum(AB .* AC,'omitnan') ./ sum(AB .* AB,'omitnan');
S     = A + P .* AB; % point S (SC being the height line of ABC)
SClen = veclen(C-S); % length of SC
ASlen = veclen(S-A);
Distance   = median(SClen,'omitnan');
Projection = median(ASlen,'omitnan');

% check
ABlen = median(veclen(AB),'omitnan');
SCvar = std(SClen);
ASvar = std(ASlen);
Threshold = 0.01;
if SCvar / ABlen > Threshold || ASvar / ABlen > Threshold 
    figure; hold on;plot(SClen);plot(ASlen); legend('SClen','ASlen')
    title('highly variable result in getDistanceAndProjectionToVector')
end
end

function [Len] = veclen(Vector3D,Dim)
if nargin==1, Dim=1; end
X=1;Y=2;Z=3;
if Dim==1
    Len = sqrt(Vector3D(X,:).^2 + Vector3D(Y,:).^2 + Vector3D(Z,:).^2);
else
    Len = sqrt(Vector3D(:,X).^2 + Vector3D(:,Y).^2 + Vector3D(:,Z).^2);
end
end
