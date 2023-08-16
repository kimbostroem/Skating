function IsReactionForces = isReactionForce(Forces, LoadedThreshold)
%% check if it is plausible that the force plates register REACTION forces
% (i.e., loaded Z is positive)
% 
% SYNTAX
% IsReactionForces = isReactionForce(Forces, LoadedThreshold);
%
% INPUT
%     Forces          (NPlates x 3D x NSamples double) force data, or:
%                     (3D x NSamples double) force data
%     LoadedThreshold (double) optional threshold to decide of the plate(s) are loaded [N]
%            
% OUTPUT
%     IsReactionForces (logical) retruns true only if the forces are clearly reactional (i.e. if Z positive)
%
% EXAMPLE 
% IsReactionForces = isReactionForce(Forces, LoadedThreshold);
% IsReactionForces = isReactionForce(Forces);
%
% Local functions: none
%
% See also: kraftAusQualisys, computeLoadedPhases
% 
% (c) 2023 by Predimo GmbH
% Author: Marc de Lussanet
% Version 230519 (MdL) first version
% Version 230522 (MdL) fixed problems with dimensionality checks

% init
narginchk(1,2);
if nargin<2, LoadedThreshold = 10; end % N
IsReactionForces = false;
Z=3;

% Plausibility check
if ~isempty(Forces)
    Size = size(Forces);
    if length(Size) == 2
        Force = Forces;
    elseif length(Size) == 3
        % sum across all plates
        Force = squeeze(sum(Forces,'omitnan'));
    else % length of size must be 2 or 3
        error('isReactionForce: Dimensionality of Forces not allowed (length(size(Forces))==%d)',length(Size))
    end
    % check: the first dimension of Force should now be 3D
    Dim = size(Force);
    if Dim(1) == 1
        ForceZ = Force;
    elseif Dim(1) == 3
        ForceZ = Force(Z,:);
    else
        error('isReactionForce: Dimensionality of Forces not allowed (Dim(Force))==%d)',Dim(1))
    end

    % sign of vertical, Z-force component
    IsMeanPositive = mean(ForceZ,'omitnan') > 0;
    IsMaxPositive = max(ForceZ,[],'omitnan') > -min(ForceZ,[],'omitnan');
    IsLoaded = any(abs(ForceZ) > LoadedThreshold);
    % if both are positive and if there is load on the plates, we (probably) have reaction forces
    if IsLoaded && IsMeanPositive && IsMaxPositive
        IsReactionForces = true;
    end
end
end