function IsLoaded = computeLoadedPhases(LoadZ, LoadThresholds, MinimalLoadDuration, PlateNo)
%% Generate logical array if current plate is loaded at each sample
% LoadZ can be positive or negative (the sign is detected automatically).
% 
% Since we are dealing with bit-noise, the noise depends on amplification. Therefore, it is best to 
% apply the threshold on the analog z-channels, if possible (e.g. for Kistler type force plates)
% 
% The beginning and end of loaded periods can be difficult to discern in the presence of noise. 
% - The first-order periods are all that are above the high threshold. 
% - Then, the loaded phases are expanded to the periods for which they are still above the lower 
%   threshold.
% - Finally, loaded periods that are shorter than MinimalLoadDuration are omitted.
% 
% SYNTAX
% IsLoaded = computeLoadedPhases(LoadZ, LoadThresholds, MinimalLoadDuration);
% IsLoaded = computeLoadedPhases(LoadZ, [], MinimalLoadDuration);
% IsLoaded = computeLoadedPhases(LoadZ);
%
% INPUT
%     LoadZ               (double array) this can be the z-force or the analog signal
%                         representing the z-force
%     LoadThresholds      (2 double) first- and second order threshold (default [20 10] N)
%     MinimalLoadDuration (double) the minimal no of samples that a loaded period counts as loaded
%     PlateNo             (double) number of the plate
%            
% OUTPUT
%     IsLoaded  (logical array) for each sample: is the plate loaded
%
% EXAMPLE 
% IsLoaded = computeLoadedPhases(Z_component_of_force, [20 10], 0.1*Frequency);
% IsLoaded = computeLoadedPhases(Sum_of_analog_channels_in_Z, [1st_order_threshold 2nd_order_threshold], 0.1*Frequency);
%
% See also: getCOPfromAnalog
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Last revision: 221222 (MdL) minor changes review

%% Init & checks
narginchk(1,4);
if nargin<2 || isempty(LoadThresholds),      LoadThresholds = [40 10]; end
if nargin<3 || isempty(MinimalLoadDuration), MinimalLoadDuration = 100; end
if nargin<4 || isempty(PlateNo),             PlateNo = 0; end
NSamples = length(LoadZ);

% check if LoadThresholds are correct
if any(LoadThresholds<0)
    LoadThresholds = -LoadThresholds;
end
if ~all(LoadThresholds>0)
    error('LoadThresholds (%f, %f) must all be >0',LoadThresholds(1),LoadThresholds(2));
end
if LoadThresholds(1) < LoadThresholds(2)
    Tmp = LoadThresholds(1);
    LoadThresholds(1) = LoadThresholds(2);
    LoadThresholds(2) = Tmp;
end

% Initialize IsLoaded and return if the force plate is never loaded
IsLoaded = false(size(LoadZ));
IsNeverLoaded = all(abs(LoadZ)<LoadThresholds(1));
if IsNeverLoaded
    return;
end

% Make LoadZ positive
SignOfLoaded  = (abs(min(LoadZ)) < abs(max(LoadZ))) * 2 - 1;
LoadZ = SignOfLoaded * LoadZ;

% Return if the forceplate is continuously loaded
IsLoaded = LoadZ > LoadThresholds(1);
IsCompletelyLoaded = all(IsLoaded);
if IsCompletelyLoaded
    return;
end

%% expand the loaded period to the 2nd-order threshold
% the low threshold might sometimes include noise, so follow the more conservative expansion method
% find first and last samples of loaded periods
LoadSt  = find(diff(IsLoaded)== 1)+1; % 1st samples of loaded periods
LoadEnd = find(diff(IsLoaded)==-1);   % last samples of loaded periods
% expand backwards from first samples
for i=1:length(LoadSt)
    fl = LoadSt(i);  Lw = find(LoadZ(1:fl)   < LoadThresholds(2),1,'last');
    if ~isempty(Lw), LoadSt(i)=Lw; end
end
% expand forwards from last samples
for i=1:length(LoadEnd)
    fl = LoadEnd(i);  Hg = find(LoadZ(fl:end) < LoadThresholds(2),1,'first')+fl-1;
    if ~isempty(Hg), LoadEnd(i) =Hg; end
end

%% Omit very brief loaded periods
% (these are not plausible and are problematic for smoothing)
%
% beginnings and ends of loaded periods
if IsLoaded(1)==true,   LoadSt  = [false LoadSt]; end %#ok<*AGROW> % if already loaded at sample 1
if IsLoaded(end)==true, LoadEnd = [LoadEnd   NSamples]; end
% Durations of loaded periods
Duration = LoadEnd-LoadSt;
if length(LoadSt) ~= length(LoadEnd), error('found inconsistent loaded periods in plate %d',PlateNo); end
if any(Duration<0),                   error('found negative loaded periods in plate %d',PlateNo); end
% Omit very short loadings
if LoadSt(1)==0 && Duration(1)<MinimalLoadDuration
    Init = 0;
else
    Init = IsLoaded(1);
end
if LoadSt(1)==0
    LoadSt(1) = 1;
end
LoadSt( Duration<MinimalLoadDuration) = [];
LoadEnd(Duration<MinimalLoadDuration) = [];
% remove duplicate loaded periods
for ll=1:length(LoadSt)-1
    if LoadSt(ll) == LoadSt(ll+1) || LoadEnd(ll) == LoadEnd(ll+1)
        LoadSt(ll) = nan;
        LoadEnd(ll) = nan;
    end
end
LoadSt(isnan(LoadSt)) = [];
LoadEnd(isnan(LoadEnd)) = [];

% Regenerate the loaded list with omitted short loadings
IsLoaded = zeros(1,NSamples);
IsLoaded(LoadSt)  = 1;
IsLoaded(LoadEnd) = IsLoaded(LoadEnd) - 1;
IsLoaded = [Init cumsum(IsLoaded)];
IsLoaded(end)=[];
% plausibility check
if max(IsLoaded)>1 || min(IsLoaded)<0
    error('unexpected value for Loaded');
end
IsLoaded = logical(IsLoaded);

% at this point, it can happen that the plate is no longer loaded (if the loading was very brief)
end



% conditions tested for validation:
%           1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
% Loaded = [1 1 1 1 1 0 1 1 1 0 1 0 1 1 0 0 0 0 1]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [1 1 1 0 0 0 0 1 0 0 1 0 1 1 0 0 0 0 1]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [1 0 1 0 0 0 0 1 1 0 1 0 1 1 0 0 0 1 0]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [0 0 1 0 0 0 0 1 0 0 1 0 1 1 0 0 0 1 0]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;

