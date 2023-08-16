function IsLoaded = computeLoadedPhases(LoadZ, LoadThresholds, MinimalLoadDuration, PlateNo)
%% Generate logical array if current plate is loaded at each sample
% We expect ground ACTION forces, i.e. negative during loaded phases
% 
% The beginning and end of loaded periods can be difficult to discern in the presence of noise. 
% - The first-order periods are all that are above the high threshold. 
% - Then, the loaded phases are expanded to the periods for which they are still above the lower 
%   threshold.
% - Finally, loaded periods that are shorter than MinimalLoadDuration are omitted. Such brief 
%   artificial loadings can occur at sudden landings and offsets, due to noise or "slipping" 
%   or "tapping" the force plate surface. 
% 
% SYNTAX
% IsLoaded = computeLoadedPhases(LoadZ, LoadThresholds, MinimalLoadDuration, PlateNo);
%
% INPUT
%     LoadZ               (double array) this can be the z-force or the analog signal
%                         representing the z-force
%     LoadThresholds      (2 double) first- and second order threshold (default [20 10] N)
%     MinimalLoadDuration (double) the minimal no of samples that a loaded period counts as loaded
%     PlateNo             (double) optional number of the plate (only for warning message)
%            
% OUTPUT
%     IsLoaded  (logical array) for each sample: is the plate loaded
%
% EXAMPLE 
% - to call with Rough and Fine threshods of 20 and 10 N respectively; minimal loading of 100 ms:
% IsLoaded = computeLoadedPhases(Z_component_of_force, [20 10], 0.1*Frequency);
% IsLoaded = computeLoadedPhases(LoadZ, [], MinimalLoadDuration);
% IsLoaded = computeLoadedPhases(LoadZ);
%
% See also: getCOPfromAnalog
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Last revision: 221222 (MdL) minor changes review
% Last revision: 230516 (MdL) do not dectect the sign automatically; updated header; 
%                             correction of -1 sample for LoadSt
% Last revision: 230522 (MdL) fixed unknown variable "Forces" (->LoadZ); omit non-used parameter DetectSign 
% Last revision: 230601 (MdL) - fixed error in addressed empty array
%                             - prevent identification of negative non-loaded periods
%                             - fixed exclusion of duplicate phases
%                             - header and comments
%                             - optional test plot

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

% check the sign of the forces
IsReactionForces = isReactionForce(LoadZ, LoadThresholds(2));
if IsReactionForces
    % warning('computeLoadedPhases: We expect ground ACTION force, but it looks like REaction force: correcting...');
    FReactZ = LoadZ;
else
    FReactZ = -LoadZ;
end

% Return if the forceplate is continuously loaded
IsLoaded = FReactZ > LoadThresholds(1);
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
    Flank = LoadSt(i);
    Low = find(FReactZ(1:Flank) < LoadThresholds(2),1,'last') - 1;
    if ~isempty(Low)
        LoadSt(i) = Low;
    end
end
% expand forwards from last samples
for i=1:length(LoadEnd)
    Flank = LoadEnd(i);
    High = find(FReactZ(Flank:end) < LoadThresholds(2),1,'first')+Flank-1;
    if ~isempty(High)
        LoadEnd(i) = High;
    end
end

%% Omit very brief loaded periods
% (these are not plausible and are problematic for smoothing)
%
% beginnings and ends of loaded periods
if IsLoaded(1)==true,   LoadSt  = [false LoadSt]; end %#ok<*AGROW> % if already loaded at sample 1
if IsLoaded(end)==true, LoadEnd = [LoadEnd   NSamples]; end
% Durations of loaded periods
Duration = LoadEnd-LoadSt;
% safety check for programming errors
if length(LoadSt) ~= length(LoadEnd), error('found inconsistent loaded periods in plate %d',PlateNo); end
if any(Duration<0),                   error('found negative loaded periods in plate %d',PlateNo); end

% deal with loaded periods that begin with the first sample
if isempty(LoadSt)
    Init = 0;
elseif LoadSt(1)==0 
    % plate is already loded on the first sample
    LoadSt(1) = 1;
    % remove the initial loading if it lasts very briefly
    if Duration(1)<MinimalLoadDuration
        Init = 0;
    else
        Init = 1;
    end
else
    Init = IsLoaded(1);
end

% remove brief loaded periods
LoadSt( Duration<MinimalLoadDuration) = [];
LoadEnd(Duration<MinimalLoadDuration) = [];

% remove duplicate loaded periods
if ~isempty(LoadSt)
    % case 1: the beginning and end are duplicate
    for ll=1:length(LoadSt)-1
        if LoadSt(ll) == LoadSt(ll+1) && LoadEnd(ll) == LoadEnd(ll+1)
            LoadSt(ll) = nan;
            LoadEnd(ll) = nan;
        end
    end
    % case 2: the beginning is duplicate (retain the last occurrence, i.e. the longest
    % of duplicate phases)
    for ll=1:length(LoadSt)-1
        if LoadSt(ll) == LoadSt(ll+1)
            LoadSt(ll) = nan;
            LoadEnd(ll) = nan;
        end
    end
    % case 3: the end is duplicate (count backwards and retain the first occurrence, i.e. the longest
    % of duplicate phases)
    for ll=length(LoadEnd) : -1 : 2
        if LoadEnd(ll) == LoadEnd(ll-1)
            LoadSt(ll) = nan;
            LoadEnd(ll) = nan;
        end
    end
    LoadSt(isnan(LoadSt)) = [];
    LoadEnd(isnan(LoadEnd)) = [];
end

% check for brief non-loaded periods
MinimalNonLoadDuration = 10; % samples (low value in order not to falsely detect loaded phases)
if length(LoadSt) >= 2
    DurationNL = LoadSt(2:end)-LoadEnd(1:end-1);
    ShortNonLoadIdx = find(DurationNL<MinimalNonLoadDuration);
    LoadSt( ShortNonLoadIdx +1) = [];
    LoadEnd(ShortNonLoadIdx)    = [];
    if any(DurationNL<1)
        warning('There seem to be periods of extremely rapid loading/unloading in the force: please check your data')
    end
end

% Recreate the loaded list with omitted short loadings
% at this point, it can happen that the plate is no longer loaded (if the loading was very brief)
IsLoaded0= IsLoaded; %#ok<NASGU>
IsLoaded = zeros(1,NSamples);
IsLoaded(LoadSt)  = 1;
IsLoaded(LoadEnd) = IsLoaded(LoadEnd) - 1;
IsLoaded      = [Init cumsum(IsLoaded)];
IsLoaded(end) = [];
% plausibility check
if max(IsLoaded)>1 || min(IsLoaded)<0
    error('unexpected value for Loaded');
end
IsLoaded = logical(IsLoaded);

% % optional plotting for 
% if any(Duration<MinimalLoadDuration)
%     Time = (1:length(FReactZ)) / (10*MinimalLoadDuration);
%     Tim0 = Time;Tim0(~IsLoaded0)=nan;
%     Tim1 = Time;Tim1(~IsLoaded) =nan;
%     figure;hold on; title('detection of loaded periods from force')
%     plot(Time,FReactZ);
%     plot(Tim0,FReactZ,'.-');
%     plot(Tim1,FReactZ,'linewidth',1);
%     legend('ForceZ','above threshold','loaded detected')
% end
end



% conditions tested for validation:
%           1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
% Loaded = [1 1 1 1 1 0 1 1 1 0 1 0 1 1 0 0 0 0 1]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [1 1 1 0 0 0 0 1 0 0 1 0 1 1 0 0 0 0 1]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [1 0 1 0 0 0 0 1 1 0 1 0 1 1 0 0 0 1 0]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;
% Loaded = [0 0 1 0 0 0 0 1 0 0 1 0 1 1 0 0 0 1 0]; NS=length(Loaded);no=1;MinDuration=2;Freq=1;no=1;

