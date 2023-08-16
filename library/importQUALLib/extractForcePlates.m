function   [Stream,Msgs,Flag] = extractForcePlates(Stream, Data, DefBody, Segments, Info,Msgs,Flag)
%% Get the force plate data from the Qualisys data
% assign the forces to the force items (foot, hand, etc), on the basis of the markers 
% that are assigned to the item. Compute the COP in the local coordinates of the Segments.
% Output them to Stream.
%
% SYNTAX
%   Stream = extractForcePlates(Stream, Data, RefBody, Segments, Info)
%
% INPUT
%     Stream 		- (struct) MNData format
%     Data  		- (struct) QTM data format, containing ground forces and marker data
%     DefBody		- (struct) imported def_body table
%     Segments		- (struct) structure of Coordinate Systems of body segments
%     Info  		- (struct) parameters from preferences file
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs          - (char) notifications, warning and error messages
% 
% OUTPUT
%     Stream 		- (struct) MNData format with added measured contact forces
%     Msgs          - (char) notifications, warning and error messages
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%
% Local functions: omitTrailingNanValues
% 
% See also: importQUAL
% 
% (c) 2021 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230125 (MdL) improved variable names
% version 230222 (MdL) fixed output parameter list for getCOPfromAnalog
% version 230306 (MdL) improvements in plotting function
% version 230310 (MdL) handle Units
% version 230516 (MdL) provide force plate Type to getCOPfromAnalog
% version 230616 (MdL) changed function name addWorldObject to structuredGroundLevel
% version 230622 (MdL) Plot worldObjects with COP and 3D plotting of COP with markers
% version 230804 (MdL) crop after applying the world objects

%% init
DoInsert_IsContact       = false; % flag for adding the detected contact time samples
NSintQTM                 = Info.NSintQTM;
ThresholdCOPDist         = Info.ThresholdCOPDist;
ThresholdForceItemHeight = Info.ThresholdForceItemHeight;
QTMTimeSRint             = Info.QTMTimeSRint;
DEBUGplot                = false;
WorldObjects             = makeWorldObjectsStruct(DefBody.worldObjects);
X=1;Z=3;

% get force data form QTM structure
if isfield(Data,'Units')
    if strcmp(Data.Units,'mm')
        Scale = 0.001;
    elseif strcmp(Data.Units,'m')
        Scale = 1;
    else
        Flag = -1;
        Msgs = [Msgs, {sprintf('loadInputFile: non-defined units %s',Data.Units)}];
        return;
    end
else
    Scale = 1;
end
[ForceRaw, ~, COP0, ~, Location, Analog,Type] = kraftAusQualisys(Data, Scale, Info.ForcePlateLabels);

% if any forces have been measured
IsPresentForcesAndMarkers = ~isempty(ForceRaw) && isfield(Stream,'markers') && ~isempty(Stream.markers);
if IsPresentForcesAndMarkers
    NPlates = size(Data.Force,2);
    FFreq   = size(NPlates);
    for pl=1:NPlates
        FFreq(pl) = Data.Force(pl).Frequency;
    end
    ForcePlateNames = {Data.Force(:).ForcePlateName};
    IsForcePlate    = contains(ForcePlateNames,Info.ForcePlateLabels,'IgnoreCase',true);
    
    %% Calc ground forces for each foot separately
    %
    % first interpolate the force data on the time-vector tkin of the kinematic data
    % second rotate the COP and the Force-vector if necessary
    nForceplates  = sum(IsForcePlate);
    nSampForce    = size(ForceRaw,3);
    Forces        = cell(nForceplates,1);
    COPs          = cell(nForceplates,1);
    TimeForce     = (0:(nSampForce-1))/FFreq(1);
    
    % sometimes, the QTM data end with a few NaNs: omit these
    [ForceRaw, COP0, Analog] = omitTrailingNanValues(ForceRaw, COP0, Analog);

    % Recalculate the COP, if the Analog data have been exported, too.
    % From those, the COP can be calculated much more accurately
    [COP,~,ForceRaw,Loaded,~,Msgs,Flag] = getCOPfromAnalog(Analog,ForceRaw,[],Location,FFreq(1),[],[],[],COP0,Type);

    %% Changing height here!!
    if max([WorldObjects.objectID]) > 0
        % apply the ground level to the COPs
        for pl = 1:nForceplates
            COP(pl,:,:) = getGroundLevel(squeeze(COP(pl,:,:)), WorldObjects, squeeze(ForceRaw(pl,:,:)));
        end
        % draw the cop positions in the object-world
        Positions = COP;
        for pl = 1:nForceplates
            Positions(pl,:,~Loaded(pl,:)) = nan; 
        end
        Positions = permute(Positions,[2 3 1]);
        Positions = reshape(Positions,[3, nForceplates * nSampForce]);
        drawWorld(DefBody,Positions);
    end
    % crop the COP and force where the COP is outside the plate (and at ground level)
    [COP,ForceRaw] = cropForcesToForcePlates(COP,ForceRaw,Location);

    % Resample to internal sampling rate (internalSR)
    ForceSRi = zeros(nForceplates,3,NSintQTM);
    COP_SRi  = nan(nForceplates,3,NSintQTM);
    IsLoaded = false(nForceplates,NSintQTM);
    for pl=1:nForceplates 
        COP(pl,:,~Loaded(pl,:))      =nan;
        ForceRaw(pl,:,~Loaded(pl,:))=0;
        for cc=1:3
            % use interp1 rather then resample in order to preserve the
            % zero loadings in case of jump transients
            COPs{pl}(cc,:)  = interp1(TimeForce,squeeze(COP(      pl,cc,:)),QTMTimeSRint);
            Forces{pl}(cc,:)= interp1(TimeForce,squeeze(ForceRaw(pl,cc,:)),QTMTimeSRint);
        end
        COP_SRi(pl,:,:)  = COPs{pl};
        ForceSRi(pl,:,:) = Forces{pl};
        IsLoaded(pl,:)   = ~isnan(squeeze(COP_SRi(pl,1,:)))';
    end

    %% Select the force item id (e.g. right_hand_force=1, left_foot_force=3, etc) from the markers for each force vector
    ForceItems     = [DefBody.markers.id_forceItems; nan];
    IdOfForceItems = nan(size(ForceItems));
    for iForceItem=1:length(ForceItems)
        if ~isnan(ForceItems(iForceItem))
            IdOfForceItems(DefBody.markers.id(iForceItem)) = (ForceItems(iForceItem));
        end
    end
    MarkerIDs    = [Stream.markers.id];
    ForceItemIDs = IdOfForceItems(MarkerIDs);

    % check from the DefBody (i.e., the def_body.xlsx) which force items are present in the data
    if isempty(Segments)
        return
    end
    ValidForceItemIDs = DefBody.forces.id;
    nForceItems       = length(ValidForceItemIDs);
    PresentForceItems = false(size(ValidForceItemIDs));
    for iValid = 1:length(ValidForceItemIDs)
        if any(ForceItemIDs == ValidForceItemIDs(iValid))
            PresentForceItems(iValid) = true;
        end
    end
    PresentForceItemIDs = find(PresentForceItems)';
    ForceItemLabels     = DefBody.forces.name;

    % Analyse the distances of the markers of each defined force item (feet, hands, etc) to
    % each of the force plate vectors
    ItemMinHeight = nan(nForceItems,NSintQTM);
    MinItemDistSq = nan(nForceItems,nForceplates,NSintQTM);
    if DEBUGplot
        FH = figure;hold on; title('x-y-z plot of markers and COP');
        MarkerNames = {Stream.markers.name};
    end
    for iForceItem=PresentForceItemIDs 
        % trajectories of all markers of this item
        NoOfForceItem  = ForceItemIDs==iForceItem;
		IdsOfForceItem = MarkerIDs(NoOfForceItem)';
        Present = any(IdsOfForceItem==MarkerIDs,2);
        IdsOfForceItem(~Present) = [];
        if isempty(IdsOfForceItem) % continue if no markers represent this force item
            continue;
        end
        ForceItemMarkers = nan(length(IdsOfForceItem), 3, NSintQTM);
        for ii=1:length(IdsOfForceItem)
            no = IdsOfForceItem(ii)==MarkerIDs;
            ForceItemMarkers(ii,:,:) = Stream.markers(no).dynamic.pos.data';
        end

        ItemMinHeight(iForceItem,:)= min(squeeze(ForceItemMarkers(:,Z,:)),[],1); % the lowest position of the item
        % for each force item, loop al force plates to get the shortest distance
        for pl=1:nForceplates
            COPpl=COPs{pl};
            % compute the horizontal distance over time for each of the
            % markers in the force item
            ItemMarkersDistSquared  = nan(size(ForceItemMarkers,1),NSintQTM);
            for mm = 1:size(ForceItemMarkers,1)
                ItemMarkersDistSquared(mm,:) = sum((squeeze(ForceItemMarkers(mm,:,:))-COPpl).^2,1);
            end
            ItemMarkersDistThresholded = sqrt(ItemMarkersDistSquared);
            % ignore the values when all markers are too high (above threshold)
            IsHeightMoreThanThreshold = ItemMinHeight(iForceItem,:) - COPpl(Z,:) > ThresholdForceItemHeight;
            ItemMarkersDistThresholded(:,IsHeightMoreThanThreshold) = nan; 
            % ignore the values when all markers are too far from the COP horizontally
            IsDistanceMoreThanThreshold = ItemMarkersDistThresholded > ThresholdCOPDist;
            ItemMarkersDistThresholded(IsDistanceMoreThanThreshold) = nan;
            % get the minimum distance
            MinItemDistSq(iForceItem,pl,:) = min(ItemMarkersDistThresholded,[],1);
            if DEBUGplot %&& (iForceItem==3 || iForceItem == 4)
                % if any(sqrt(ItemMarkersDistSquared)' < ThresholdCOPDist,'all')
                %     figure; 
                %     subplot(1,2,1)
                %     plot(sqrt(ItemMarkersDistSquared)'); title(sprintf('distance from COP %d, %s',pl,strrep(ForceItemLabels{iForceItem},'_','\_')));
                %     MarkerLabels = regexprep(MarkerNames,'_',' '); 
                %     legend(MarkerLabels{sum(IdsOfForceItem==MarkerIDs,'native')})
                % 
                %     subplot(1,2,2);hold on; title(sprintf('x-y plot of markers and COP: pl %d, %s',pl,strrep(ForceItemLabels{iForceItem},'_','\_')));
                %     plot(squeeze(ForceItemMarkers(:,1,IsLoaded(pl,:)))',squeeze(ForceItemMarkers(:,2,IsLoaded(pl,:)))')
                %     plot(COPpl(1,:),COPpl(2,:),'LineWidth',2); 
                %     IsNearest = ~isnan(ItemMarkersDistThresholded(1,:));
                %     plot(squeeze(ForceItemMarkers(:,1,IsNearest))',squeeze(ForceItemMarkers(:,2,IsNearest))','.');
                %     plot(COPpl(1,IsNearest),COPpl(2,IsNearest),'LineWidth',4); axis('equal')
                %     axis('equal')
                % end
                if any(~isnan(COPpl),'all')
                    figure(FH); hold on;
                    title('Debug plot for extractForcePlates:','COP traces on all plates, and loaded markers');
                    FPLoc = squeeze(Location(pl,:,:));
                    FPLoc = [FPLoc(4,:); FPLoc]; %#ok<AGROW>
                    IsNearest = ~isnan(ItemMarkersDistThresholded(1,:));
                    plot3(...
                        squeeze(ForceItemMarkers(:,1,IsLoaded(pl,:)))',...
                        squeeze(ForceItemMarkers(:,2,IsLoaded(pl,:)))',...
                        squeeze(ForceItemMarkers(:,3,IsLoaded(pl,:)))')
                    plot3(...
                        squeeze(ForceItemMarkers(:,1,IsNearest))',...
                        squeeze(ForceItemMarkers(:,2,IsNearest))',...
                        squeeze(ForceItemMarkers(:,3,IsNearest))','.')
                    plot3(COPpl(1,:),COPpl(2,:),COPpl(3,:),'LineWidth',2); axis('equal')
                    plot3(FPLoc(:,1),FPLoc(:,2),FPLoc(:,3),'LineWidth',2); axis('equal')
                end
            end
        end
    end
    % select the item that has the marker that is nearest for each COP
    NearestItem = nan(nForceplates,NSintQTM);
    for pl=1:nForceplates
        Tmp = squeeze(MinItemDistSq(:,pl,:));
        [Val,Nearest] = min(Tmp);
        Nearest(isnan(Val)) = nan;
        NearestItem(pl,:) = Nearest;
    end
    
    % Check if there are already items present in the force stream
    if isfield(Stream,'forces')
        ForceItemOffset = length(Stream.forces);
    else
        ForceItemOffset = 0;
    end
    
    %% loop the force plates and test for each force vector, to which item it applies, if any
    IsItemOnPlate = zeros(nForceItems,NSintQTM);
    for iForceItem=PresentForceItemIDs % fi=4
        ForceItemLabel   = ForceItemLabels{iForceItem};
        ItemOnWhichPlate = NearestItem == iForceItem;
        IsItemOnPlate(iForceItem,:) = any(ItemOnWhichPlate);
        
        % Weigh the COP for all the force plates on which the foot is located.
        % (weigh each by the vertical component of the force)
        % ASSUMPTION: feet are always on different plates
        ForceIt = ForceSRi; % Force for this item initialized with the resampled (SRi) forces
        COP_It  = COP_SRi;  % COP for this item initialized with the resampled COPs
        for pl=1:nForceplates
            % set the force to zero if the item is not on force plate "pl"
            ForceIt(pl,:,~ItemOnWhichPlate(pl,:))=0;
            % weigh the COP of forceplate "pl" with the vertical component of the force for this plate
            COP_It(pl,X:Z,:) = COP_It(pl,X:Z,:) .* ForceIt(pl,Z,:);
        end
        % to return ground reaction force, we have to multiply the measured force by -1
        Force      = -squeeze(sum(ForceIt(:,X:Z,:),1))';
        WeighedCOP = squeeze(sum(COP_It,1,'omitnan') ./ sum(ForceIt(:,Z,:),1,'omitnan'));
        
        %% calc local COP from global COP
        %
        % assign global coord systems
        switch ForceItemLabel
            case 'right_foot_force', GlobalCS = Segments.rFoot;
            case 'left_foot_force',  GlobalCS = Segments.lFoot;
            case 'right_hand_force', GlobalCS = Segments.rHand;
            case 'left_hand_force',  GlobalCS = Segments.lHand;
            case 'right_knee_force', GlobalCS = Segments.rTibia;
            case 'left_knee_force',  GlobalCS = Segments.lTibia;
            case 'belt_force',       GlobalCS = Segments.thorax;
            case 'pelvis_force',     GlobalCS = Segments.pelvis;
            otherwise,               error('Unknown force %s',ForceItemLabel);
        end
        % multiply the local Coordinate System (CS) with the global COP location to get local coordinates for the COP
        COP = zeros(NSintQTM,3);
        for Sample=1:NSintQTM
            GlobalR = squeeze(GlobalCS(1:3, 1:3, Sample)); % rotation matrix for this sample
            PosSample = squeeze(GlobalCS(4, X:Z, Sample)); % joint position for this sample
            COPsample = squeeze(WeighedCOP(X:Z, Sample)); % COP position for this sample
            RelativePos = (COPsample(:) - PosSample(:)); % position of COP with respect to joint
            COP(Sample,:) = GlobalR * RelativePos; % rotate by the current orientation of the joint
        end
        COP(isnan(COP)) = 0; % if there is no contact force, set COP to the joint position
        idx = DefBody.forces.id(matches(DefBody.forces.name,ForceItemLabel));
        Stream.forces(ForceItemOffset+iForceItem).id   = idx;
        Stream.forces(ForceItemOffset+iForceItem).name = ForceItemLabel;
        Stream.forces(ForceItemOffset+iForceItem).dynamic.force.units = 'N';
        Stream.forces(ForceItemOffset+iForceItem).dynamic.force.data  = Force;
        Stream.forces(ForceItemOffset+iForceItem).dynamic.COP.units   = 'm';
        Stream.forces(ForceItemOffset+iForceItem).dynamic.COP.data    = COP;
        if DoInsert_IsContact
            Stream.forces(ForceItemOffset+iForceItem).dynamic.isContact.units = 'bool'; %#ok<*UNRCH> 
            Stream.forces(ForceItemOffset+iForceItem).dynamic.isContact.data  = IsItemOnPlate(iForceItem,:)';
        end
    end
    
end
end

%% =================================================================================================
%% =================================================================================================
%% =================================================================================================

function [ForceQual, COP0, Analog] = omitTrailingNanValues(ForceQual, COP0, Analog)
%% sometimes, the QTM data end with a few NaNs: omit these
%
% SYNTAX
% [ForceQual, COP0, Analog] = omitTrailingNanvalues(ForceQual, COP0, Analog);
%
% INPUT & OUPUT
%     ForceQual (NPlates x NDimensions x NSamples) Force data
%     COP0      (NPlates x NDimensions x NSamples) Centre of Pressure
%     Analog    (NPlates x NChannels x NSamples) Analog Channels
%
% See also: extractForcePlates
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Last revision: 27.5.2022 FIXED check for isnan on "all" dimensions

% samples for which all dimensions have a nan value
IsNanForce = squeeze(all(isnan(ForceQual),2));

% error check
ForcePlatesWithOnlyNan = all(IsNanForce,2);
if any(ForcePlatesWithOnlyNan)
    List = find(ForcePlatesWithOnlyNan);
    Format = repmat(' %d,',1,numel(List));
    error(['Force plate' Format ' has only NaN values. Please try to export the forces anew.'],List)
end

% remove samples with trailing nan values
if any(IsNanForce,"all")
    IsNanSamples = all(IsNanForce,1);
    Delete = find(~IsNanSamples,1,'Last') + 1;
    NSamp = size(ForceQual,3);
    ForceQual(:,:,Delete:end) = repmat(ForceQual(:,:,Delete-1),1,1,NSamp-Delete+1);
    COP0(     :,:,Delete:end) = repmat(COP0(     :,:,Delete-1),1,1,NSamp-Delete+1);
    Analog(   :,:,Delete:end) = repmat(Analog(   :,:,Delete-1),1,1,NSamp-Delete+1);
end
end




