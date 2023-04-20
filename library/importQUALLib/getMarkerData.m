function [Markers,Present,MarkerStream,IsHippoSet,IsFootSet,Flag,Msgs] = getMarkerData(Data,IsSkeleton,Info,UserPref,Flag,Msgs)
%% Extract the marker data
% - Perform relational gapfilling according to Rigid bodies and Segments from table "markers.xlsx"
% - Resample the Trajectories
% - Select the Labels
% - Identify the prefixes defining subjects (initial characters followed by "_" or ":")
% - Ignore non-plausible prefixes
% - Identify the "name" and "id" to each label
% - Set Qualisys "Chest" marker label to sternum_marker if clavicle_marker is already defined
% - Check that no two labels refer to the same marker
% - Report the satistics of the detected marker sets
%
% SYNTAX
%   [Markers,Present,MarkerStream] = getMarkerData(Data,MarkerTable,IsSkeleton,Info,UserPref);
%
% INPUT
%     Data          - (struct) Qualisys mat format
%     IsSkeleton    - (logical) is a skeleton present
%     Info          - (struct) information
%     UserPref      - (struct) user preferences
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs          - (char) notifications, warning and error messages
% 
% OUTPUT
%     Markers       - (struct namefields.(3,NSintQTM)) structure with all marker data (also those
%                     without valid label)
%     Present       - (struct namefields.logical) for each marker if it is present
%     MarkerStream  - (struct MNData) Markers whose label is defined in the Marker table
%     IsHippoSet    - (logical) has the hipponardo markerset been detected
%     IsFootSet     - (logical) has the extended foot markerset been detected
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs          - (char) notifications, warning and error messages
%
% Local functions: findMarkers, CheckMarkerTable, whichMarkerSet, getSortedNamesDefinedByLabels
%
% See also: importQUAL
% 
% (c) 2021 by Predimo GmbH
% version 230210 (MdL) header: difference between content of Markers and MarkerStream; workaround
% for Matlab quirk in function "all"; more robust check for selecting ambiguous clavicle marker.
% Further FIX for repair of duplicate clavicle marker; fixed name Markers.clavicula_marker (tather
% than clavicle_marker
% version 230222 (MdL) get lowest position across markers

%% init
DefBody = loadBody('def_body');
MarkerTable  = DefBody.markers;
AddFillInfo  = false; % flag for MarkerStream, to add the kind of filling that has been applied 
PrintLabels  = true;  % flag to print the list of marker labels and matched marker names
NSintQTM     = Info.NSintQTM; % no of samples
QTMTime      = Info.QTMTime;  % time array of recorded data
QTMTimeSRint = Info.QTMTimeSRint; % time array to which the data are resampled
MarkerStream = [];
Present      = [];
Markers      = [];
IsHippoSet   = false;
IsFootSet    = false;
if isfield(Data,'Force')
    ForcePlateLabels = {Data.Force.ForcePlateName};
else
    ForcePlateLabels = {};
    if UserPref.isCombine
        msg = 'No forces found with the kinematic data.';
        fprintf('%d\n',msg); 
        Msgs = [Msgs; cellstr(msg)];
        return
    end
end
if ~isfield(UserPref,'MarkerPrefix'), UserPref.MarkerPrefix = {''}; end


% Foot and/or hand, hip markers are needed to allocate force plate vectors to segments
DoGetForcesAndMarkers = any(contains(ForcePlateLabels,Info.ForcePlateLabels,'IgnoreCase',true));

% If a skeleton is requested and present and markers are not needed, the import of markers can be skipped
TryUseSkeleton = UserPref.useSkeleton;
IsDefinedSkeleton = IsSkeleton && isfield(Data,'Skeletons');
DoSkipMarkerImport = ~DoGetForcesAndMarkers && TryUseSkeleton && IsDefinedSkeleton;
if DoSkipMarkerImport
    % do not analyse markers if a Skeleton is wanted and present
    fprintf('Skipping import of marker data because skeleton data is requested and present.\n');
    return;
end

%% stop here if the data contain no trajectories
if ~isfield(Data,'Trajectories') || ~isfield(Data.Trajectories,'Labeled') || ...
        ~isfield(Data.Trajectories.Labeled,'Labels') || isempty(Data.Trajectories.Labeled.Labels)
    fprintf('No markerdata found\n');
    return;
end
fprintf('\nDetecting marker labels for marker import...\n')

% Else: markers are to be imported: initialize the structures
for QTMnr = 1:length(MarkerTable.name)
    Markers.(MarkerTable.name{QTMnr}) = nan(3,NSintQTM); 
    Present.(MarkerTable.name{QTMnr}) = false;
end

% Replace Gaps of zero values with nan (this may occur in some c3d data).
NMarkers = size(Data.Trajectories.Labeled.Data,1);
if NMarkers >1
    IsGapOfZeros = squeeze(all(Data.Trajectories.Labeled.Data(:,1:3,:)==0, 2));
else % workaround: function "all" behaves differently if there is just one marker
    IsGapOfZeros = all(squeeze(Data.Trajectories.Labeled.Data(:,1:3,:))==0);
end
for i=1:NMarkers
    Data.Trajectories.Labeled.Data(i,1:3,IsGapOfZeros(i,:)) = nan;
end

%% Fill gaps for rigid bodies and relational markers
Pref.DoFillSingleMarkers = true; % better fill at the joint level in such cases
Pref.DoRefillManualFills = true;
Pref.DoMarkerFig = false; % save a plot for each marker of which any gaps were filled
Pref.DoOverviewFig = false; % save a plot for each segment for which any gaps were filled
Pref.MaxJump = 0;
Pref.Path = 'FillFigs';
Pref.MaxSmallGap = -1; % fill no small gaps
Pref.Unit = 'm';
% Do the relational gapfilling
Data = doGapfillingRGB(Data,Pref);

% All trajectories: x-z coordinates 
Trajectories = Data.Trajectories.Labeled.Data(:,1:3,:); % (Ntraj, XYZ, NSamp)
Labels       = Data.Trajectories.Labeled.Labels;   % all markers present in data
% Omit empty trajectories
EmptyTrajectories = all(isnan(Trajectories(:,1,:)),3);
Labels        = Labels(      ~EmptyTrajectories);
Trajectories  = Trajectories(~EmptyTrajectories,:,:);
NTrajectories = size(Trajectories,1);

%% Check the MarkerTable for duplicate entries and extract a complete list of all unique valid marker names
[ValidLabels,Flag,Msgs] = CheckMarkerTable(MarkerTable,Flag,Msgs);
% If MarkerPrefix is set, only the labels containing that prefix are selected
GiveFeedback = true;
[~,LabelsWithoutPrefix] = getQTMPrefix(Labels,ValidLabels,UserPref.MarkerPrefix,GiveFeedback);

%% Check which markersets are present in the data
[IsHippoSet,IsFootSet,BestSet,Msgs] = whichMarkerSet(LabelsWithoutPrefix,Msgs);

%% Select the markers that are listed in the marker-table (from def_body)
[SelectedQTMIDs,SelectedQTMLabels,Msgs] = findMarkers(MarkerTable,LabelsWithoutPrefix,BestSet,Msgs);

%% Resample
QTMTrajResamp = zeros(NTrajectories,3,NSintQTM);
for QTMnr = 1:NTrajectories
    % resample to the internal sample rate 
    for cc=1:3
        Tmp = interp1(QTMTime,squeeze(Trajectories(QTMnr,cc,:)),QTMTimeSRint);
        QTMTrajResamp(QTMnr,cc,:) = Tmp;
    end
end

%% lookup the markerlabels from the dev_body marker table, and fill the stream
No = 0;
Label_Name_struct = struct;
MarkerStream = struct;
BottomPosition = QTMTrajResamp; % matrix for selecting the lowest body marker
for QTMnr = 1:NTrajectories % QTMnr=18
    % find the defBody index for the current marker
    LabelNo = strcmp(LabelsWithoutPrefix{QTMnr},SelectedQTMLabels);
    MarkerId = SelectedQTMIDs(LabelNo);
    if any(LabelNo) && ~all(isnan(QTMTrajResamp(QTMnr,1,:)))
        No = No + 1;
        MarkerName = MarkerTable.name{MarkerTable.id==MarkerId};
        Markers.(MarkerName) = squeeze(QTMTrajResamp(QTMnr,:,:)); % copy the trajectory
        Present.(MarkerName) = true;
        MarkerStream(No).id = MarkerId(1); %#ok<*AGROW> 
        MarkerStream(No).name = MarkerName;
        MarkerStream(No).label = LabelsWithoutPrefix{QTMnr}; % this is needed for which marker set
        MarkerStream(No).dynamic.pos.units = 'm';
        MarkerStream(No).dynamic.pos.data = squeeze(QTMTrajResamp(QTMnr,:,:))';
        % fill information of the type of filling (not implementd)
        if AddFillInfo
            MarkerStream(No).dynamic.filltype.units = 'int'; %#ok<*UNRCH> 
            MarkerStream(No).dynamic.filltype.data = Data.Trajectories.Labeled.Type(QTMnr,:)';
        end
        % Label_Name_struct: structure array of length Marker table. The label indices are thus
        % matched to the original table, which is needed for checking (see below)
        Label_Name_struct(QTMnr).LabelNo = find(LabelNo,1);
        Label_Name_struct(QTMnr).id = MarkerId(1);
        Label_Name_struct(QTMnr).Label = LabelsWithoutPrefix{QTMnr};
        Label_Name_struct(QTMnr).Name = MarkerName;
        if PrintLabels
            if QTMnr==1
                fprintf('Identified Markernames from the "markers" table:\n');
                fprintf('  nr LNo  id\t\t Label  Name\n');
            end
            fprintf('% 4d% 4d% 4d\t%14s  %s\n', QTMnr, find(LabelNo,1), MarkerId(1), ...
                LabelsWithoutPrefix{QTMnr}, MarkerName); 
        end
    else
        Label_Name_struct(QTMnr).Name = 'void';
        % ignore markers that do not belong to the body
        BottomPosition(QTMnr,:,:) = nan; 
    end
end
% Get the lowest postition across time, to estimate the floor level
[~,LowestMarker] = min(BottomPosition(:,3,:));
LowestMarker = squeeze(LowestMarker);
LowestMarkerNos = unique(LowestMarker);
Markers.bottom_position = nan(3,NSintQTM);
for i=1:length(LowestMarkerNos)
    ThisMarkerIsLowest = LowestMarker == LowestMarkerNos(i);
    Markers.bottom_position(:,ThisMarkerIsLowest) = squeeze(BottomPosition(LowestMarkerNos(i),:,ThisMarkerIsLowest));
end
Present.bottom_position = true;

%% if there are no markers with a valid label, stop here and return
if isempty(fieldnames(MarkerStream))
    MarkerStream = [];
    return;
end

%% Check for Clavicle versus Sternum marker
% the Qualisys "Chest" marker may be placed anywhere om the sternum. If there is no other sternum
% marker, it is explicitly assumed that it is the Clavicle marker. If there is a Clavicle marker,
% then assume that the "Chest" is the inferior clavicle marker
MarkerNames = getSortedNamesDefinedByLabels({MarkerStream.name});
ClavicleNo = find(matches(MarkerNames,'clavicula_marker'));
if ~isempty(ClavicleNo) && length(ClavicleNo) > 1
    Duplicates = find(matches({MarkerStream.name},MarkerNames{ClavicleNo(1)}));
    if length(Duplicates) ~= 2, error('too many chest and clavicle markers'); end
    % check which of the two markers is higher than the other (this would not work in lying or
    % upside down orientations)
    A = MarkerStream(Duplicates(1)).dynamic.pos.data(:,3);
    B = MarkerStream(Duplicates(2)).dynamic.pos.data(:,3);
    LowerMarker = 1 + (median(A-B)>0);
    UpperMarker = 3 - LowerMarker;
    MarkerStream(Duplicates(LowerMarker)).name = 'sternum_marker';
    MarkerStream(Duplicates(LowerMarker)).id   = 9;
    MarkerStream(Duplicates(UpperMarker)).name = 'clavicula_marker';
    MarkerStream(Duplicates(UpperMarker)).id   = 8;
    % Copy to the Markers structure
    Markers.sternum_marker  = MarkerStream(Duplicates(LowerMarker)).dynamic.pos.data';
    Markers.clavicula_marker = MarkerStream(Duplicates(UpperMarker)).dynamic.pos.data';
    Present.sternum_marker  = true;
end

%% Error check: a marker may have different (synonymous) labels, but only one label may refer to a marker
if ~isempty(MarkerStream)
    MarkerNames = getSortedNamesDefinedByLabels({MarkerStream.name});
    UniqueNames = unique(MarkerNames);
    if length(UniqueNames) < length(MarkerNames)
        % find the last label that still matches: this one duplicates the one that follows
        DuplicateNo = find(strcmp(MarkerNames(1:length(UniqueNames)),UniqueNames),1,'last');
        error('More than one Markers refer to the same Label (%s)',MarkerNames{DuplicateNo});
    end
end

% remove the "label" field again, since it is not defined as field for the MNDat structure
if isfield(MarkerStream,'label')
    MarkerStream = rmfield(MarkerStream,'label');
end
fprintf('... ready with label detection\n\n');
end

%% ========================================================================
%% ========================================================================

function MarkerNames = getSortedNamesDefinedByLabels(MarkerNamesCellArray)
%% get the sorted list of markers names that are defined by labels
MarkerNames = sort(MarkerNamesCellArray(~cellfun('isempty',MarkerNamesCellArray)));
end

%% ========================================================================

function [SelectedIDs,SelectedLabels,Msgs] = findMarkers(MarkerTable,LabelsWithoutPrefix,BestSet,Msgs)
%% Compare QTMLabels to the MarkerTable (from def_body.xlsx), 
% and select the labels that are listed in MarkerTable.
% 
% INPUT
%     MarkerTable         (table) "markers" table-sheet from DefBody
%     LabelsWithoutPrefix (cell array) list of all valid labels without the prefixes
%     BestSet             (char) Name of the best matching marker set
%     Msgs          - (char) notifications, warning and error messages
% OUTPUT
%     SelectedIDs    (double array) ids of the Labels that are listed in the markers table
%     SelectedLabels (cell array) names of the Labels that are listed in the markers table
%     Msgs          - (char) notifications, warning and error messages
%
% Author: Marc de Lussanet
% (c) 2021 by Predimo GmbH
% version 220926 : Use BestSet for Labels that occur for more than one marker_name
% version 221128 (MdL) on non-recognized labels be lenient on the Case

SelectedIDs    = [];
SelectedLabels = {};
NonUsedLabels  = {};

MarkerSetNames = fieldnames(MarkerTable);
IsId_Type = startsWith(MarkerSetNames,'id_');
IsNameId = matches(MarkerSetNames,{'id','anatomicID','name'});
MarkerSetNames(IsId_Type | IsNameId) = [];

% gather unique marker label names
LabelsMatrix = cell(length(MarkerTable.id),length(MarkerSetNames));
for i=1:length(MarkerSetNames)
    LabelsMatrix(:,i) = MarkerTable.(MarkerSetNames{i});
end

% compare the marker labels with the ones from the table
for i=1:length(LabelsWithoutPrefix)
    userLabel = LabelsWithoutPrefix{i};
    % - The markers HeadL and HEADL (HeadR and HEADR) of the marker sets WWU and QSkel 
    %   are at different locations, so we should check if there is a difference between 
    %   case-sensitive and case-insensitive search. 
    % - If there are differences, the case-sensitive search is preferred.
    % - Else, the case-insensitive search is preferred, to be lenient to 
    %   inconsistent case usage in marker labelling.
    idxArray           = matches(LabelsMatrix,userLabel);
    idxArrayIgnoreCase = matches(LabelsMatrix,userLabel,"IgnoreCase",true);
    idx           = any(idxArray,2);
    idxIgnoreCase = any(idxArrayIgnoreCase,2);
    if sum(idx) == 1 % there is exactly one match for this label
        SelectedIDs(end+1) = MarkerTable.id(idx);
        SelectedLabels{end+1} = userLabel;
    elseif any(idx) % there is more than one match for this label
        % select the best matching marker-label set
        BestSetCol = strcmp(MarkerSetNames,BestSet);
        % use the label reference for the best matching set
        idx = idxArray(:,BestSetCol);
        if any(idx) % the Label is present in the BestSet
            SelectedIDs(end+1) = MarkerTable.id(idx);
            SelectedLabels{end+1} = userLabel;
        else % the Label is not present in the BestSet
            Msg = sprintf('Marker label "%s" refers to different markernames but is not found in BestSet: neglect',userLabel);
            warning(Msg); %#ok<SPWRN> 
            Msgs = [Msgs; Msg];
            NonUsedLabels{end+1} = userLabel;
        end
    elseif sum(idxIgnoreCase) == 1 % there is exactly one case-independent match for this label
        SelectedIDs(end+1) = MarkerTable.id(idxIgnoreCase);
        SelectedLabels{end+1} = userLabel;
    else
        NonUsedLabels{end+1} = userLabel;
    end
end
% report which labeled markers are not used (possibly due due miss-spellings)
if ~isempty(NonUsedLabels)
    fprintf('The following labeled markers are not listed in the "markers" table:\n');
    % print the labels in 10 per line
    for i = 0 : floor(length(NonUsedLabels)/10)-1
        fprintf('\t%s\n', strjoin(NonUsedLabels(10*i + (1:10)),', '));
    end
    % the last line may be shorter
    if isempty(i)
        i=-1;
    end
    fprintf('\t%s\n\n',     strjoin(NonUsedLabels(10*(i+1) + 1: end),', '));
end
end

%% ========================================================================

function [UniqueLabels,Flag,Msgs] = CheckMarkerTable(MarkerTable,Flag,Msgs)
%% check if each label in MarkerTable is unique
% gather unique marker label names
% 
% INPUT
%     MarkerTable   - (table) "markers" table-sheet from DefBody
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs          - (char) notifications, warning and error messages
% OUTPUT
%     UniqueLabels  - (logical array)
%     Flag          - (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs          - (char) notifications, warning and error messages
%
% Version 221218 (MdL) select only marker set columns (by ignorednames list); added 'id_Segment1'

% the following columns do not contain valid marker labels
IgnoredNames = {...
    'id', 'anatomicID', 'name', ...
    'id_Rigid', ...
    'id_Segment', 'id_Segment1', 'id_Segment2', 'id_Segment3', ...
    'id_SubSegment1', 'id_SubSegment2', 'id_SubSegment3', ...
    'id_forceItems', ...
    'description', 'documentation'};
Colnames     = fieldnames(MarkerTable);
MarkerSets   = Colnames(~matches(Colnames,IgnoredNames));
% create a cell array of all labels present in the marker table
UniqueLabels = {};
for i=1:length(MarkerSets)
    UniqueLabels = [UniqueLabels; MarkerTable.(MarkerSets{i})]; 
end
UniqueLabels = unique(UniqueLabels); 
UniqueLabels(1) = [];

% Check for any marker id is higher than the length of the marker list
if any(MarkerTable.id>length(MarkerTable.id))
    Errors = find(MarkerTable.id>length(MarkerTable.id));
    for i=1:length(Errors)
        Msg = sprintf('marker id (%d) too high (%s)',MarkerTable(Errors(i)).id,MarkerTable(Errors(i)).name);
        fprintf('%s\n',Msg);
        Msgs = [Msgs; Msg];
    end
    Msg = '... there is something wrong with the listed marker ids!';
    fprintf('%s\n',Msg);
    Msgs = [Msgs; Msg];
    Flag = -1;
end
end

%% ========================================================================

function [IsHippo_Set,IsExtendedFootSet,BestSet,Msgs] = whichMarkerSet(MarkerLabels,Msgs)
%% Determine the present marker set(s) from the marker labels
%
% INPUT
%     MarkerLabels (cell)   Marker labels in Data
%     Msgs         (char)   Notifications, warning and error messages
%
% OUTPUT
%     IsHippo_Set       (logical) Are all labels present for the Hipponardo marker set
%     IsExtendedFootSet (logical) Are all labels present for the extended foot set
%     BestSet           (char)    Name of the best matching marker set
%     Msgs              (char)    Notifications, warning and error messages
%
% (c) 2020 by Predimo GmbH
% version 230120 (MdL) fixed markerset flags to singular logical
% version 230305 (MdL) make sure that BestSet is not empty

%% init
DefBody      = loadBody('def_body');
MarkerTable  = DefBody.markers;
% DefinedSets are the columns that contain Labels (except id, anatomicID, name and ones that begin 
% with id_ (i.e. segment and rigid body definitions)
DefinedSets  = fieldnames(MarkerTable);
DefinedSets(startsWith(DefinedSets,'id_')) = [];
DefinedSets(matches(DefinedSets,{'id','anatomicID','name'})) = [];
% Deterministic parameters
NNamesFound  = zeros(size(DefinedSets)); % For each set: the number of markernames detected
NLabelsFound = zeros(size(DefinedSets)); % For each set: the number of labels found in the data
NMissing     = zeros(size(DefinedSets)); % For each set: the number of missing labels
NDesired     = zeros(size(DefinedSets)); % For each set: the number of desired labels

%% the list of marker ids in the data
% PresentIds = [MarkerLabels.id];

%% communicate the result of the marker set that is present
fprintf('Markerset______ NReqrd Label  Name Missing\n')
for SetNo = 1:length(DefinedSets)
    % List of the desired Labels for the current marker set
    DesiredLabels   = MarkerTable.(DefinedSets{SetNo});
    % Logical: which of these labels are present
    PresentLabelNos = matches(DesiredLabels, MarkerLabels);
    % List of the present required labels (without possible supernumerary or miss-spelled labels)
    PresentLabels   = DesiredLabels(PresentLabelNos);
    % Desired ids of the current marker set
    DesiredIds      = MarkerTable.id;  % all ids
    DesiredIds(strcmp('',DesiredLabels)) = []; % exclude non-desired ones
    % Desired Ids that are not found
    MissingLabels   = setdiff(DesiredLabels,PresentLabels); % returns the data in A that is not in B
    MissingLabels(cellfun(@isempty,MissingLabels)) = [];
    % Statistics for all markersets
    NDesired(SetNo)     = length(DesiredIds);    % For each set: the number of desired labels
    NLabelsFound(SetNo) = length(PresentLabels); % For each set: the number of labels found in the data
    NMissing(SetNo)     = length(MissingLabels); % For each set: the number of missing labels
    NNamesFound(SetNo)  = NDesired(SetNo) - NMissing(SetNo); % : the number of markernames detected
    fprintf('%15s\t%6d%6d%6d%6d\n',DefinedSets{SetNo},NDesired(SetNo),NLabelsFound(SetNo),NNamesFound(SetNo),NMissing(SetNo));
end

%% report the result
if any(NMissing==0)
    CompleteMarkerSetsMsg = sprintf('Complete Markersets: "%s"',strjoin(DefinedSets(NMissing==0),'", "'));
else
    CompleteMarkerSetsMsg = sprintf('No complete markersets detected');
end
% Criterium: 
% - half the markers must be present (so the foot set is also detected if just one
% foot is labeled); 
% - minus 1, to be lenient;
% - at least 10 from the set to prevent that alternative labelnames are chosen
IsHalfMinOnePresent = NLabelsFound(:)>(NDesired(:)/2)-1 & NLabelsFound(:)>10;
if any(IsHalfMinOnePresent)
    MatchingMarkerSets = DefinedSets(NLabelsFound(:)>(NDesired(:)/2)-1);
    MatchingMsg = sprintf('Matching Markersets: "%s" (at least half-1 the required labels is present)', ...
        strjoin(MatchingMarkerSets,'", "'));
    % The BestSet is used if a label refers to different marker_names in different sets.
    % It is the one with most Present Labels from MatchingMarkerSets: Choose the set with
    % most present markers in order to exclude sets with alternative names. 
    [~,BestSet] = max(NLabelsFound(:) & IsHalfMinOnePresent);
    BestSet = DefinedSets{BestSet};
else
    MatchingMarkerSets = {};
    % The BestSet is used if a label refers to different marker_names in different sets.
    % It is the one with most Present Labels from MatchingMarkerSets: Choose the set with
    % most present markers in order to exclude sets with alternative names. 
    [~,BestSet] = max(NLabelsFound(:));
    MatchingMsg = sprintf('No matching markersets detected');
end
Msgs = [Msgs; CompleteMarkerSetsMsg; MatchingMsg];
fprintf('%s\n%s\n\n', CompleteMarkerSetsMsg, MatchingMsg);
IsHippo_Set       = any(contains(MatchingMarkerSets,'Hipponardo'));
IsExtendedFootSet = any(contains(MatchingMarkerSets,'GhentFoot'));

end
