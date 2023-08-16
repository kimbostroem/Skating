function [DataFilled,Prefix,LabelsClean] = doGapfillingRGB(Data,Pref)
%% Select files for gap-filling and save the filled data as well as a list of
%% filled gaps and the estimated quality.
%
% Marker-based trajectories with gaps are interpolated hierarchically.
% Gaps are interpolated by polynomials (5th order) and extrapolated (0th order) using
% tails of 30 or more samples (depending on the frequency).
% 0. Analyse the marker names on the basis of the "MarkerTable"
% 1. Remove small jumps in the markers ("spikes" in velocity) due to occlusions
% 2. Remove "small" gaps of "MaxSmallGap" or fewer samples
% 3. Interpolate rigid bodies:
%   a) Gaps with zero markers present
%   b) Gaps with one marker present
%   c) Gaps with two markers present
%   d) Gaps with more than two markers present
% 4. Interpolate relational gaps, i.e. with respect to markers on the same segment
%   a) With respect to one marker
%   b) With respect to two markers
% 5. Fill all further gaps
% 6. Plot schemas of the filled gaps by type of filling and the estimated quality
% 7. Export/return the data
%
% Needs:  MarkerTable : this is an Excel Table, listing the allowed Marker names,
%                       the segments and the rigid bodies
%
% Inputs:
%     Data        - (Qualisys struct) imported data from Qualysis type mat file
%     Pref        - (struct) various settings
% Outputs:
%     DataFilled  - (Qualisys struct) filled data from Qualysis type mat file
%     Prefix      - (char) Subject-specific prefix to the marker labels (e.g. 'S1_')
%     LabelsClean - (cell) The marker labels without the prefix
%
% Local functions: checkMarkerLabels, fill6dof0, fill6dof1, fill6dof2, fill6dof3, makeMarkerFigures,
%                  makeOverviewFigure, datenAusQualisys, veclen, saveFigure, saveCurrentFigure

%
% (c) 2020 by Movement Science, WWU University of Muenster
% Marc de Lussanet
% version 22 230207 (MdL & LK) fixed bug that caused unfilling of only the last gap; reporting of
% all-gap trajectory
% version 230418 (MdL) fixed error in fill6dof2 and fill6dof3
% version 230630 (MdL) fixed error that would make erroneous gaps in a filled marker

%% init
narginchk(1,2);
if nargin < 2, Pref = struct; end % preferences
Dim  = 2; % QTM data have time series in 2nd dimension

%% Preferences
% Override manual fills made in QTM?
if ~isfield(Pref,'DoRefillManualFills'),  Pref.DoRefillManualFills = false; end
% Do not fill again files that have already been processed?
if ~isfield(Pref,'DoSkipAlreadyFilled'),  Pref.DoSkipAlreadyFilled = true; end
if ~isfield(Pref,'DoFillSingleMarkers'),  Pref.DoFillSingleMarkers = true; end %
% Figures: 1=overview for all markers; 2 create and save fig for each marker
if ~isfield(Pref,'DoMarkerFig'),          Pref.DoMarkerFig   = false;     end
if ~isfield(Pref,'DoOverviewFig'),        Pref.DoOverviewFig = false;     end
if ~isfield(Pref,'MaxSmallGap'),          Pref.MaxSmallGap   = 0.033; end % [s] max length of small gaps
if ~isfield(Pref,'MaxJump'),              Pref.MaxJump       = 0;     end % [m] warn at jumps larger than this
if ~isfield(Pref,'Unit'),                 Pref.Unit          = 'mm';  end % [mm] unit of Data
% where to save figures of fills (do not save if empty)
if ~isfield(Pref,'File'),                 Pref.File          = '';    end % first part of name for saved figures
if ~isfield(Pref,'Path'),                 Pref.Path          = '';    end % where to save the figs
if ~isfield(Pref,'Extrapolate'),          Pref.Extrapolate   = false; end % do not extrapolate missing vals
[~,File,~]=fileparts(Pref.File);
Path=Pref.Path;

%% read the marker table (this can be def_body.xlsx)
DefBody     = loadBody('def_body');
MarkerTable = DefBody.markers;
VariableNames    = fieldnames(MarkerTable);
% get valid marker labels
ValidLabels  = checkMarkerLabels(MarkerTable);

%% organize Data parameters
Labels     = Data.Trajectories.Labeled.Labels;   % all markers present in data
Freq       = Data.FrameRate;
DataFilled = Data;
NSamples   = Data.Frames;
NMarkers   = length(Labels);
switch Pref.Unit
    case 'm',  Scale = 1;
    case 'cm', Scale = 0.01;
    case 'mm', Scale = 0.001;
    otherwise, error('undefined unit of data');
end
% parameters for polyfillgaps of single markers
Tail       = max(30,round(.15*Freq)); % 30 30 45 60 75 for 100 200 ... 500 Hz (60 ms)
Degree     = 8;
% Get the QTM prefixes of markers
[Prefix,LabelsClean] = getQTMPrefix(Labels,ValidLabels);
% if there are none, make the prefix string of length 0 (otherwise, the filling is skipped)
if isempty(Prefix)
    Prefix = {''};
else
    for i=1:length(Prefix)
        Prefix{i} = [Prefix{i} '_'];
    end
end
Nprefix = length(Prefix); % one prefix for each subject

%% initialize
% Data type for Qualisys files: 0=missing, 1=present, 2=filled manually, 3=virtual, 4=???
IsTypeDimension = size(Data.Trajectories.Labeled.Data,2)==4;
[~,~,FileExtension] = fileparts(Data.File);
IsQualisys      = strcmp(FileExtension,'.qtm') || ...
    isfield(Data,'Meta') && isfield(Data.Meta,'Company') && ...
    ~isempty(Data.Meta.Company) && contains(Data.Meta.Company,'Qualisys');
if isfield(Data.Trajectories.Labeled,'Type') && IsQualisys
    FilledType  = double(Data.Trajectories.Labeled.Type);
elseif IsTypeDimension && IsQualisys
	% filled and virtual parts have residual==0
	FilledType  = 1+double(squeeze(Data.Trajectories.Labeled.Data(:,4,:)==0));
	IsMissing   = isnan(squeeze(Data.Trajectories.Labeled.Data(:,4,:)));
	FilledType(IsMissing) = 0;
else
    % there is no information about the gap-filling status of the trajectories
	FilledType  = ones(NMarkers,NSamples);
	IsMissing   = isnan(squeeze(Data.Trajectories.Labeled.Data(:,1,:)));
	FilledType(IsMissing) = 0;
end
IsGapAll        = false(NMarkers,NSamples);
% The structure GapAll gathers all information to the gaps of the current marker
GapAll          = struct;
GapAll.id       = zeros(1,NMarkers);
GapAll.isRigid  = zeros(1,NMarkers);

%% Get Information from the marker table
% Select the IDs of the rigid bodies in the MarkerTable
RigidBodies = unique(MarkerTable.id_Rigid);
RigidBodies(isnan(RigidBodies)) = [];
% Read from the marker table which markers are located on the same segment
SegCols          = find(contains(VariableNames,'id_Segment'));
GapAll.isSegment = zeros(length(SegCols),NMarkers);
SegmentTable     = [];
for i=1:length(SegCols)
    Columns = SegCols(i);
    VariableName = VariableNames{Columns};
    SegmentTable = [SegmentTable MarkerTable.(VariableName)]; %#ok<AGROW>
end
Segments = unique(SegmentTable);
Segments(isnan(Segments)) = [];
% Read from the marker table which markers are located on the same subsegment
SubSegCols       = find(contains(VariableNames,'id_SubSegment'));
GapAll.isSubSegm = zeros(length(SubSegCols),NMarkers);
SubSegmTable     = [];
for i=1:length(SubSegCols)
    Columns = SubSegCols(i);
    VariableName = VariableNames{Columns};
    SubSegmTable = [SubSegmTable MarkerTable.(VariableName)]; %#ok<AGROW>
end
SubSegmts = unique(SubSegmTable);
SubSegmts(isnan(SubSegmts)) = [];

%% find the marker to the current label and check if it belongs to a RigidBody or Segment
for La = 1:NMarkers
    for c=3:length(VariableNames)
        % skip non-string type columns
        if ~iscell(MarkerTable.(VariableNames{c})), continue; end
        Columns = MarkerTable.(VariableNames{c});
        if any(strcmp(Columns,LabelsClean(La)))
            GapAll.row(La)         = find(strcmp(Columns,LabelsClean(La)));
            GapAll.id(La)          = MarkerTable.id(GapAll.row(La));
            GapAll.isRigid(La)     = MarkerTable.id_Rigid(   GapAll.row(La));
            for i=1:length(SubSegCols)
                Columns = SubSegCols(i);
                VariableName = VariableNames{Columns};
                GapAll.isSubSegm(Columns,La) = MarkerTable.(VariableName)(GapAll.row(La));
            end
            for i=1:length(SegCols)
                Columns = SegCols(i);
                VariableName = VariableNames{Columns};
                GapAll.isSegment(Columns,La) = MarkerTable.(VariableName)(GapAll.row(La));
            end
        end
    end
    % if marker labels are miss-spelled, they are not assigned to their segment
    if GapAll.id(La)==0
        % fprintf('doGapfillingRGB: Label "%s" is not defined. Miss-spelled?\n',Labels{La});
    end
end
% the marker table has some nan values. Replace these with zeros to avoid errors
GapAll.isRigid(isnan(GapAll.isRigid))     = 0;
GapAll.isSubSegm(isnan(GapAll.isSubSegm)) = 0;
GapAll.isSegment(isnan(GapAll.isSegment)) = 0;

%% 1 %%
%% loop the markers (labels): remove jumps and small gaps
for La = 1:NMarkers
    %% The tracked position
    Pos = datenAusQualisys(Data,Labels(La),1:3,Scale);
    % remove the manually filled gaps if desired, and if the gap does not last the entire marker
    ManFilled = FilledType(La,:)==2;
    if Pref.DoRefillManualFills && ~all(ManFilled)
        Pos(:,ManFilled) = nan;
    end

    %% Remove track jumps and check if there are any gaps
    IsGap = isnan(Pos);
    GapAll.IsEmpty(La) = false;
    SaveLargeJumpFigAs = fullfile(Path,sprintf('%s_%s_LargeJumps',File,Labels{La}));
    Filled  = Pos;
    if ~any(IsGap,'all')
        % % if there are no gaps: only remove the jumps
        % Filled = removeTrackJumps(Pos,Freq,[],Pref.MaxJump,SaveLargeJumpFigAs);
    elseif all(IsGap,'all')
        warning('marker "%s" is never present',Labels{La});
        GapAll.IsEmpty(La) = true;
    else
        % Filled = Pos;
        % % fill the gaps, remove the jumps and fill the small gaps again
        % Filled = polyfillgaps(Pos,Freq,[],[],[],[],[],[],Dim);
        % Filled = removeTrackJumps(Filled,Freq,[],Pref.MaxJump,SaveLargeJumpFigAs);
        % % remove the filled jumps again
        % Filled(IsGap) = nan;
        % fill the small gaps
        [Filled,GapSE] = polyfillgaps(Filled,Freq,[],[],Pref.MaxSmallGap,[],[],[],Dim);

        %% analyse the quality of the filling
        IsGap = false(1,NSamples);
        if ~isempty(GapSE)
            NGaps = length(GapSE(:,1));
        else
            NGaps = 0;
        end
        for i=1:NGaps
            Rng = GapSE(i,1):GapSE(i,2);
            IsGap( Rng) = true;
            FilledType(La,Rng) = 10; % autofill with respect to no other markers
        end
        IsGapAll(La,:)   = IsGap;
    end

    %% write current marker back to data structure
    DataFilled.Trajectories.Labeled.Data(La,1:3,:) = Filled/Scale;
end

%% 2 %%
%% go through the rigid bodies
IsGapAll    = squeeze(DataFilled.Trajectories.Labeled.Data(:,1,:));
IsGapAll    = isnan(IsGapAll); % find all nan values in all markers
% Check if there are further rigid bodies defined in the data
if isfield(Data,'RigidBodies')
    % create  new, unique IDs for these rigid bodies
    RigidBodies = [RigidBodies; ((1:length(Data.RigidBodies.Name))+length(Labels))'];
    for rgb = 1:length(Data.RigidBodies.Name)
        % ASSUMPTION: the marker labels contain the name of the rigid body
        IsRGB = contains(Labels,Data.RigidBodies.Name{rgb});
        IsRGB = IsRGB & GapAll.id==0; % Only if the label is not listed in MarkerTable
        GapAll.isRigid(IsRGB) = rgb + length(Labels); % assign the current ID
    end
end
% loop all rigid bodies
for r=1:length(RigidBodies)
    % loop all prefixes
    for ww = 1:Nprefix
        IsOnRGB = false(1,NMarkers);
        %% create list of markers that belong to the current RGBody
        for La=1:NMarkers
            if ~GapAll.IsEmpty(La) % the marker is not completely empty
                if GapAll.isRigid(La)==RigidBodies(r) % if on the current RGB
                    % check the prefix
                    HasCurrentPrefix = strncmp(Labels{La}, Prefix{ww}, length(Prefix{ww}));
                    NoPrefixPresent  = length(Prefix)==1 && isempty(Prefix{ww});
                    if NoPrefixPresent || HasCurrentPrefix
                        IsOnRGB(La) = true;
                    end
                end
            end
        end
        %% only if there are any gaps in any of these markers
        if ~any(IsGapAll(IsOnRGB, :),'all'), continue; end
        % ... and if enough markers of the rgb are present
        if sum(IsOnRGB)<3,                   continue; end

        %% fill the rigid body
        FilledRGB = datenAusQualisys(DataFilled,Labels(IsOnRGB),1:3,Scale);
        [FilledRGB,Gaps0] = fill6dof0(FilledRGB,Freq);
        [FilledRGB,Gaps1] = fill6dof1(FilledRGB);
        [FilledRGB,Gaps2] = fill6dof2(FilledRGB,Freq);
        [FilledRGB,Gaps3] = fill6dof3(FilledRGB);
        % Gaps0 can be 20 (good), 24 (dubious) or 25 (probably bad)
        Gaps0 = max(Gaps0, FilledType(IsOnRGB,:));
        FilledTypeRGB = max(Gaps0, 21*Gaps1 + 22*Gaps2 + 23*Gaps3);

        %% write markers of current rigid body back to data structure
        DataFilled.Trajectories.Labeled.Data(IsOnRGB,1:3,:) = FilledRGB/Scale;
        FilledType(IsOnRGB,:) = FilledTypeRGB;

        if Pref.DoOverviewFig
            if any(FilledTypeRGB==0,'all')
               Msg = whereAreTheGaps(FilledTypeRGB==0);
               fprintf('not all gaps were filled in RGB #%d (filled later): %s\n', RigidBodies(r),Msg);
            end
            Title = sprintf('%s Rigid body no %d',File,RigidBodies(r));
            makeOverviewFigure(FilledTypeRGB,Labels(IsOnRGB),Path,Title);
        end
    end
end

%% 3a %%
%% Relational filling: go through the SUB-segments
IsGapAll = squeeze(DataFilled.Trajectories.Labeled.Data(:,1,:));
IsGapAll = isnan(IsGapAll);
for ww = 1:Nprefix
    for r=1:length(SubSegmts) % r=2 r=3
        %% create list of markers that belong to the current SubSegment
        IsOnSeg = false(1,NMarkers);
        for La=1:NMarkers
            if ~GapAll.IsEmpty(La) % the marker is not completely empty
                if any(GapAll.isSubSegm(:,La)==SubSegmts(r)) % if on the current segment
                    % check the prefix
                    HasCurrentPrefix = strncmp(Labels{La}, Prefix{ww}, length(Prefix{ww}));
                    if HasCurrentPrefix
                        IsOnSeg(La) = true;
                    end
                end
            end
        end
        %% only if there are any gaps in any of these markers
        if ~any(IsGapAll(IsOnSeg, :),'all'), continue; end
        if sum(IsOnSeg)<2, continue; end

        %% relational filling of the gaps
        FilledSeg         = datenAusQualisys(DataFilled,Labels(IsOnSeg),1:3, Scale);
        [FilledSeg,Gaps0] = fill6dof0(FilledSeg,Freq);
        [FilledSeg,Gaps1] = fill6dof1(FilledSeg);
        [FilledSeg,Gaps2] = fill6dof2(FilledSeg,Freq);
        [FilledSeg,Gaps3] = fill6dof3(FilledSeg);

        FilledTypeRGB = FilledType(IsOnSeg,:);
        FilledTypeRGB(Gaps0>0) = 30;
        FilledTypeRGB(Gaps1) = 31;
        FilledTypeRGB(Gaps2) = 32;
        FilledTypeRGB(Gaps3) = 33;

        %% write markers of current rigid body back to data structure
        DataFilled.Trajectories.Labeled.Data(IsOnSeg,1:3,:) = FilledSeg/Scale;
        FilledType(IsOnSeg,:) = FilledTypeRGB;

        if Pref.DoOverviewFig
            if any(FilledTypeRGB==0,'all')
                Msg = whereAreTheGaps(FilledTypeRGB==0);
                fprintf('Not all gaps were filled in SubSegment #%d. %s\n',SubSegmts(r),Msg);
            end
            Title = sprintf('%s Subsegment no %d',File,SubSegmts(r));
            makeOverviewFigure(FilledTypeRGB,Labels(IsOnSeg),Path,Title);
        end
    end
end

%% 3b %%
%% Relational filling: go through the segments
IsGapAll = squeeze(DataFilled.Trajectories.Labeled.Data(:,1,:));
IsGapAll = isnan(IsGapAll);
for ww = 1:Nprefix
    for r=1:length(Segments) % r=2 r=3
        %% create list of markers that belong to the current Segment
        IsOnSeg = false(1,NMarkers);
        for La=1:NMarkers
            if ~GapAll.IsEmpty(La) % the marker is not completely empty
                if any(GapAll.isSegment(:,La)==Segments(r)) % if on the current segment
                    % check the prefix
                    HasCurrentPrefix = strncmp(Labels{La}, Prefix{ww}, length(Prefix{ww}));
                    NoPrefixPresent  = length(Prefix)==1 && isempty(Prefix{ww});
                    if NoPrefixPresent || HasCurrentPrefix
                        IsOnSeg(La) = true;
                    end
                end
            end
        end
        %% only if there are any gaps in any of these markers
        if ~any(IsGapAll(IsOnSeg, :),'all'), continue; end
        if sum(IsOnSeg)<2, continue; end

        %% relational filling of the gaps
        Segment           = datenAusQualisys(DataFilled,Labels(IsOnSeg),1:3, Scale);
        [FilledSeg,Gaps1] = fill6dof1(Segment);
        [FilledSeg,Gaps2] = fill6dof2(FilledSeg,Freq);
        [FilledSeg,Gaps3] = fill6dof3(FilledSeg);

        FilledTypeRGB = FilledType(IsOnSeg,:);
        FilledTypeRGB(Gaps1) = 31;
        FilledTypeRGB(Gaps2) = 32;
        FilledTypeRGB(Gaps3) = 33;

        %% write markers of current rigid body back to data structure
        DataFilled.Trajectories.Labeled.Data(IsOnSeg,1:3,:) = FilledSeg/Scale;
        FilledType(IsOnSeg,:) = FilledTypeRGB;

        if Pref.DoOverviewFig
            % report non-filled gaps
            if any(FilledTypeRGB==0,'all')
                % find the fist and last samples that are not gaps
                DataStart = find(all(FilledTypeRGB,1),1,'first');
                DataEnd   = find(all(FilledTypeRGB,1),1,'last');
                % do not warn for non-filled leading and trailing gaps
                if any(FilledTypeRGB(:,DataStart:DataEnd)==0,'all')
                    Msg = whereAreTheGaps(FilledTypeRGB==0);
                    fprintf('Not all gaps were filled in Segment #%d. %s\n',Segments(r),Msg);
                end
            end
            Title = sprintf('%s Segment no %d',File,Segments(r));
            makeOverviewFigure(FilledTypeRGB,Labels(IsOnSeg),Path,Title);
        end
    end
end

%% 4 %%
%% loop the markers (labels) and fill all remaining gaps, and estimate the quality
if Pref.DoFillSingleMarkers
    for La = 1:NMarkers 
        % The tracked position
        Pos = datenAusQualisys(DataFilled,Labels(La),1:3, Scale);
        IsGap = isnan(Pos);
        if any(IsGap,'all') && ~GapAll.IsEmpty(La)
            % GapSE is a matrix of NGaps x 2, giving first and last sample of each gap
            [Filled,GapSE,EstErrors] = polyfillgaps(Pos,Freq,Tail,Degree,[],[],[],[],Dim);
            % classify the filled gaps
            GapLn  = GapSE(:,2)-GapSE(:,1)+1; % array with the no of samples in each gap
            Unfill = []; % logical array of the samples that should be unfilled again
            for g=1:length(GapSE(:,1))
                Rng = GapSE(g,1):GapSE(g,2);
                IsGap( Rng) = true;
                IsExtrapolated = Rng(1) == 1 || Rng(end) == NSamples;
                if Pref.Extrapolate==false && IsExtrapolated
                    % the extrapolated range is to be omitted again (after removing the jumps)
                    Unfill = [Unfill Rng]; %#ok<AGROW> 
                elseif EstErrors(g)<.01 && (GapLn(g)/Freq) <= .5
                    FilledType(La,Rng) = 40; % 'good';
                elseif EstErrors(g)<.02 && (GapLn(g)/Freq) <= 10.
                    FilledType(La,Rng) = 41; % 'dubious';'dubious';
                else
                    FilledType(La,Rng) = 42; % 'bad';
                end
            end
            % Remove track jumps 
            Filled = removeTrackJumps(Filled,Freq,[],Pref.MaxJump,SaveLargeJumpFigAs);
            Filled(:,Unfill) = nan;
            % write current marker back to data structure
            DataFilled.Trajectories.Labeled.Data(La,1:3,:) = Filled/Scale;
        end
    end
else
    %% Remove track jumps
    for La = 1:NMarkers
        Filled = datenAusQualisys(DataFilled,Labels(La),1:3, Scale);
        IsGap = isnan(Filled);
        GapAll.IsEmpty(La) = false;
        SaveLargeJumpFigAs = fullfile(Path,sprintf('%s_%s_LargeJumps',File,Labels{La}));
        if ~any(IsGap,'all')
            % if there are no gaps: only remove the jumps
            Filled = removeTrackJumps(Filled,Freq,[],Pref.MaxJump,SaveLargeJumpFigAs);
        elseif all(IsGap,'all')
            warning('marker "%s" is never present',Labels{La});
            GapAll.IsEmpty(La) = true;
        else
            % fill the gaps, remove the jumps and fill the small gaps again
            Filled = polyfillgaps(Pos,Freq,[],[],[],Pref.Extrapolate,[],[],Dim);
            Filled = removeTrackJumps(Filled,Freq,[],Pref.MaxJump,SaveLargeJumpFigAs);
            % remove the filled jumps again
            Filled(IsGap) = nan;
        end
        % write current marker back to data structure
        DataFilled.Trajectories.Labeled.Data(La,1:3,:) = Filled/Scale;
    end
end
DataFilled.Trajectories.Labeled.Type = uint8(FilledType);

if all(FilledType<10,'all')
    fprintf('... No gaps were filled\n\n');
elseif any(FilledType == 0,'all')
    fprintf('There are remaining gaps of length (samples):\nTrajectory Start  Mid  End Label\n');
    [~,St,Md,En,HasGap] = whereAreTheGaps(FilledType==0);
    for i=1:length(HasGap)
        if HasGap(i)
            fprintf('%10d%6d%5d%5d %s\n',i,St(i),Md(i),En(i),Labels{i});
        end
    end
    fprintf('... Finished gap filling\n\n')
else
    fprintf('... All gaps were filled.\n\n');
end
makeMarkerFigures(DataFilled,Pref.DoMarkerFig,Path,File,GapAll,Scale);

end







%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================
%% specific functions
%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================







%% ========================================================================


function [Msg,StartSamples,EndSamples,MidSamples,HasGap] = whereAreTheGaps(IsGap)
%% count how many samples at beginning, end and middle are missing

% Init
NTraj = size(IsGap,1);
Msg='';
StartSamples = zeros(NTraj,1);
EndSamples = zeros(NTraj,1);
MidSamples = zeros(NTraj,1);

% check if there are any gaps
HasGap = any(IsGap,2);
if ~any(HasGap)
    return; % return if there are no gaps
else
    % construct a message string to report the gaps in each of the trajectories
    Msg = sprintf('Gap sizes in %d of %d traj: ',sum(HasGap),length(HasGap));
    for i=1:NTraj
        if HasGap(i)
            if all(IsGap(i,:))
                Msg = sprintf('%s Traj %d: entirely empty.   ',...
                    Msg,i);
            else            
                StartSamples(i) = find(~IsGap(i,:),1)-1;
                EndSamples(i) = length(IsGap(i,:)) - find(~IsGap(i,:),1,'last');
                MidSamples(i) = sum(IsGap(i,:))-StartSamples(i)-EndSamples(i);
                Msg = sprintf('%s Traj %d: %d, %d, %d (Start, Mid, End).   ',...
                    Msg,i,StartSamples(i),EndSamples(i),MidSamples(i));
            end
        end
    end
end
end


%% ========================================================================


function [UniqueLabels] = checkMarkerLabels(MarkerTable)
%% check that each label is unique for one marker definition
% gather unique marker label names
ColNames     = fieldnames(MarkerTable);
UniqueLabels = {};
for i=2:length(ColNames)
    if ~iscell(MarkerTable.(ColNames{i})), continue; end % skip non-string type columns
    UniqueLabels = [UniqueLabels; MarkerTable.(ColNames{i})]; %#ok<AGROW>
end
UniqueLabels = unique(UniqueLabels); UniqueLabels(1) = [];
end


%% ========================================================================


% function [Prefix,LabelsClean] = getQTMPrefix(Labels,LabelNames,MarkerPrefix)
% %
% % Version 2 (210601). Optional parameter MarkerPrefix
%
% if nargin<3 || strcmp(MarkerPrefix, '*')
%    MarkerPrefix = '';
% end
%
% Prefixes    = cell(size(Labels));
% LabelsClean = Labels;
% % The prefix is set before the label and separated by an underscore '_'
% for i=1:length(Labels)
% 	if ~isempty(strfind(Labels{i},'_')) % if there is at least one underscore
% 		[Pr,Label] = strtok(Labels{i},'_'); % separate by the underscore
% 		if any(strcmp(LabelNames,Label(2:end))) % the rest (without the '_') must be a valid label name
% 			Prefixes{i} = Pr;
% 		end
% 	end
% end
%
% % If a specific MarkerPrefix is required, then only remove this one
% if isempty(MarkerPrefix)
%    RemovedPrefixes = Prefixes;
% else
%    RemovedPrefixes = {MarkerPrefix};
% end
%
% % Get the prefix and check that there is just one prefix (Qualisys allows
% % several subjects in a single recording ut this is currently not supported)
% Prefix = {}; Nprx = 0;
% for i=1:length(Labels)
% 	if ~isempty(RemovedPrefixes{i})
% 		if isempty(Prefix)
% 			Nprx      =  Nprx+1;
% 			Prefix(1) = RemovedPrefixes(i);
% 		elseif ~any(strcmp(Prefix, RemovedPrefixes{i}))
% 			Nprx      =  Nprx+1;
% 			Prefix(Nprx) = RemovedPrefixes(i); %#ok<AGROW>
% 		end
% 		% get the labelname without the prefix
% 		[~,LabelsClean{i}] = strtok(Labels{i},'_');
% 		LabelsClean{i}(1) = '';
% 	end
% end
% end


%% ========================================================================


function [Filled,FillType] = fill6dof0(RigidBody,Freq)
%% fill gaps in rigid body while no marker is present
% Seek the marker with the shortest gap and fill that one by polyfillgaps.
% This should be followed by fill6dof1.
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 1, 210326
% Version 2, 210516 : create FillType from EstErrors of polyfillgaps
%                     This involves estimating the quality of the fill:
%                     20=good;24=dubious;25=bad
% Version 3, 210629 : bug in assignment of the quality of the filled gap

%% initialize
[Nm,Nc,Ns] = size(RigidBody); % Nmarkers, Ncoordinates, Nsamples
Filled     = RigidBody;
Dim        = 2;

%% some test gaps (third RGB of first file)
% RigidBody(3,:,10:5000)=nan;  % testgap
% RigidBody(4,:,250:350)=nan;  % testgap
% RigidBody(4,:,1990:2000)=nan;  % testgap

%% Find samples with zero markers present
IsGapAll = reshape(any(isnan(RigidBody),2),Nm,Ns);
NPresent = Nm - sum(IsGapAll,1);
Gaps0    = NPresent==0; % no marker present
FillType = zeros(Nm,Ns);
%figure;hold on; for i=1:Nm, plot(IsGapAll(i,:)+i*1.1); end

%% stop if there are no gaps to fill
if ~any(Gaps0)
    return;
end

%% analyse the gap intervals for all markers
StEn    = diff([false(Nm,1) IsGapAll false(Nm,1)],[],2);
St      = StEn(:,1:end-1)== 1; % start sample of all gaps
En      = StEn(:,2:end)  ==-1; % end sample of all gaps

%% Start and end of gaps with zero markers present
StEn0   = diff([false Gaps0 false]);
St0     = find(StEn0==1);
En0     = find(StEn0==-1)-1;
Sz      = size(St0);
if length(Sz)==1, NGaps0=1;
else,             NGaps0=Sz(2);
end

%% find for each Gap0 the marker whose gap is shortest (MFill)
ILenAll  = ones(Nm,NGaps0)*Ns;
for m=1:Nm
    for i=1:NGaps0
        if IsGapAll(m,St0(i))
            Stmi = find(St(m,1:St0(i)),  1,'last');
            Enmi = find(En(m,St0(i):end),1,'first')+St0(i)-1;
            ILenAll(m,i) = Enmi-Stmi+1; % the length of the gap for each marker that coincides with the current Gap0
        end
    end
end
[~,MFill] = min(ILenAll); % the marker that ist to be filled for each of the gaps

%% loop the gaps
FillErr   = nan(1,Ns);
GapLength = zeros(1,Ns);
FillType0 = zeros(1,Ns);
for i=1:NGaps0
    Fill = reshape(Filled(MFill(i),:,:),Nc,Ns); % the marker that is to be filled for the present gap
    [Fill,Gaps,EstErrors] = polyfillgaps(Fill,Freq,[],[],[],[],[],[],Dim); % fill all gaps of this marker
    Filled(MFill(i),:,St0(i):En0(i)) = Fill(:,St0(i):En0(i)); % copy the fill into the Rigid Body
    j = Gaps(:,1)<=St0(i) & Gaps(:,2)>=En0(i); % find which of the filled gaps applies to the current one
    if ~any(j) || sum(j)>1, error('ask marc'); end
    FillErr(St0(i):En0(i)) = EstErrors(j); % copy the estimated error of the current fill
    GapLength(St0(i):En0(i)) = (En0(i)-St0(i)+1) / Freq;
end

% classify the quality of the fills
FillType0(                FillErr > .02 | GapLength > 10)  = 25; % 'bad';
FillType0(FillType0==0 & (FillErr > .01 | GapLength >0.5)) = 24; % 'dubious';'dubious';
FillType0(FillType0==0 &  FillErr <=.01)                   = 20; % 'good';
for j=1:Nm
    FillType(j,:) = FillType0;
end

end


%% ========================================================================


function [Filled,IsGapsFilled] = fill6dof1(RigidBody)
%% fill gaps in rigid body while just one marker is present
% Seek the marker with the shortest gap and fill it.
% Vector to present marker before and after gap; get rotation about the
% normal vector and interpolate linearly.
% This function follows upon and fill6dof0 precedes fill6dof2.
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 1, 210326
% Version 2, 210629 bug in En1 (was one sample too late)
%                   bug in MPresent which could be shorter than NGaps

%% initialize
[Nm,Nc,Ns] = size(RigidBody); % Nmarkers, Ncoordinates, Nsamples
Filled     = RigidBody;

%% Find samples with one marker present
IsGapAll = reshape(any(isnan(RigidBody),2),Nm,Ns);
NPresent = Nm - sum(IsGapAll,1);
Gaps1    = NPresent==1; % one marker present
IsGapsFilled = false(size(IsGapAll));

%% stop if there are no gaps to fill
if ~any(Gaps1)
    return;
end

%% Start and end of gaps with one marker present
StEn   = diff([false Gaps1 false]);
St     = find(StEn==1);
En     = find(StEn==-1)-1;
GapLen = En-St+1;
Sz     = size(St);
if length(Sz)==1, NGaps=1;
else,             NGaps=Sz(2);
end

%% find for each Gap1 the marker whose gap is shortest (MFill)
%  and the marker that is present
ILenAll  = ones(Nm,NGaps)*Ns;
MPresent = zeros(1,NGaps);
for m=1:Nm
    for i=1:NGaps
        if ~any(IsGapAll(m,St(i):En(i)))
            MPresent(i) = m; % the marker that is present
        end
    end
end
[~,MFill] = min(ILenAll); % the marker with the shortest gap will be filled

%% If a gap contains samples with no markers present, it cannot be filled (see fill6dof0)
if any(MPresent==0)
    warning('one gap cannot be filled')
    Omit = find(MPresent==0);
    MPresent(Omit) = [];
    GapLen(Omit) = [];
    MFill(Omit) = [];
    St(Omit) = [];
    En(Omit) = [];
    NGaps = length(St);
end

%% loop the gaps
for i=1:NGaps
    Ref   = reshape(RigidBody(MPresent(i),:,:),Nc,Ns); % the reference marker
    Miss  = reshape(RigidBody(MFill(i),   :,:),Nc,Ns); % the missing marker
    Edges = [max(St(i)-1,1) min(En(i)+1,Ns)]; % the sample before and after the gap (edges)
    Vec   = Miss(:,Edges) - Ref(:,Edges); % the vector between them at the edges
    if St(i)==1,  Vec(:,1)=Vec(:,2);  end % if gap at beginning, use constant length
    if En(i)==Ns, Vec(:,2)=Vec(:,1);  end % if gap at end, use constant length
    Len   = veclen(Vec); % the length of the vector
    Norm  = cross(Vec(:,1),Vec(:,2),1); % the norm vector
    Norm  = Norm / veclen(Norm);        % scaled to unit length
    Angle = acos(dot(Vec(:,1),Vec(:,2))/(Len(1)*Len(2))); % the angle around the norm vector
    if St(i)==1 || En(i)==Ns
        TAngle = (0 : GapLen(i)+1) * 0.0; % gap at beginning or end: assume no rotation
    else
        TAngle = (0 : GapLen(i)+1) * Angle / (GapLen(i)+1); % linear interpolation
    end
    Scale = 1 + (0 : GapLen(i)+1) * (Len(2)/Len(1) - 1) / (GapLen(i)+1); % compensation for small changes of the vector
    CosPhi= cos(TAngle);
    SinPhi= sin(TAngle);
    NormxVec = cross(Norm,Vec(:,1)); % cross product
    NormDotVec = dot(Norm,Vec(:,1)); % dot product
    if any(isnan(NormxVec)),   NormxVec=0;   end % if there is no rotation, e.g., if gap at beginning or end
    if any(isnan(Norm)),       Norm=0;       end % if there is no rotation, e.g., if gap at beginning or end
    if any(isnan(NormDotVec)), NormDotVec=0; end % if there is no rotation, e.g., if gap at beginning or end
    %% Rodrigues' rotation formula (See Wikipedia)
    TVec = zeros(3,length(TAngle));
    for t=1:length(TAngle)
        TVec(:,t) = Vec(:,1) * CosPhi(t) + NormxVec * SinPhi(t) + Norm * NormDotVec * (1-CosPhi(t));
    end
    %% fill with the translation of the reference vector + the rotating vector times the length scaling compensation
    Filled(MFill(i),:,St(i):En(i)) = Ref(:,St(i):En(i)) + TVec(:,2:end-1).*Scale(2:end-1);
    % figure;plot(squeeze(Filled(MFill(i),:,St(i)-25:En(i)+25))');
end
IsGapsFilled = reshape(any(isnan(Filled),2),Nm,Ns);
IsGapsFilled = IsGapAll & ~IsGapsFilled;

end


%% ========================================================================


function [Filled,IsGapsFilled] = fill6dof2(RigidBody,Freq,IsRelational)
%% fill the gap of markers with respect to two present markers,
% in a three-point 6-dof / rigid body (previously named fill6dof)
% The normal vector to ABC before and after disppearance of C is determined. The normal vector is
% interpolated and C is then the normal to the normal vector and vector SB.
% Extrapolation with constant value
%
% RGDBody      [Markers,XYZ,Ns]
% Freq
% IsRelational (Optional) If true, then the first marker is to be filled
%              with respect to the others (at least two)
% A,B          are the reference axis
% C            is filled with respect to A,B and the pre-and post orientation of the 6DOF, using alpha
%
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 2 : 190520 - new way to compute, without the problematic quadratic solution of Version 1
% Version 3 : 190626 - various corrections
% Version 4 : 190629 - still some more
% Version 5 : 190927 - check the tails for gaps in the other markers
% Version 6 : 210321 - Freq,Gaps,EstErrors
% Version 7 : 210326 - renamed to fill6dof2
%                    - more than three markers can be provided.
%                    - loop all markers and all gaps.
%                    - logical GapsFilled true the gaps that are filled for each marker
% Version 8 : 210516 - solved a bug which made that only gaps in the first marker were filled
% Version 9 : 210629 bug in End (was one sample too late)
% Version 10: 230418 bug in End (was one sample too early) -> file 'DHB_Vis_0094.mat'

if nargin<3, IsRelational=false; end

%% initialize
[Nm,Nc,Ns] = size(RigidBody); % Nmarkers, Ncoordinates, Nsamples
Filled     = RigidBody;
% use a long tail, for the rotation of the 6DOF is interpolated and
% that is usually much slower than the position
Tail       = max(30,round(.15*Freq)); % 30 30 45 60 75 for 100 200 ... 500 Hz (60 ms)
Degree     = 5;
Dim        = 2;

%% Find samples with one marker present
IsGapAll     = reshape(any(isnan(RigidBody),2),Nm,Ns);
NPresent     = Nm - sum(IsGapAll,1);
Gaps2        = NPresent==2; % two markers present
IsGapsFilled = false(size(IsGapAll));

%% stop if there are no gaps to fill
if ~any(Gaps2) || Nm<3
    return;
end

%% analyse the gap intervals for all markers
StEn     = diff([false(Nm,1) IsGapAll false(Nm,1)],[],2);
St       = StEn(:,1:end-1)== 1;
En       = StEn(:,2:end)  ==-1;

%%
if IsRelational
    Nfill=1;
else
    Nfill=Nm;
end

%% loop the markers
for m=1:Nfill
    CGap  = reshape(RigidBody(m,:,:),Nc,Ns); % gaps of point C of triangle ABC are filled
    Start = find(St(m,:)==1);   % first sample of each gap
    End   = find(En(m,:)==1);   % last sample of each gap
    Gaps2 = NPresent==2;        % two markers present (this must be refreshed each loop)
    Gaps2 = IsGapAll(m,:) .* Gaps2; % samples during which exactly two other markers are present
    GapsC = isnan(CGap(1,:)); % all nan values in marker C
    %% loop the gaps
    a=0;b=0; % index of A and B
    for s=1:length(Start)
        if Gaps2(Start(s))==1 && Gaps2(End(s))==1 % only 2-markers-present-gaps
            Edges = [Start(s)-1 End(s)+1]; % the sample before and after
            Edges(Edges<1 ) = 1;
            Edges(Edges>Ns) = Ns;
            Gap   = Start(s) : End(s);
            % if the same A and B markers are present, as in the previous gap,
            % then the existing interpolation can be used
            if s>1 && a && b && ...
                    ~IsGapAll(a,Start(s)) && ~IsGapAll(a,End(s))   && ...
                    ~IsGapAll(a,Edges(1)) && ~IsGapAll(a,Edges(2)) && ...
                    ~IsGapAll(b,Start(s)) && ~IsGapAll(b,End(s))   && ...
                    ~IsGapAll(b,Edges(1)) && ~IsGapAll(b,Edges(2))
                Filled(m,:,Gap) = Cint(:,Gap);
                continue
            end
            % else: find the two markers, A and B that are present
            a=0;b=0; % index of A and B
            for ab=1:Nm
                if  ab~=m && ...
                        ~IsGapAll(ab,Start(s)) && ~IsGapAll(ab,End(s)) && ...
                        ~IsGapAll(ab,Edges(1)) && ~IsGapAll(ab,Edges(2))
                    if     ~a,  a=ab;
                    elseif ~b,  b=ab;
                    end
                end
            end
            if a && b
                A    = reshape(RigidBody(a,:,:),Nc,Ns);
                B    = reshape(RigidBody(b,:,:),Nc,Ns);

                %% 6DOF is a Triangle ABC, C is the one that gets missing
                AC   = CGap-A;
                AB   = B   -A;

                %% get point S (SC being the height line of ABC)
                % P is the relative location of S on AB
                Ct   = 1:Ns;
                P_t  = sum((AB .* AC),'omitnan') ./ sum((AB .* AB),'omitnan');
                % linear interpolation of the gaps (the length P should be constant)
                % (using the median of P would be more accurate, but this leads to small jumps at the edges)
                P_t(GapsC) = interp1(Ct(~GapsC),P_t(~GapsC),Ct(GapsC));
                if Start(s)==1, P_t(Start(s):End(s)) = P_t(End(s)  +1); end % extrapolation with constant value
                if End(s)==Ns,  P_t(Start(s):End(s)) = P_t(Start(s)-1); end % ... same
                S          = A + P_t .* AB; % point S (SC being the height line of ABC)
                SCvec      = CGap - S;      % vector SC
                SClen_t    = veclen(SCvec); % length of SC
                SClen_t(GapsC) = interp1(Ct(~GapsC),SClen_t(~GapsC),Ct(GapsC)); % again linear interpolation of gaps
                if Start(s)==1, SClen_t(Start(s):End(s)) = SClen_t(End(s)  +1); end % extrapolation with constant value
                if End(s)==Ns,  SClen_t(Start(s):End(s)) = SClen_t(Start(s)-1); end % ... same
                NormVecGap = cross(SCvec,AB,1);                % the norm vector on ABC
                NormVecGap = NormVecGap./veclen(NormVecGap,1); % normalized length

                %% now SC is the norm vector to vector AB and NormVec
                NormVfill  = polyfillgaps(NormVecGap,Freq,Tail,Degree,[],[],[],[],Dim);
                NormVfill  = NormVfill ./ veclen(NormVfill);     %

                %% get SC for the filled
                SCfill     = cross(AB./veclen(AB),NormVfill);
                SCfill     = SClen_t .* SCfill ./ veclen(SCfill);
                Cint       = S + SCfill;
                Filled(m,:,Gap) = Cint(:,Gap);
                %figure;hold on; plot(Cint');plot(CGap','linewidth',1);
            end
        end
    end
    %figure;hold  on; plot(squeeze(Filled(m,:,:))','linewidth',1);plot(CGap','linewidth',1);
end
IsGapsFilled = reshape(any(isnan(Filled),2),Nm,Ns);
IsGapsFilled = IsGapAll & ~IsGapsFilled;

end


%% ========================================================================


function [Filled,IsGapsFilled] = fill6dof3(RigidBody)
%% fill gaps while 3 or more markers present (see documentation)
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 1: 210406
% Version 2: 210518 bug that made that a previously complete marker became a gap
% Version 3: 210629 bug in Enm (was one sample too late)
% Version 4: 230418 bug in End (was one sample too early) -> see fill6dof3

%% initialize
[Nm,Nc,Ns] = size(RigidBody); % Nmarkers, Ncoordinates, Nsamples
Filled     = RigidBody;

%% Analyse the gaps
IsGapAll     = reshape(any(isnan(RigidBody),2),Nm,Ns); % samples with any gap for each marker
NPresent     = Nm - sum(IsGapAll,1);     % n markers present for each sample
IsGaps3      = NPresent>2 & NPresent<Nm; % there is a gap and three or more markers are present
IsGapsFilled = false(size(IsGapAll));

% stop here if there are no gaps to fill or fewer than four markers
if ~any(IsGaps3) || Nm<4
    return;
end

%% analyse the gap intervals for all markers
IsGap3All = IsGapAll;
IsGap3All(:,~IsGaps3) = false;
StEn     = diff([false(Nm,1) IsGap3All false(Nm,1)],[],2);
St       = StEn(:,1:end-1)== 1;
En       = StEn(:,2:end)  ==-1;

% seek markers that are present throughout (i.e. maximum gap of zero)
% if there ar no three markers with zero gaps, select the ones that have the shortest gaps)
MaxGaplen= maxGaplength(IsGapAll,true);

for m=1:Nm
    if MaxGaplen(m)~=0 % if the marker contains gap(s)
        Start = find(St(m,:)==1); % beginning of the gaps with three others present
        End   = find(En(m,:)==1); % end of the gaps with three others present
        %% loop the gaps of the current marker
        for gg=1:length(Start)
            %% get three reference markers (select by shortest gaps)
            MaxGapSel = MaxGaplen;
            for i=1:Nm % exclude markers that have a gap inside the current gap
                if any(IsGapAll(i,Start(gg):End(gg)))
                    MaxGapSel(m) = Ns;
                end
            end
            MaxGapSel(m) = Ns+1; % exclude the current marker from search
            [~,Idx] = sort(MaxGapSel); % sort by max gaplength
            Mfill = zeros(4,Nc,Ns);
            %% try a fill for each combination of the three reference markers,
            % to get the best fill quality
            for i=1:3
                a = Idx(1+mod(i-1,3));
                b = Idx(1+mod(i+0,3));
                c = Idx(1+mod(i+1,3));

                MGap = reshape(RigidBody(m,:,:),Nc,Ns); % this is a marker that is absent in the current interval i
                A    = reshape(RigidBody(a,:,:),Nc,Ns);
                B    = reshape(RigidBody(b,:,:),Nc,Ns);
                C    = reshape(RigidBody(c,:,:),Nc,Ns);

                %% 6DOF is a Triangle ABC, M is the one that gets missing
                AB   = B-A;
                AC   = C-A;
                AM   = MGap-A;

                %% get point S (SC being the height line of ABC); Srel is the relative location of S on AB
                Srel   = sum((AB .* AC)) ./ sum((AB .* AB));
                S      = A + Srel .* AB; % point S (SC being the height line of ABC)
                SC     = C - S;          % vector SC
                N_ABC  = cross( SC,AB,1);         % the norm vector on ABC

                %% get point U (UM being the height line of ABM); Urel is the relative location of T on AB
                Urel   = sum((AB .* AM)) ./ sum((AB .* AB)); % U should be constant...
                Urel   = median(Urel,'omitnan');
                U      = A + Urel .* AB; % point U (UM being the height line of ABM)
                UM     = MGap - U;      % vector UM

                %% V with VM the height line of [US]CM, the triangle projected along Vector AB
                Vrel   = sum((SC .* UM)) ./ sum((SC .* SC)); % V should be constant...
                Vrel   = median(Vrel,'omitnan');
                UV     = Vrel .* SC;
                V      = U + UV; % point U (UM being the height line of projected UCM)

                VM     = MGap - V;
                VMrel  = median(veclen(VM)./veclen(N_ABC),2,'omitnan');
                VMsign = sign(median(median(VM.*N_ABC,2,'omitnan'),'omitnan'));

                Mfill(i,:,:) = A  +  Urel .* AB +  Vrel .* SC + VMsign * VMrel .* N_ABC;

            end

            %% select the best fit from all three, including the weighed mean
            Error = zeros(1,4);
            for i=1:3
                Error(i) = mean(veclen(squeeze(Mfill(i,:,:))-MGap),'omitnan');
            end
            Weight = 1./(Error(1:3)*sum(1./Error(1:3))); % sum(Weight)
            for i=1:3
                Mfill(4,:,:) = Mfill(4,:,:) + Mfill(i,:,:)*Weight(i);
            end
            Error(4) = mean(veclen(squeeze(Mfill(4,:,:))-MGap),'omitnan');
            [~,Sel] = min(Error);
            % assign the best fill to the current marker and gap
            Filled(m,:,Start(gg):End(gg)) = Mfill(Sel,:,Start(gg):End(gg));
        end
    end
end
IsGapsFilled = reshape(any(isnan(Filled),2),Nm,Ns);
IsGapsFilled = IsGapAll & ~IsGapsFilled;

end


%% ========================================================================


function	makeMarkerFigures(DataFilled,DoMarkerFig,ResultDir,File,GapAll,Scale)

if ~DoMarkerFig
    return; 
end

[~,Freq,Labels,Ns] = datenAusQualisys(DataFilled);
Time       = (0:Ns-1)/Freq; % time
Nm         = length(Labels);
M2CM       = 100; % cm
ScreenSize = get(0,'ScreenSize');Width = 560; Height = 420;
[Colors,Lgnd] = farbkodeErstellen();
% loop all markers
for La = 1:Nm
    NameStr    = [File '_' Labels{La}];
    Pos        = datenAusQualisys(DataFilled,Labels(La),1:3, Scale) * M2CM;
    FilledType = double(DataFilled.Trajectories.Labeled.Type(La,:));
    % recode to avoid type value of zero
    FilledType(FilledType==1)=-1; % present is not interesting
    FilledType(FilledType==0)=1;  % set missing to 1
    Types = unique(FilledType);
    % do not plot some types
    Types(Types==-1)=[]; % do not plot present
    if isempty(Types), continue; end

    Pos2= [ScreenSize(3)/2,1,Width,Height];
    Fmn = mean(Pos,2,'omitnan');
    Fg = figure('Position',Pos2); hold on;
    set(Fg,'PaperSize',[Width/72,Height/72]);
    plot(Time,Pos'   -Fmn');%,'linewidth',1);


for Tp = Types' %Tp=23
    Fills  = Pos;
    Fills(FilledType~=Tp)=nan;
    plot(Time,Fills' -Fmn','.','color',Colors(Tp,:));
end
[~,ObjH]=legend(Lgnd([1 Types]),'location','eastoutside');
%// set size as desired
ObjH = findobj(ObjH, 'type', 'line'); %// objects of legend of type line
set(ObjH, 'Markersize', 12); %// set marker size as desired



%     Red    = Fills;
%     Orange = Fills;
%     Other  = Fills;
%     Fills(:, FilledType~=40)=nan;
%     Orange(:,FilledType~=41 & FilledType~=24)=nan;
%     Red(:,   FilledType~=42 & FilledType~=25)=nan;
%     NotOther = (FilledType==1)  | (FilledType==24) | (FilledType==25) | ...
%         (FilledType==40) | (FilledType==41) | (FilledType==42);
%     Other(:, NotOther) =nan;
% 
%     if (any(~isnan(Red(  1,:))) || any(~isnan(Orange(1,:))) || ...
%    		any(~isnan(Fills(1,:))) || any(~isnan(Other( 1,:))) )
%         Pos2= [ScreenSize(3)/2,1,Width,Height];
%         Fmn = mean(Pos,2,'omitnan');
%         Fg = figure('Position',Pos2); hold on;
%         set(Fg,'PaperSize',[Width/72,Height/72]);
%         plot(Time,Pos'   -Fmn');%,'linewidth',1);
%         plot(Time,Fills' -Fmn','g.');
%         plot(Time,Orange'-Fmn','y.');
%         plot(Time,Red'   -Fmn','r.');
%         plot(Time,Other' -Fmn','m.');
        NameInPlot = strrep(['File "' File '", Label "' Labels{La} '"'],'_','\_');
        SegString = ' ';
        if ~isnan(GapAll.isRigid(La)) && GapAll.isRigid(La)>0
            SegString = sprintf(' RGB: %d ',GapAll.isRigid(La));
        end
        Segments = GapAll.isSegment(:,La);
        Segments(isnan(Segments) | Segments==0) = [];
        if ~isempty(Segments)
            SegString = sprintf('%s. Seg: ', SegString);
            for i=1:length(Segments)
                SegString = sprintf('%s %d',SegString,Segments(i));
            end
        end
        SubSegments = GapAll.isSubSegm(:,La);
        SubSegments(isnan(SubSegments) | SubSegments==0) = [];
        if ~isempty(SubSegments)
            SegString = sprintf('%s. SubSeg: ', SegString);
            for i=1:length(SubSegments)
                SegString = sprintf('%s %d',SegString,SubSegments(i));
            end
        end
        title('x,y,z pos, jumps removed and filled.',[NameInPlot SegString]);
        xlabel('time (s)');
        ylabel('centered position (cm)')
        if ~isempty(NameStr)
            saveCurrentFigure(Fg,NameStr,ResultDir);
            close(Fg);
        end
%     end
end
end

%% ========================================================================


function	makeOverviewFigure(FilledType,Labels,ResultDir,File)

Ns         = length(FilledType);
Nm         = length(Labels);
ScreenSize = get(0,'ScreenSize');
%% figure
% recode to avoid type value of zero
FilledType(FilledType==1)=-1; % present is not interesting
FilledType(FilledType==0)=1;  % set missing to 1
Types = unique(FilledType);
% do not plot some types
Types(Types==-1)=[]; % do not plot present
if isempty(Types), return; end

% setup the plot
Width = round(min(400 + Ns/10, ScreenSize(3))); % 560;
Height= round(min(150 + 10*Nm, ScreenSize(4)));
Pos2  = [ScreenSize(3)/2,1,Width,Height];
Fg    = figure('Position',Pos2); hold on;
%set(Fg,'PaperSize',[200+Width/1500,Height/72]);
NameInPlot = strrep(['File: "' File '"'],'_','\_');
title('Gaps filled in each of the markers, by kind of filling',NameInPlot);
[Colors,Lgnd] = farbkodeErstellen();

Markers = ones(Nm,Ns).*(1:Nm)';
for Tp = Types' %Tp=23
    ThisType = ones(Nm,Ns).*(1:Ns);
    for La = 1:Nm
        ThisType(FilledType~=Tp)=nan;
    end
    Ms       = reshape(Markers, 1,Nm*Ns);
    ThisType = reshape(ThisType,1,Nm*Ns);
    plot(ThisType,Ms,'.','color',Colors(Tp,:));
end
[~,ObjH]=legend(Lgnd(Types),'location','eastoutside');
%// set size as desired
ObjH = findobj(ObjH, 'type', 'line'); %// objects of legend of type line
set(ObjH, 'Markersize', 12); %// set marker size as desired

%% Print the label names on the y-axis
set(gca,'yticklabel',[]) %Remove tick labels
%% Get tick positions
axis([1 Ns 1 length(Labels)]);
Ticks  = strrep(Labels,'_','\_');
% Adjust the offset based on the size of figure
xTicks = get(gca,'xtick');
HorizontalOffset = xTicks(end)/8;
ySpacing = 1:length(Labels);
for yy = 1:length(Ticks)
    text(- HorizontalOffset, ySpacing(yy), Ticks{yy})
end
if ~isempty(File)
    NameStr = [File '_AllGaps'];
    saveCurrentFigure(Fg,NameStr,ResultDir);
    close(Fg);
end
end


%% ========================================================================


function [Colors, Lgnd] = farbkodeErstellen()
NTypes = 42;
Farben = parula(15);
Colors = ones(NTypes,3);
Colors( 1,:)=Farben( 1,:); % 'measured'
Colors( 2,:)=Farben( 2,:); % 'missing'
Colors( 3,:)=Farben( 3,:); % 'manually
Colors( 4,:)=Farben( 4,:); % 'virtual
Colors( 5,:)=Farben(12,:); % 'virtual
Colors(10,:)=Farben( 5,:); % 'small filled gaps
Colors(20,:)=Farben( 6,:); % 'rigid body 0
Colors(21,:)=Farben( 7,:); % 'rigid body 1
Colors(22,:)=Farben( 8,:); % 'rigid body 2
Colors(23,:)=Farben( 9,:); % 'rigid body 3
Colors(24,:)=[1 1   0.6];  % 'rigid body 0 (dubious)
Colors(25,:)=[1 0.6 0.6];  % 'rigid body 0 (bad)
% Gaps0 can be 20 (good), 24 (dubious) or 25 (probably bad)
Colors(30,:)=Farben(15,:); % 'relational 1
Colors(31,:)=Farben(10,:); % 'relational 1
Colors(32,:)=Farben(11,:); % 'relational 2
Colors(33,:)=Farben(13,:); % 'relational 1 refilled
Colors(34,:)=Farben(14,:); % 'relational 2 refilled
Colors(40,:)=[0 1 0]; % 'g'
Colors(41,:)=[1 1 0]; % 'y'
Colors(42,:)=[1 0 0]; % 'r'
Lgnd{1}  = ' 1. measured';
Lgnd{2}  = ' 0. missing';
Lgnd{3}  = ' 2. manually';
Lgnd{4}  = ' 3. virtual';
Lgnd{5}  = ' 4. unknown QTM type';
Lgnd{10} = '10. small';
Lgnd{20} = '20. RGB 0 present';
Lgnd{21} = '21. RGB 1 present';
Lgnd{22} = '22. RGB 2 present';
Lgnd{23} = '23. RGB 3 present';
Lgnd{24} = '24. RGB 0 Dubious';
Lgnd{25} = '25. RGB 0 =BAD=';
Lgnd{30} = '30. relational 0';
Lgnd{31} = '31. relational 1';
Lgnd{32} = '32. relational 2';
Lgnd{33} = '33. relational 3';
Lgnd{40} = '40. probably good';
Lgnd{41} = '41. Dubious';
Lgnd{42} = '42. possiblly BAD';
end


%% ========================================================================


function [Len,Idx] = maxGaplength(Data,Gapvalue)
%% return the longest gap of nan values in Data

if nargin<2, Gapvalue = nan; end
Size = size(Data);
Len = zeros(1,Size(1));
Idx = zeros(1,Size(1));
for i=1:Size(1)
    if isnan(Gapvalue)
        NoGapPos = find(~isnan([0 Data(i,:) 0]));
    else
        NoGapPos = find(([0 Data(i,:) 0])~=Gapvalue);
    end
    [~, GrpIdx] = max(diff(NoGapPos));
    Len(i) = length(Data(i,NoGapPos(GrpIdx):NoGapPos(GrpIdx+1)-2));
    if Len(i)
        Idx(i) = NoGapPos(GrpIdx);
    end
end
end






%% ========================================================================
%% ========================================================================
%% ========================================================================
%% General usage functions
%% ========================================================================
%% ========================================================================
%% ========================================================================








%% ========================================================================

function [Pos,Freq,Labels,Nsam] = datenAusQualisys( Data, MarkerName, XYZ, Scale)
%% Examples:
%   [~,Freq,Labels] = datenAusQualisys(Data)         % return only measurement frequency and list of labels
%   Pos = datenAusQualisys(Data,{'Lab1','Lab2'},1:3) % return 3D coordinates of 2 Markers (2 x 3 x nSamp)
%   Pos = datenAusQualisys(Data,'',0)                % return 4D coordinates of all Markers (nMark x 4 x nSamp)

%% this function gets the data structure (Dat) from mat exported from QTM (Qualisys track manager)
%% Output:
%   Pos    (optional) matrix of [nLabels , nCoordinates(x,y,z,res) , nSamples]  , in mm
%   Freq   (optional) Measurement frequency
%   Labels (optional) cell list of all labels
%% Input:
%   Dat        Structure from QTM-exported .mat file
%   MarkerName (optional) cellstring of label names (if empty, then data of all labels are returned)
%   XYZ        (optional) which coordinates (x=1,y=2,z=3,rms=4; if empty, then all are returned)
% - Singleton dimensions are removed from Pos (so if just one coordinate of one marker is requested, an array is returned)
%   Scale      (optional) scaling factor of the data (default mm2m)
%
% Marc de Lussanet, Movement Science, WWU University of Muenster
% Version 2 (4.10.2017)
%     (squeeze only if necessary) and sort the dimensions and Case-Insensitive Label-selection
% Version 3 (15.4.2019)     Label selection using strcmpi (case-insensitive and same length)
% Version 4 (26.8.2019)     return if no trajectories exist
% Version 5 (25.5.2020)     return labels that are actually read
% Version 6 (11.3.2021)     return Nsam
% Version 7 (19.7.2022)     optional scale

narginchk(1,4);
if nargin<3 || isempty(XYZ),        XYZ        = 1:4;    end
if nargin<4 || isempty(Scale),      Scale      = 0.001;  end % default XYZ

% Initialise
Repair = 1;
Pos    = [];
Freq   = nan;
Labels = {};

% return if no trajectories exist
if ~isfield(Data, 'Trajectories'),    return; end
Pos    = Data.Trajectories.Labeled.Data * Scale;
Labels = Data.Trajectories.Labeled.Labels;
Freq   = Data.FrameRate;
Nsam   = Data.Frames;
if nargin<2 || isempty(MarkerName), MarkerName = Labels; end

% if just one argument then return Labels and freq
if nargin==1,                        return; end
% if the function is called without a label list, then return everything
Nlab   = size(MarkerName,2);
Nxyz   = size(XYZ       ,2);
Pos    = NaN([Nlab Nxyz Nsam]);

% copy by label, to asign the correct labels to each other
for i = 1 : Nlab
    Index  = strcmpi(Labels, MarkerName{i}); % 190415
    if sum(Index)>0
        if sum(Index)>1
            disp(['datenAusQualisys: Label ' MarkerName{i} ' multiply defined']);
            % if a label exists more than once, take the first occurrence
            if Repair==1
                x = find(Index,1);
                Index(:) = 0;
                Index(x) = 1;
            end
        end
        Pos(i,XYZ,:) = Data.Trajectories.Labeled.Data(Index,XYZ,:) * Scale;
    end
end

% remove singleton dimensions, but not if some labels are not found
Sz = size(Pos);
if Nxyz == 1 && Nlab == 1
    Pos    = squeeze(Pos);
elseif Nlab == 1
    Pos    = reshape(Pos,Sz(2:3));
elseif Nxyz == 1
    Pos    = reshape(Pos,[Sz(1) Sz(3)]);
end
Labels    = MarkerName;  % return labels that are actually read

end


%% ========================================================================


function [Len] = veclen(Vector3D,Dim)
if nargin==1, Dim=1; end
X=1;Y=2;Z=3;
if Dim==1
    Len = sqrt(Vector3D(X,:).^2 + Vector3D(Y,:).^2 + Vector3D(Z,:).^2);
else
    Len = sqrt(Vector3D(:,X).^2 + Vector3D(:,Y).^2 + Vector3D(:,Z).^2);
end
end


%% ========================================================================
%% ========================================================================
%% ========================================================================


function saveFigure(fig,filePath,fileType,varargin)
% Kim J Bostroem. Movement Science, WWU Muenster

% if (length(varargin)==2)
% 	width = varargin{1};
% 	height = varargin{2};
% 	scrsz = get(0,'ScreenSize');
% 	set(fig,'Position',[scrsz(1),scrsz(4),width,height]);
% 	set(fig,'Position',[1,1,width,height]);
% else
% 	figPos = get(fig,'Position');
% 	width  = figPos(3);
% 	height = figPos(4);
% end
% set(fig,'PaperPositionMode','auto');
% set(fig,'PaperUnits', 'points');
% set(fig,'PaperSize', [width height]);
% set(fig,'renderer','painters');

% print resolution
% 72 : monitor resolution
% 300 : printer resolution
printResolution = 300;

switch fileType
    case 'fig'
        saveas(fig,sprintf('%s.fig',filePath));
    case 'png'
        print( fig,sprintf('%s',filePath),'-dpng',sprintf('-r%.0f',printResolution));
    case 'pdf'
        print( fig,sprintf('%s',filePath),'-dpdf',sprintf('-r%.0f',printResolution));
    case 'eps'
        print( fig,sprintf('%s',filePath),'-deps',sprintf('-r%.0f',printResolution));
end

end


function saveCurrentFigure(varargin)
% saveCurrentFigure(fig,fileName,fileDir)
% Kim J Bostroem. Movement Science, WWU Muenster

fig = gcf;
fileName = 'Figure';
fileDir = '';
if length(varargin)    ==1, fig = varargin{1};
elseif length(varargin)==2, fig = varargin{1}; fileName = varargin{2};
elseif length(varargin)>=3, fig = varargin{1}; fileName = varargin{2}; fileDir = varargin{3};
end

figPos    = get(fig,'Position');
width     = figPos(3);
height    = figPos(4);
fileTypes = {'png'};                              % fileTypes = {'fig','pdf','png'};

%set(fig,'PaperPositionMode','auto');
%set(fig,'PaperUnits', 'points');
%set(fig,'PaperSize', [width height]);
%set(fig,'renderer',  'painters');

for i=1:length(fileTypes)
    if ~exist(fileDir, 'dir') && ~isempty(fileDir) %fullfile(fileDir,fileTypes{i}), 'dir')
        mkdir(fileDir);                             %fullfile(fileDir,fileTypes{i}));
    end
    filePath = fullfile(fileDir,fileName);         % fullfile(fileDir,fileTypes{i},fileName);
    saveFigure(fig,filePath,fileTypes{i},width,height);
end

end


%% ========================================================================





