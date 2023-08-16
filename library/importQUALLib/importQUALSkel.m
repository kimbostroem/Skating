function [SkelData,IsSkeleton,SimInMN,GlobalMN] = importQUALSkel(Data,UserPref,SimInMN,Info,PreBody)
%% Qualisys MAT file with Skeleton (source) to Myonardo stream (target)
% Import the sports skeleton from the QTM-exported kinematic file and write
% the angles to stream.
% Also deal with events and forces that are  contained in the QTM file.
%
% SYNTAX
% [SkelData,SimInMN,IsSkeleton] = importQUALSkel(Data,UserPref,SimInMN,Info)
%
% INPUT
%     Data     (Qualisys struct) structure containing the imported Qualisys data
%     UserPref (struct) structure containing user preferences
%     SimInMN  (MNData struct) structure containing body parameters
%     Info 	   (struct) structure containing further information (e.g. the time)
%     PreBody  (struct) User-specific parameters updated to tables of def_body (see prepareSimParam)
%
% OUTPUT
%     SkelData (MNData struct) structure containing skeleton data
%     SimInMN  (MNData struct) structure containing body parameters
%     IsSketon (logical) were skeleton data imported
%     GlobalMN (struct) orientation of the segments in global coordinates
%
% See also: importQUAL, importTheiaSkel
%
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author:  Kim Bostroem and Marc de Lussanet, Movement Science, WWU Muenster
% version 230216 (MdL) FIXED set IsSkeleton flag also if it is a body file
% version 230606 (MdL) handle marker global positions

%% init
QTMTime = Info.QTMTime;
InternalSR = Info.QTMTimeSRint; % internal sample rate
SkelData = [];
IsSkeleton = false;
GlobalMN = [];
MM2M = 1/1000;

%% load skeleton data if present and required
if  UserPref.useSkeleton && ~UserPref.isCombine && ...
    isfield(Data,'Skeletons') && ~isempty(Data.Skeletons)
    Skeleton = Data.Skeletons;
else
    fprintf('No Qualisys Skeleton imported.\n');
    return
end

%% Adapt conversion table
SkelData.srcType = 'qual_sports';
Def_body = loadBody('def_body');
ConversionTable = Def_body.conversion;
% copy the Qualisys fields to the generic fields (for convertCS)
ConversionTable.srcName = ConversionTable.([SkelData.srcType,'_name']);
ConversionTable.srcOrientation = ConversionTable.(SkelData.srcType);

%% Extract skeleton data
SegmentNames = Skeleton.SegmentLabels;
nSegmentNames = length(SegmentNames);
for srcJoint = 1:nSegmentNames
    Quat = squeeze(Skeleton.RotationData(:,srcJoint,:));   % orientations as quaternions
    Pos = squeeze(Skeleton.PositionData(:,srcJoint,:)) * MM2M; % positions as 3-vectors
    Quat = [Quat(4,:); Quat(1,:); Quat(2,:); Quat(3,:)]; % re-arrange components to comply to Matlab standard
    SkelData.signals.joints(srcJoint).name = SegmentNames{srcJoint};
    SkelData.signals.joints(srcJoint).dynamic.rot.units = '1';
    SkelData.signals.joints(srcJoint).dynamic.rot.data = zeros(size(Quat,2),9);
    for ss=1:size(Quat,2)
        Rotmat = kbquat2rotm(Quat(:,ss));
        SkelData.signals.joints(srcJoint).dynamic.rot.data(ss,:) = Rotmat(:);
    end
    SkelData.signals.joints(srcJoint).dynamic.pos.units = 'm';
    SkelData.signals.joints(srcJoint).dynamic.pos.data = Pos';
end

%% if SimInMN is not defined, we are reading for the body file: stop here
if isempty(SimInMN)
    % flags and messages
    IsSkeleton = true;
    fprintf('Imported Qualisys Skeleton.\n');
    return;
end

%% Make joint stream structure of myonardo-consistent joints and the internal sampling rate
Time = ((1:Data.Frames)-1)/Data.FrameRate;

% convert the spatial orientations of the segments into joint angles, according to the Myonardo
% definition (preBody: see prepareSimParam)
[Jointstream,GlobalStream] = convertCS(SkelData,Time,ConversionTable,UserPref.isGlobalPos,PreBody);

% write to stream and resample
if ~isempty(Jointstream)
    % SimInMN, GlobalMN are both MN data structures; GlobalMN is to create virtual markers and to 
    % apply measured forces
    SimInMN.signals.joints = Jointstream.signals.joints;
    GlobalMN.signals.joints = GlobalStream.signals.joints;

    % Resample
    NJoints = length(SimInMN.signals.joints);

    % The local rotations and translations
    for jj = 1:NJoints
        JointRot = Jointstream.signals.joints(jj).dynamic.rot.data;
        % resample to the internal sample rate (internalSR)
        TmpRot = interp1(QTMTime,JointRot,InternalSR); % interp1 accepts NSamples x dim data
        SimInMN.signals.joints(jj).dynamic.rot.data = TmpRot;

        % The local translations of specific joints
        if isfield(Jointstream.signals.joints(jj).dynamic,'trans')
            JointTrans = Jointstream.signals.joints(jj).dynamic.trans.data;
            TmpTrans = interp1(QTMTime,JointTrans,InternalSR); % interp1 accepts NSamples x dim data
            SimInMN.signals.joints(jj).dynamic.trans.data = TmpTrans;
        end

        % The local positions of specific joints
        if isfield(Jointstream.signals.joints(jj).dynamic,'pos')
            JointTrans = Jointstream.signals.joints(jj).dynamic.pos.data;
            TmpTrans = interp1(QTMTime,JointTrans,InternalSR); % interp1 accepts NSamples x dim data
            SimInMN.signals.joints(jj).dynamic.pos.data = TmpTrans;
        end
    end

    % The global positions and rotations
    for jj = 1:length(GlobalMN.signals.joints)

        % The global positions of GRF-related joints
        JointPos = GlobalStream.signals.joints(jj).dynamic.pos.data;
        JointDistal = GlobalStream.signals.joints(jj).dynamic.distal.data;
        JointGlobalR = GlobalStream.signals.joints(jj).dynamic.globalR.data;
        % Resample to the internal sample rate (internalSR)
        TmpPos = interp1(QTMTime,JointPos,InternalSR);
        TmpDistal = interp1(QTMTime,JointDistal,InternalSR);
        JointGlobalR = reshape(JointGlobalR,[length(QTMTime), 9]); % bring to NSamples x dim shape for interp1
        TmpGlobalR = interp1(QTMTime,JointGlobalR,InternalSR);
        TmpGlobalR = reshape(TmpGlobalR,[length(InternalSR), 3,3]); % reshape to rotation matrix
        GlobalMN.signals.joints(jj).dynamic.pos.data = TmpPos;
        GlobalMN.signals.joints(jj).dynamic.distal.data = TmpDistal;
        GlobalMN.signals.joints(jj).dynamic.globalR.data = TmpGlobalR;

        % save the segment length, in order to get the endpoint (preBody: see prepareSimParam)
        iSeg = PreBody.joints.followerSegment(jj);
        GlobalMN.signals.joints(jj).static.length.data = PreBody.segments.size(iSeg) * PreBody.general.height;
    end

    % flags and messages
    IsSkeleton = true;
    fprintf('Imported Qualisys Skeleton.\n');
else
    fprintf('No Qualisys Skeleton found\n');
end
end

