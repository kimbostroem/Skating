function UserBody = extractBodyFromQUALSkel(SkelData)
%% Estimate body segment lengths ("size") from Qualisys MAT file with Skeleton
% The segment size is scaled to the bodyHeight and to the default nextPos (see defBody),
% except for the pelvis. 
%
% SYNTAX
% estBody = extractBodyFromQUALSkel(SkelData,defBody)
%
% INPUT
%     SkelData        - (struct) containing the skeleton segments
%
% OUTPUT
%    UserBody	(struct) MN data structure
%
% Local functions: estimateBodyHeight
%
% See also: importQUALSkel
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Kim J Bostroem, Marc de Lussanet
% Version 221230 (MdL) comment
% Version 230418 (MdL) removed the calculations for nexPos for hip and lumbar joints, for this
% destroyes the muscle positions; correct scaling for the RBM and Pelvis segments (according to 
% marker import (extractBodyFromSegments); correct scale factor
% Version 230609 (MdL) improved comments, better variable names for spine
% version 230613 (MdL) do not calculate ribcage length: we use the default from def_body (see also convertCS)

%% init
UserBody = struct;
SkelSegData = SkelData.signals.joints;
SkelSegNames = {SkelSegData.name};
Y=2;

%% Adapt conversion table
SourceType = 'qual_sports';
DefBody = loadBody('def_body');
SizeDefault     = DefBody.segments.size; % default segment sizes
NextPosDefault  = DefBody.segments.nextPos;
ConversionTable = DefBody.conversion;
SegmentNames    = DefBody.segments.name;
% copy the Qualisys fields to the generic fields (for estimateBodyHeight)
ConversionTable.srcName = ConversionTable.([SourceType,'_name']);
ConversionTable.srcOrientation = ConversionTable.(SourceType);
SkelJntNames = {SkelData.signals.joints.name};

%% Calc body height
BodyHeight = estimateBodyHeight(ConversionTable, SkelSegData, SkelSegNames, DefBody);

%% Calc segment sizes
% init size array
SegmentSizes = nan(size(DefBody.segments.size));
for JointNo = 1:length(DefBody.joints.id) % No : the row number of the id
    JointId = DefBody.joints.id(JointNo);
    SegmentId = DefBody.joints.followerSegment(JointNo);
    SegmentNo = DefBody.segments.id == SegmentId;
    ConversionNo = (ConversionTable.id == JointId);
    thisNames = strsplit(ConversionTable.srcName{ConversionNo},';');
    switch JointId

        case 21 % RBM_joint -> (same position and size as) pelvis
            % Pelvis (at mid-asi) and RBM (at mid-hip) have the same orientation and size
            % the size is the "default asi distance" but scaled to the actual hip distance
            PlvNo = startsWith(SegmentNames,'pelvis_');
            RBMNo = startsWith(SegmentNames,'RBM_');
            LFemurID = 27;
            RFemurID = 26;
            LHipNextNo = DefBody.segments.next(PlvNo,:) == LFemurID; % which entry of nextpos gives the left hip
            RHipNextNo = DefBody.segments.next(PlvNo,:) == RFemurID; % which entry of nextpos gives the right hip
            % pelvis default width scaled to the body height
            PelvisWidthDefault = SizeDefault(PlvNo) * BodyHeight; % [m]
            % The position of the hip joints is given by nextPos, expressed in the pelvis width
            NormalizedLHipPosDefault = squeeze(NextPosDefault(PlvNo, LHipNextNo, :));
            NormalizedRHipPosDefault = squeeze(NextPosDefault(PlvNo, RHipNextNo, :));
            NormalizedHipWidthDefault = abs(NormalizedLHipPosDefault(Y) - NormalizedRHipPosDefault(Y));
            HipDistanceDefault = PelvisWidthDefault * NormalizedHipWidthDefault;
            LHipNo = matches({SkelData.signals.joints.name},'LeftUpLeg');
            RHipNo = matches({SkelData.signals.joints.name},'RightUpLeg');
            LHipPos = SkelData.signals.joints(LHipNo).dynamic.pos.data;
            RHipPos = SkelData.signals.joints(RHipNo).dynamic.pos.data;
            HipDistance = mean(vecnorm(LHipPos - RHipPos, 2,2),'omitnan');
            
            % Pelvis scaling factor in X
            % PelvisScalingFactor_X = HipDistance / HipDistanceDefault;
            % Pelvis scaling factor in Y
            PelvisScalingFactor_Y = HipDistance / HipDistanceDefault;
            % Pelvis scaling factor in Z
            % PelvisScalingFactor_Z = HipDistance / HipDistanceDefault;
            
            PelvisScalingFactor = PelvisScalingFactor_Y;
            PelvisWidth = PelvisWidthDefault * PelvisScalingFactor;
            % set the sizes
            SegmentSizes(PlvNo) = PelvisWidth / BodyHeight;
            SegmentSizes(RBMNo) = PelvisWidth / BodyHeight;

        case 12 % lumbar joint -> thorax
            % The Qualisys sports skeleton has three back segments, Spine, Spine1 and Spine2.
            % The joint between Spine1 and Spine2 is stiff (see Sports-Marker-Set.pdf).
            % The translation of the lumbar joint is dealt with in convertCS.

        case {2,3} % right/left clavicle
            % The length of the clavicle is determined by the distance between the clavicle joint and the shoulder joint
            followerJoint = JointId + 4; % right/left shoulder joint
            followerName  = ConversionTable.srcName{ConversionTable.id == followerJoint};
            SkelThis1No   = strcmp(thisNames{1}, SkelJntNames);
            SkelFollowNo  = strcmp(followerName, SkelJntNames);
            IsFlankerSegmentsDefined = any(SkelThis1No) && any(SkelFollowNo);
            if IsFlankerSegmentsDefined
                % the length in [m]
                SegmentSizes(SegmentNo) = vecnorm(...
                    SkelData.signals.joints(SkelFollowNo).dynamic.pos.data(1,:)...
                    - SkelData.signals.joints(SkelThis1No).dynamic.pos.data(1,:)); % first sample
                % ignore the nextpos translations from defBody
                % 2nd follower is shoulder (1st is scapula), and from it the z-component only
                nextPosZ       = vecnorm(DefBody.segments.nextPos(SegmentNo,2,3)); 
                ThoraxSegment  = DefBody.joints.baseSegment(JointId);
                 % 2nd follower to thorax is clavicle (1st is neck), and from it the y-component only
                nextPosBaseY   = abs(DefBody.segments.nextPos(ThoraxSegment,2,2));
                nextPosIgnore  = nextPosZ - nextPosBaseY;
                % normalize the segment length to the body height
                SegmentSizes(SegmentNo) = SegmentSizes(SegmentNo) / (BodyHeight * nextPosIgnore);
            end

        otherwise % all other joints
            if ~any(ConversionNo) || any(cellfun(@isempty,thisNames))
                % { 4, 5} scapula
                % C1-8, skull, L1-S1, pelvis
                continue
            end
            followerJoint = DefBody.joints.id(DefBody.joints.baseSegment==SegmentId);
            if isempty(followerJoint)
                % {1} neck (head)
                % {10,11} wrist (hand)
                % {19,20) mt (toes)
                continue
            end
            % { 6, 7} shoulder (humerus)
            % { 8, 9} elbow (radius)
            % {13,14} hip (femur)
            % {15,16} knee (tibia)
            % {17,18} ankle (foot)
            followerNames = strsplit(ConversionTable.srcName{ConversionTable.id == followerJoint(1)},';');
            FollowerNo    = strcmp(followerNames{1}, SkelJntNames);
            SkelThis1No   = strcmp(thisNames{1}, SkelJntNames);
            IsFlankerSegmentsDefined = any(FollowerNo) && any(SkelThis1No);
            if IsFlankerSegmentsDefined
                % the length in [m]
                SegmentSizes(SegmentNo) = vecnorm(SkelData.signals.joints(FollowerNo).dynamic.pos.data(1,:) - SkelData.signals.joints(SkelThis1No).dynamic.pos.data(1,:));
                % scale size to the z-distance to next segment normalized to size
                if all(~isnan(DefBody.segments.nextPos(SegmentNo,1,:))) % only if nextPos is defined
                    % ignore the nextpos translations from defBody
                    nextPosZ = abs(DefBody.segments.nextPos(SegmentNo,1,3));
                    % normalize the segment length to the body height
                    SegmentSizes(SegmentNo) = SegmentSizes(SegmentNo) / (BodyHeight * nextPosZ);
                end
            end
    end
    % if the size is not defined, it's zero and should be set to NaN
    if SegmentSizes(SegmentNo) == 0
        SegmentSizes(SegmentNo) = NaN;
    end
end

%% Post-process sizes
% calc scaling factor to corresponding size in default body
ScaleToDefault = SegmentSizes ./ SizeDefault;

%% save into estimated body structure
UserBody.general.height = BodyHeight;
UserBody.segments.id = DefBody.segments.id;
UserBody.segments.name = DefBody.segments.name;
UserBody.segments.size = SegmentSizes;
UserBody.segments.scale = ScaleToDefault;

% % print the lengths:
% fprintf('\n==================================================================================\n');
% fprintf('%5s %25s %8s %7s\n', 'id','changed sizes','len','scale')
% for Id = 1:length(UserBody.segments.id)
%     Idx = UserBody.segments.id == Id;
%     if any(Idx) && ~isnan(UserBody.segments.size(Idx))
%         fprintf('%5d %25s %5.1f cm %5.0f %%\n', Id, UserBody.segments.name{Idx}, ...
%             100*UserBody.segments.size(Idx)*BodyHeight, 100*UserBody.segments.scale(Idx) );
%     end
% end
% fprintf('====================================================================================\n\n');
end

% ==================================================================================================

function BodyHeight = estimateBodyHeight(ConversionTable, SkelSegData, SkelSegNames, DefBody)
%% Estimate body height

upperJoints  = [21, 12, 1]; % RBM, lumbar, neck
lowerJointsR = [21, 13, 15, 17]; % RBM, right hip, right knee, right ankle
lowerJointsL = [21, 14, 16, 18]; % RBM, left hip, left knee, left ankle

% check if all the joints that are flanking the necessary segments are present
DoCalcHeight = true; % default: calculate the body height
allJoints = unique([upperJoints lowerJointsR lowerJointsL]);
JointNames = {};
for j = 1:length(allJoints)
    JointNames = [JointNames, strsplit(ConversionTable.srcName{ConversionTable.id == allJoints(j)}, ';')];
end

% if any of the required JointNames are missing, the height cannot be calculated
MissingSegmentIndex = ~matches(JointNames, SkelSegNames);
if any(MissingSegmentIndex)
    Fmt = ['Some segments needed for calculating body height are missing:' repmat(' "%s"', 1, sum(MissingSegmentIndex)) '\n'];
    warning( Fmt, JointNames{MissingSegmentIndex});
    DoCalcHeight = false;
end

% if all necessary segments are present calculate the height
if DoCalcHeight
    %% upper body
    JointIds = upperJoints;
    upperHeight = 0;
    % create a cell array of the the segment names, splitting the combined segments of the conversion table
    JointNames = {};
    for j = 1:length(JointIds)
        JointNames = [JointNames, strsplit(ConversionTable.srcName{ConversionTable.id == JointIds(j)}, ';')];
    end
    % RBM-Lumbar + Lumbar-Neck (i.e. ignoring the head segment which is not included in the description)
    for k = 1:length(JointNames)-1
        ProxJointIndex = strcmp(SkelSegNames, JointNames{k});
        DistJointIndex = strcmp(SkelSegNames, JointNames{k+1});
        Len = median(vecnorm(SkelSegData(ProxJointIndex).dynamic.pos.data - SkelSegData(DistJointIndex).dynamic.pos.data, 2, 2), 'omitnan');
        if k<length(JointNames)-1
            upperHeight = upperHeight + Len;
        else
            upperHeight = upperHeight + 1*Len;
        end
        % fprintf('dist between %s-%s = %.2f\n',SkelSegData(ProxJointIndex).name,SkelSegData(DistJointIndex).name,Len);
    end
    % Add default head height
    HeadHeight = DefBody.segments.size(1) * DefBody.general.height; %  default head size
    upperHeight = upperHeight + HeadHeight;

    %% Right legs
    JointIds = lowerJointsR;
    lowerHeightR = 0;
    JointNames = {};
    for j = 1:length(JointIds)
        myNames = strsplit(ConversionTable.srcName{ConversionTable.id == JointIds(j)}, ';');
        JointNames = [JointNames, myNames{1}];
    end
    for k = 1:length(JointNames)-1
        ProxJointIndex = strcmp(SkelSegNames, JointNames{k});
        DistJointIndex = strcmp(SkelSegNames, JointNames{k+1});
        Len = median(vecnorm(SkelSegData(ProxJointIndex).dynamic.pos.data - SkelSegData(DistJointIndex).dynamic.pos.data, 2, 2), 'omitnan');
        lowerHeightR = lowerHeightR + Len;
        % fprintf('dist between %s-%s = %.2f\n',SkelSegData(ProxJointIndex).name,SkelSegData(DistJointIndex).name,Len);
    end

    %% left legs
    JointIds = lowerJointsL;
    lowerHeightL = 0;
    JointNames = {};
    for j = 1:length(JointIds)
        JointNames = [JointNames, strsplit(ConversionTable.srcName{ConversionTable.id == JointIds(j)}, ';')]; %#ok<*AGROW>
    end
    for k = 1:length(JointNames)-1
        ProxJointIndex = strcmp(SkelSegNames, JointNames{k});
        DistJointIndex = strcmp(SkelSegNames, JointNames{k+1});
        Len = median(vecnorm(SkelSegData(DistJointIndex).dynamic.pos.data - SkelSegData(ProxJointIndex).dynamic.pos.data, 2, 2), 'omitnan');
        lowerHeightL = lowerHeightL + Len;
        % fprintf('dist between %s-%s = %.2f\n',SkelSegData(ProxJointIndex).name,SkelSegData(DistJointIndex).name,Len);
    end

    % calc body height
    BodyHeight = upperHeight + max([lowerHeightL, lowerHeightR]);
    BodyHeight = round(BodyHeight, 2); % round to cm
    fprintf('Estimated height: %.2f m (estimated from Qualisys Skeleton)\n', BodyHeight);
else
    BodyHeight = DefBody.general.height;
    warning('height could not be calculated from Skeleton. Using default value (%.2fm)', BodyHeight);
end
end



