function UserBody = extractBodyFromQUALSkel(SkelData)
%% Estimate body segment lengths ("size") from Qualisys MAT file with Skeleton
% The segment size is scaled to the bodyHeight and to the default nextPos (see defBody),
% except for the pelvis. 
% Pelvis : origin of the pelvis and the relative positions of hips and l5 are QTM-specific, and so
% the nextPos is changed for the UserBody
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

%% init
UserBody = struct;
SkelSegData = SkelData.signals.joints;
SkelSegNames = {SkelSegData.name};

%% Adapt conversion table
srcType = 'qual_sports';
DefBody = loadBody('def_body');
ConversionTable = DefBody.conversion;
% copy the Qualisys fields to the generic fields (for estimateBodyHeight)
ConversionTable.srcName = ConversionTable.([srcType,'_name']);
ConversionTable.srcOrientation = ConversionTable.(srcType);
SkelJntNames = {SkelData.signals.joints.name};

%% Calc body height
BodyHeight = estimateBodyHeight(ConversionTable, SkelSegData, SkelSegNames, DefBody);

%% Calc segment sizes
% init size array
mySizes = nan(size(DefBody.segments.size));
myNextPos = nan(size(DefBody.segments.nextPos));
for JointNo = 1:length(DefBody.joints.id) % No : the row number of the id
    JointId = DefBody.joints.id(JointNo);
    SegmentId = DefBody.joints.followerSegment(JointNo);
    SegmentNo = DefBody.segments.id == SegmentId;
    ConversionNo = (ConversionTable.id == JointId);
    thisNames = strsplit(ConversionTable.srcName{ConversionNo},';');
    switch JointId

        case 21 % RBM_joint -> (same position and size as) pelvis
            % The next joints of the hip joint are the lumbar joint and left/right hip joints.
            SkelThis1No  = strcmp(thisNames{1}, SkelJntNames);
            PelvisId      = 20;
            PelvisNo      = DefBody.segments.id == PelvisId;
            nextJoints    = DefBody.joints.id(DefBody.joints.baseSegment==20); % lumbar joint, left/right hip joint, L5_S1_joint (latter not needed)
            followerJoint = nextJoints(1); % lumbar joint
            followerNames = strsplit(ConversionTable.srcName{ConversionTable.id == followerJoint},';'); % 3 lumbar source joints
            LumbarNo      = strcmp(followerNames{1}, SkelJntNames);
            rHipName      = ConversionTable.srcName{ConversionTable.id == nextJoints(2)};
            lHipName      = ConversionTable.srcName{ConversionTable.id == nextJoints(3)};
            rHipNo        = strcmp(rHipName, SkelJntNames);
            lHipNo        = strcmp(lHipName, SkelJntNames);
            % position and orientation of the pelvis
            PelvisBase    = SkelData.signals.joints(SkelThis1No).dynamic.pos.data(1,:); % first sample
            PelvisRot     = SkelData.signals.joints(SkelThis1No).dynamic.rot.data(1,:,:); % rotation matrix
            % initialize PelvisSize to default
            PelvisSize    = DefBody.segments.size(SegmentNo);
            % if the legs are defined, calculate the hip positions and the segment length
            if any(rHipNo) && any(lHipNo)
                % position of the hips in global coordinates (at first sample)
                LeftHip  = SkelData.signals.joints(lHipNo).dynamic.pos.data(1,:); % first sample
                RightHip = SkelData.signals.joints(rHipNo).dynamic.pos.data(1,:); % first sample
                % pelvis size in [m]
                HipWidth = vecnorm(LeftHip - RightHip); % distance between left and right hip
                % normalize the segment length to the body height
                PelvisSize = HipWidth / BodyHeight;
                % The X(front) and Z(rostral) with respect to the pelvis coordinate system
                % get the position in local coordinates
                pelvisCS       = zeros(4,3);
                pelvisCS(1:3,:)= reshape(PelvisRot,[3 3]);
                pelvisCS(4,  :)= PelvisBase;
                LHipLocal = getLocalPosOnSegment(pelvisCS,LeftHip,'invert');
                RHipLocal = getLocalPosOnSegment(pelvisCS,RightHip,'invert');
                % nextPos for left and right hip: recalculate these conform the skeleton definition
                myNextPos(PelvisNo,2,:) = RHipLocal([2 1 3]); % offset right hip
                myNextPos(PelvisNo,3,:) = LHipLocal([2 1 3]); % offset left hip
            end
            % calculate the lumbar joint
            if any(LumbarNo)
                LumbarPos = SkelData.signals.joints(LumbarNo).dynamic.pos.data(1,:); % first sample
                LumbarVerticalOffset = vecnorm(PelvisBase - LumbarPos);
                myNextPos(PelvisNo,1,:) = [0 0 LumbarVerticalOffset]; % offset L5
            end
            % the RBM and pelvis segments are welded together and have the same length. The RBM is a
            % "mock-segment" with minimal mass, no gyr and no nextPos and the same position as the
            % pelvis
            mySizes(SegmentNo) = PelvisSize;
            mySizes(PelvisId) = PelvisSize;
            % normalize the nextPos to the pelvis size
            myNextPos(PelvisId,:,:) = myNextPos(PelvisId,:,:) / (BodyHeight * PelvisSize);

        case 12 % lumbar joint -> thorax
            % The lumbar joint corresponds 3 subsequent joints in the Skeleton, named Spine, Spine1 and Spine2. The distance
            % between the 1st and 2nd spine joints corresponds to the translation of the lumbar joint. The distance between the 2nd
            % and 3rd spine joint plus the distance between the 3rd spine joint and the neck joint corresponds to the length of the
            % thorax segment.
            SkelThis2No  = strcmp(thisNames{2}, SkelJntNames);
            SkelThis3No  = strcmp(thisNames{3}, SkelJntNames);
            SkelThisEndNo = strcmp(thisNames{end}, SkelJntNames);
            cervicalNames = strsplit(ConversionTable.srcName{ConversionTable.id == 1},';'); % Neck and Head
            SkelCervId    = strcmp(cervicalNames{1}, SkelJntNames);
            IsFlankerSegmentsDefined = any(SkelThis2No) && any(SkelThis3No) && any(SkelThisEndNo) && any(SkelCervId);
            if IsFlankerSegmentsDefined
                % the length in [m]
                Length = vecnorm(...
                    SkelData.signals.joints(SkelThis2No).dynamic.pos.data(1,:)...
                    - SkelData.signals.joints(SkelThis3No).dynamic.pos.data(1,:)) + ...
                    vecnorm( SkelData.signals.joints(SkelThisEndNo).dynamic.pos.data(1,:)...
                    - SkelData.signals.joints(SkelCervId).dynamic.pos.data(1,:));
                mySizes(SegmentNo) = Length;
                % ignore the nextpos translations from defBody
                nextPosZ = abs(DefBody.segments.nextPos(SegmentNo,1,3));
                % normalize the segment length to the body height
                mySizes(SegmentNo) = mySizes(SegmentNo) / (BodyHeight * nextPosZ);
            end

        case {25,26,27,28,29,30,31} % cervical vertebra joints -> cervical vertebra
            % The distance between the neck and head joint corresponds to the sum of the sizes of the 7 cervical vertebral
            JointNames = strsplit(ConversionTable.srcName{ConversionTable.id == 1},';');
            Skel1No  = strcmp(JointNames{1}, SkelJntNames);
            Skel2No  = strcmp(JointNames{2}, SkelJntNames);
            IsFlankerSegmentsDefined = any(Skel1No) && any(Skel2No);
            if IsFlankerSegmentsDefined
                % the length in [m]
                mySizes(SegmentNo) = vecnorm(SkelData.signals.joints(Skel2No).dynamic.pos.data(1,:) - SkelData.signals.joints(Skel1No).dynamic.pos.data(1,:))/7; % first sample
                % ignore the nextpos translations from defBody
                nextPosZ = abs(DefBody.segments.nextPos(SegmentNo,1,3));
                % normalize the segment length to the body height
                mySizes(SegmentNo) = mySizes(SegmentNo) / (BodyHeight * nextPosZ);
            end

        case {33,34,35,36,37} % lumbar vertebra joints -> lumbar vertebra
            % The distance between the 1st and 2nd lumbar source joint corresponds to the sum of the sizes of the 5 lumbar vertebral
            JointNames = strsplit(ConversionTable.srcName{ConversionTable.id == 12},';');
            Skel1No  = strcmp(JointNames{1}, SkelJntNames);
            Skel2No  = strcmp(JointNames{2}, SkelJntNames);
            IsFlankerSegmentsDefined = any(Skel1No) && any(Skel2No);
            if IsFlankerSegmentsDefined
                % the length in [m]
                mySizes(SegmentNo) = vecnorm(SkelData.signals.joints(Skel2No).dynamic.pos.data(1,:) - SkelData.signals.joints(Skel1No).dynamic.pos.data(1,:))/5; % first sample
                % ignore the nextpos translations from defBody
                nextPosZ = abs(DefBody.segments.nextPos(SegmentNo,1,3));
                % normalize the segment length to the body height
                mySizes(SegmentNo) = mySizes(SegmentNo) / (BodyHeight * nextPosZ);
            end

        case {2,3} % right/left clavicle
            % The length of the clavicle is determined by the distance between the clavicle joint and the shoulder joint
            followerJoint = JointId + 4; % right/left shoulder joint
            followerName  = ConversionTable.srcName{ConversionTable.id == followerJoint};
            SkelThis1No   = strcmp(thisNames{1}, SkelJntNames);
            SkelFollowNo  = strcmp(followerName, SkelJntNames);
            IsFlankerSegmentsDefined = any(SkelThis1No) && any(SkelFollowNo);
            if IsFlankerSegmentsDefined
                % the length in [m]
                mySizes(SegmentNo) = vecnorm(...
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
                mySizes(SegmentNo) = mySizes(SegmentNo) / (BodyHeight * nextPosIgnore);
            end

        otherwise % all other joints
            if ~any(ConversionNo) || any(cellfun(@isempty,thisNames))
                % { 4, 5} scapula
                continue
            end
            followerJoint = DefBody.joints.id(DefBody.joints.baseSegment==SegmentId);
            if isempty(followerJoint)
                % {1} neck (head)
                % {19,20) mt (toes)
                continue
            end
            % { 6, 7} shoulder (humerus)
            % { 8, 9} elbow (radius)
            % {10,11} wrist (hand)
            % {13,14} hip (femur)
            % {15,16} knee (tibia)
            % {17,18} ankle (foot)
            followerNames = strsplit(ConversionTable.srcName{ConversionTable.id == followerJoint(1)},';');
            FollowerNo    = strcmp(followerNames{1}, SkelJntNames);
            SkelThis1No   = strcmp(thisNames{1}, SkelJntNames);
            IsFlankerSegmentsDefined = any(FollowerNo) && any(SkelThis1No);
            if IsFlankerSegmentsDefined
                % the length in [m]
                mySizes(SegmentNo) = vecnorm(SkelData.signals.joints(FollowerNo).dynamic.pos.data(1,:) - SkelData.signals.joints(SkelThis1No).dynamic.pos.data(1,:));
                % scale size to the z-distance to next segment normalized to size
                if all(~isnan(DefBody.segments.nextPos(SegmentNo,1,:))) % only if nextPos is defined
                    % ignore the nextpos translations from defBody
                    nextPosZ = abs(DefBody.segments.nextPos(SegmentNo,1,3));
                    % normalize the segment length to the body height
                    mySizes(SegmentNo) = mySizes(SegmentNo) / (BodyHeight * nextPosZ);
                end
            end
    end
    % if the size is not defined, it's zero and should be set to NaN
    if mySizes(SegmentNo) == 0
        mySizes(SegmentNo) = NaN;
    end
end

%% Post-process sizes
% calc scaling factor to corresponding size in default body
scales = nan(size(DefBody.segments.scale));
for i = 1:length(scales)
    scales(i) = mySizes(i) * BodyHeight / (DefBody.segments.size(i) * DefBody.general.height);
end

%% save into estimated body structure
UserBody.general.height = BodyHeight;
UserBody.segments.id = DefBody.segments.id;
UserBody.segments.name = DefBody.segments.name;
UserBody.segments.size = mySizes;
UserBody.segments.nextPos = myNextPos;
UserBody.segments.scale = scales;

% % print the lengths and nextPos:
% fprintf('%5s %25s %5s %5s  np: %6s %6s %6s xx %6s %6s %6s xx %6s %6s %6s\n', ...
%     'id','name','len','scale','np1-x','np1-y','np1-z','np2-x','np2-y','np2-z','np3-x','np3-y','np3-z')
% for Id = 1:length(UserBody.segments.id)
%     No = find(UserBody.segments.id == Id);
%     if ~isempty(No)
%         fprintf('%5d %25s %5.2f %5.2f  np: %6.2f %6.2f %6.2f xx %6.2f %6.2f %6.2f xx %6.2f %6.2f %6.2f\n', ...
%             Id, UserBody.segments.name{No}, ...
%             UserBody.segments.size(No)*BodyHeight, UserBody.segments.scale(No), ...
%             UserBody.segments.nextPos(No,1,1), UserBody.segments.nextPos(No,1,2), UserBody.segments.nextPos(No,1,3), ...
%             UserBody.segments.nextPos(No,2,1), UserBody.segments.nextPos(No,2,2), UserBody.segments.nextPos(No,2,3), ...
%             UserBody.segments.nextPos(No,3,1), UserBody.segments.nextPos(No,3,2), UserBody.segments.nextPos(No,3,3) ...
%             );
%     end
% end

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


