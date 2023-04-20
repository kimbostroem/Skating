function  JointStream = calcJointStream(DefBody, Segments, Info, IsFootSet)
%% convert the coordinate systems into the stream.joints
% 
% SYNTAX
%   JointStream = calcJointStream(DefBody, Segments, Info, IsFootSet)
%
% INPUT
%     DefBody       - (table)
%     Segments      - (struct) Coordinate systems of the segments
%     Info          - (struct) information
%     IsFootSet     - (logical) has the extended foot markerset been detected
%
% OUTPUT
%    JointStream    - (struct) Containing estimated segment length, body height and weight (user_body)
%
% Local functions: extractBodyFromSegments, estimateBodyWeight, calcSegmentSize, scaleToMeter
%
% See also: importQUAL
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230214 (MdL) further correction to rotations of clavicle
% version 230227 (MdL) rotate the foot to compensate nextpos in the toes

%% Init
NSintQTM = Info.NSintQTM;

%% Base joint RBM
PelvisJointOri = createAngles(repmat(eye(3),[1 1 NSintQTM]),Segments.pelvis,'xyz','RBM'); % eye() returns the identity matrix
MidHipPos      = squeeze(Segments.lFemur(4,:,:) + Segments.rFemur(4,:,:)) / 2;

% lumbar_joint : the pelvis-to thorax joint is currently modeled as a 3d joint + translation
% The translation reflects the series of joints of the lumbar vertebrae
% 1. calculate the joint rotation
thorax_pelvis  = createAngles(Segments.pelvis,Segments.thorax,'xyz','Lumbar');
% 2. The translatory component, i.e. from the end of the pelvis segment (pos 5) and
%    the base of the thorax segment (pos 4).
if isfield(Segments,'lumbar') && size(Segments.lumbar,1) >= 5
    LumbarJointTranslation = squeeze(Segments.lumbar(5,:,:) - Segments.lumbar(4,:,:));
else
    LumbarJointTranslation = [];
end

%% neck_joint
skull_thorax = createAngles(Segments.thorax, Segments.skull,'xyz','Neck');
% 2. The translatory component, i.e. from the end of the pelvis segment (pos 5) and
%    the base of the thorax segment (pos 4).
if isfield(Segments,'neck') && size(Segments.neck,1) >= 5
    NeckJointTranslation = squeeze(Segments.neck(5,:,:) - Segments.neck(4,:,:));
else
    NeckJointTranslation = [];
end

%% clavicle_joint
% The default Myonardo has the clavicle to the side and 20 deg backwards. Get the orientation with
% respect to to the default orientation
% - Rotate the z-axis about 90 deg about the thorax x-axis
ThoraxSegmentLeft  = rotateCS(Segments.thorax,   -90, Segments.thorax,  'x'); % rotate 90 deg left
ThoraxSegmentRight = rotateCS(Segments.thorax,    90, Segments.thorax,  'x'); % rotate 90 deg right
% - Compensate for the Myonardo.slx hard-coded rotation about the vertical axis
ThoraxSegmentLeft  = rotateCS(ThoraxSegmentLeft,  20, Segments.thorax, 'z');
ThoraxSegmentRight = rotateCS(ThoraxSegmentRight,-20, Segments.thorax, 'z');
% - Turn around the left clavicle about its z-axis such that the y points up again
ThoraxSegmentLeft  = rotateCS(ThoraxSegmentLeft, 180, ThoraxSegmentLeft, 'z'); % invert about long axis
% Calculate the clavicle joint angle
left_clavicle_joint  = createAngles(ThoraxSegmentLeft, Segments.lClavicle,'xyz','lClav');
right_clavicle_joint = createAngles(ThoraxSegmentRight,Segments.rClavicle,'xyz','rClav');

%% shoulder_joint
% The default upper arm is to the side, i.e. rotated 20 deg forwards with respect to the clavicle
lClavicleCSrot = Segments.lClavicle;
% Bring the right clavicle in the correct orientation by a turn about its long axis
rClavicleCSrot = Segments.rClavicle;
rClavicleCSrot = rotateCS(rClavicleCSrot, 180,   rClavicleCSrot, 'z');
% - Compensate for the Myonardo.slx hard-coded rotation about the vertical axis
lClavicleCSrot = rotateCS(lClavicleCSrot, -20,   Segments.thorax,'z');
rClavicleCSrot = rotateCS(rClavicleCSrot,  20,   Segments.thorax,'z');

% the shoulder joint with respect to the clavicle
right_shoulder = createAngles(rClavicleCSrot,    Segments.rHumerus,'xyz','rShoulder');
left_shoulder  = createAngles(lClavicleCSrot,    Segments.lHumerus,'xyz','lShoulder');
right_elbow    = createAngles(Segments.rHumerus, Segments.rRadius, 'xyz','rElbow'); %yxz
left_elbow     = createAngles(Segments.lHumerus, Segments.lRadius, 'xyz','lElbow'); %yxz
right_wrist    = createAngles(Segments.rRadius,  Segments.rHand,   'xyz','rWrist');
left_wrist     = createAngles(Segments.lRadius,  Segments.lHand,   'xyz','lWrist');

%% hip_joint
rFemurCSrot = rotateCS(Segments.rFemur,180,Segments.rFemur,   'y');
lFemurCSrot = rotateCS(Segments.lFemur,180,Segments.lFemur,   'y');
right_hip   = createAngles(Segments.pelvis, rFemurCSrot,'xyz','rHip');
left_hip    = createAngles(Segments.pelvis, lFemurCSrot,'xyz','lHip');
right_hip   = right_hip .* repmat([-1 1 -1]',[1, NSintQTM]);   % MdL/KB: we do not know what this old fix is good for
left_hip    = left_hip  .* repmat([-1 1 -1]',[1, NSintQTM]);   % MdL/KB: we do not know what this old fix is good for

%% knee, ankle
right_knee  = createAngles(Segments.rFemur,     Segments.rTibia,  'xyz','rKnee');
left_knee   = createAngles(Segments.lFemur,     Segments.lTibia,  'xyz','lKnee');
% Rotate the joint such that with an angle of 0, the foot segment is horizontally forward
AnkleAngle  = -90; % deg
% Correct the ankle angle for the translation by nextPos in the toes (this fix should be undone for the extended foot)
RFootNo     = contains([DefBody.segments.name],'right_foot_segment');
X           = 1;
AnkleAngle  = AnkleAngle + atand(DefBody.segments.nextPos(RFootNo,1,X));
rTibiaCSrot = rotateCS(    Segments.rTibia, AnkleAngle, Segments.rTibia,  'y');
lTibiaCSrot = rotateCS(    Segments.lTibia, AnkleAngle, Segments.lTibia,  'y');
right_ankle = createAngles(rTibiaCSrot, Segments.rFoot,   'xyz','rAnkle');
left_ankle  = createAngles(lTibiaCSrot, Segments.lFoot,   'xyz','lAnkle');

%% Copy joint angles into the joints stream
JointStream( 1) = makeMNJoint(DefBody,'neck_joint',          skull_thorax,  NeckJointTranslation);
JointStream( 2) = makeMNJoint(DefBody,'right_clavicle_joint',right_clavicle_joint);
JointStream( 3) = makeMNJoint(DefBody,'right_shoulder_joint',right_shoulder);
JointStream( 4) = makeMNJoint(DefBody,'right_elbow_joint',   right_elbow);   
JointStream( 5) = makeMNJoint(DefBody,'right_wrist_joint',   right_wrist);   
JointStream( 6) = makeMNJoint(DefBody,'left_clavicle_joint', left_clavicle_joint);
JointStream( 7) = makeMNJoint(DefBody,'left_shoulder_joint', left_shoulder); 
JointStream( 8) = makeMNJoint(DefBody,'left_elbow_joint',    left_elbow);         
JointStream( 9) = makeMNJoint(DefBody,'left_wrist_joint',    left_wrist);   
JointStream(10) = makeMNJoint(DefBody,'lumbar_joint',        thorax_pelvis, LumbarJointTranslation);
JointStream(11) = makeMNJoint(DefBody,'right_hip_joint',     right_hip);    
JointStream(12) = makeMNJoint(DefBody,'right_knee_joint',    right_knee);   
JointStream(13) = makeMNJoint(DefBody,'right_ankle_joint',   right_ankle);  
JointStream(14) = makeMNJoint(DefBody,'left_hip_joint',      left_hip);     
JointStream(15) = makeMNJoint(DefBody,'left_knee_joint',     left_knee);    
JointStream(16) = makeMNJoint(DefBody,'left_ankle_joint',    left_ankle);   
JointStream(17) = makeMNJoint(DefBody,'RBM_joint',           PelvisJointOri, MidHipPos); 

%% Rotate the Neck and Lumbar segments into local coordinates
% (Neck, Lumbar (and RBM) segments now have their orientation in world coordinates)
% Neck segment
if isfield(JointStream( 1).dynamic,'trans') && ~isempty(JointStream( 1).dynamic.trans.data)
    WorldTranslation = JointStream( 1).dynamic.trans.data;
    BaseR = squeeze(Segments.thorax(1:3,:,:)); % time series of the rotation matrix
    LocalTranslation = WorldTranslation;
    for i = 1:NSintQTM
        Local = squeeze(BaseR(:,:,i)) * LocalTranslation(i,:)';
        JointStream( 1).dynamic.trans.data(i,:) = Local;
    end
end
% Lumbar segment
if isfield(JointStream(10).dynamic,'trans') && ~isempty(JointStream(10).dynamic.trans.data)
    WorldTranslation = JointStream(10).dynamic.trans.data;
    BaseR = squeeze(Segments.pelvis(1:3,:,:));
    LocalTranslation = WorldTranslation;
    for i = 1:NSintQTM
        Local = squeeze(BaseR(:,:,i)) * LocalTranslation(i,:)';
        JointStream(10).dynamic.trans.data(i,:) = Local;
    end
end

%% Extended (5-segment) foot 
if IsFootSet

    rShankCSrot = rotateCS(Segments.rShank,180,Segments.rFemur,   'y'); % RDub
    lShankCSrot = rotateCS(Segments.lShank,180,Segments.lFemur,   'y'); % RDub

    %% Alternative ankle joint calculations
    right_ankleV3D = createAngles(rShankCSrot, Segments.rFootV3D,   'xyz','rAnkleV3D');    % 
    left_ankleV3D  = createAngles(lShankCSrot, Segments.lFootV3D,   'xyz','lAnkleV3D');    % 
    right_tarsalShank = createAngles(rShankCSrot,   Segments.rHindfoot,    'xyz','rTarsShank');% 
    left_tarsalShank  = createAngles(lShankCSrot,   Segments.lHindfoot,    'xyz','lTarsShank');% 
    
    %%  Ghent foot model
    % alternative calculation: overwrite right and left ankle with hindfoot joint
    right_tarsalTibia = createAngles(rTibiaCSrot,   Segments.rHindfoot,    'xyz','rTarsTibia');% RDub jointID = 17
    left_tarsalTibia  = createAngles(lTibiaCSrot,   Segments.lHindfoot,    'xyz','lTarsTibia');%      jointID = 18
    right_ankle       = right_tarsalTibia;
    left_ankle        = left_tarsalTibia;
    right_midtarsal   = createAngles(Segments.rHindfoot,   Segments.rMidfoot,    'xyz','rMidTars');%
    left_midtarsal    = createAngles(Segments.lHindfoot,   Segments.lMidfoot,    'xyz','lMidTars');%
    right_midtarsalmtmed = createAngles(Segments.rMidfoot,    Segments.rMedforefoot,'xyz','rTarsMed');% (tarsal - medial metatarsal)
    left_midtarsalmtmed  = createAngles(Segments.lMidfoot,    Segments.lMedforefoot,'xyz','lTarsMed');% (tarsal - medial metatarsal)
    right_midtarsalmtlat = createAngles(Segments.rMidfoot,    Segments.rLatforefoot,'xyz','rTarsML'); % (tarsal - lateral metatarsal)
    left_midtarsalmtlat  = createAngles(Segments.lMidfoot,    Segments.lLatforefoot,'xyz','lTarsML'); % (tarsal - lateral metatarsal)
    right_mtmedmtlat  = createAngles(Segments.rMedforefoot,Segments.rLatforefoot,'xyz','rMTarsML');% (medial metatarsal - lateral metatarsal)
    left_mtmedmtlat   = createAngles(Segments.lMedforefoot,Segments.lLatforefoot,'xyz','lMTarsML');% (medial metatarsal - lateral metatarsal)
    right_mtp1        = createAngles(Segments.rMedforefoot,Segments.rHallux,     'xyz','rMTP1');   % (medial metatarsal- hallux)
    left_mtp1         = createAngles(Segments.lMedforefoot,Segments.lHallux,     'xyz','lMTP1');   % (medial metatarsal- hallux)
    
    %% copy the angles into the stream
    % overwrite the original ankle definition
    JointStream(13) = makeMNJoint(DefBody,'right_ankle_joint',    right_ankle);% (tibia - hindfoot)
    JointStream(16) = makeMNJoint(DefBody,'left_ankle_joint',    left_ankle);  % (tibia - hindfoot)
    % seek the counter
    st = length(JointStream);
    JointStream(st+ 1) = makeMNJoint(DefBody,'right_midtarsal_joint',        right_midtarsal);  % (midfoot-hindfoot)
    JointStream(st+ 2) = makeMNJoint(DefBody,'right_midtarsal_medmeta_joint',right_midtarsalmtmed);% (midfoot-medial forefoot)
    JointStream(st+ 3) = makeMNJoint(DefBody,'right_midtarsal_latmeta_joint',right_midtarsalmtlat);% (midfoot-lateral forefoot)
    JointStream(st+ 4) = makeMNJoint(DefBody,'right_medmeta_latmeta_joint',  right_mtmedmtlat); % (medial forefoot - lateral forefoot)
    JointStream(st+ 5) = makeMNJoint(DefBody,'right_medmeta_p1_joint',       right_mtp1);       % (medial forefoot - hallux)
    JointStream(st+ 6) = makeMNJoint(DefBody,'left_midtarsal_joint',         left_midtarsal);   % RDub ... not sure if these numbers now need to change
    JointStream(st+ 7) = makeMNJoint(DefBody,'left_midtarsal_medmeta_joint', left_midtarsalmtmed);
    JointStream(st+ 8) = makeMNJoint(DefBody,'left_midtarsal_latmeta_joint', left_midtarsalmtlat);
    JointStream(st+ 9) = makeMNJoint(DefBody,'left_medmeta_latmeta_joint',   left_mtmedmtlat);
    JointStream(st+10) = makeMNJoint(DefBody,'left_medmeta_p1_joint',        left_mtp1);
    % alternative ankle joints
    JointStream(st+11) = makeMNJoint(DefBody,'right_ankleV3D_joint',   right_ankleV3D);    % (Shank- footV3D )
    JointStream(st+12) = makeMNJoint(DefBody,'left_ankleV3D_joint',    left_ankleV3D);     % (Shank- footV3D )
    JointStream(st+13) = makeMNJoint(DefBody,'right_tarsalShank_joint',right_tarsalShank); % (shank - hindfoot)
    JointStream(st+14) = makeMNJoint(DefBody,'left_tarsalShank_joint', left_tarsalShank);  % (tibia - hindfoot)

end
end
