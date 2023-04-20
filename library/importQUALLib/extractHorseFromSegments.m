function [SubjectBody] = extractHorseFromSegments(Segments, JointOffsets, DefBody, Markers, Present, MaxNSamplesToUse)
%% Estimate body segment lengths from marker-defined segments to fill the segments table
% Estimates for segments table:
% - the lengths
% - the nextpos
% .. and BodyHeight
%
% SYNTAX
% SubjectBody = extractBodyFromSegments(Segments, DefBody, Markers, Present, MaxNSamplesToUse, IsHorse)
%
% INPUT
%     Segments         (struct) Coordinate systems of the segments
%     DefBody          (table) sheet conversion from def_body.xlsx
%     Markers          (struct of 3 x N) time series for each marker that is present
%     Present          (struct of logical) for each marker: is it present
%     MaxNSamplesToUse (double) number of (non-NaN) samples to select for the stable pose at the 
%                      beginning of the time series
%
% OUTPUT
%     SubjectBody     (struct) MN data structure
%
% Local functions: getNextNames, estimateBodyHeightFromMarkers, scaleNextPos2UserVertColumn
%
% See also: extractBodyFromQUALSkel, extractBodyFromTheia, estimateSegLengths
% 
% (c) 2022 by Predimo GmbH

%% fill the segments sheet for the user_body table
% 
% version 221005 (MdL) provide MaxNSamplesToUse as argument for calcSegmentSize
% version 221011 (MdL) lengths for all segments
% version 221015 (MdL) fix zero-sized segments, fix name of neck and lumber
% version 221109 (MdL) outsourced from estimateSegLengths; completely reviewed
% version 221117 (MdL) checks; added JointOffsets

%% Params
P            = true; % do plot on error in calcSegmentSize
IsHorse      = false;

%% Init
NextPos       = nan(size(DefBody.segments.nextPos));
MxNS          = MaxNSamplesToUse;
DefaultLengths= [DefBody.segments.size];
MeanLengths   = nan(size(DefaultLengths));
SegmentNames  = DefBody.segments.name;

%% Estimate body height
DefaultHeight = DefBody.general.height;
SubjectHeight = estimateBodyHeightFromMarkers(Markers, Segments, DefaultHeight, MaxNSamplesToUse, IsHorse);

%% Get the segment lengths from Segment CSs

% Pelvis
PlvNo = startsWith(SegmentNames,'pelvis_');
PelvisWidth = JointOffsets.PelvisWidth;
MeanLengths(PlvNo) = PelvisWidth;

% RBM (the position of the body in the world) 
RBMNo = startsWith(SegmentNames,'RBM_');
MeanLengths(RBMNo) = MeanLengths(PlvNo)/1000; % [make its length negligible]

% Lumbar Back
% Most correct would be to scale the segment lengths and the radius of gyration. However this turns
% out tricky. Since the radius of gyration is not a very accurate representation any way, we do not
% scale the lumber segments, but rather the nextpos. For this we need the total lumbar length
LumbarNos = startsWith(SegmentNames,'L');
L1No      = startsWith(SegmentNames,'L1_');
LumbarLength = calcSegmentSize(L1No,Segments.lumbar, [], DefBody,MxNS,P);

% Ribcage
RibNo = startsWith(SegmentNames,'ribcage_');
ThoraxLength = calcSegmentSize(RibNo,Segments.thorax, [], DefBody,MxNS,P);
MeanLengths(RibNo) = ThoraxLength;

% Neck
% calculate the length of the individual segments, using the relative lengths from def_body
NeckNos = startsWith(SegmentNames, 'C'); % row no in segments table
C1No    = startsWith(SegmentNames,'C1_');
NeckLength = calcSegmentSize(C1No,Segments.neck, [], DefBody,MxNS,P);

% Head
Hd = startsWith(SegmentNames,'skull_');
MeanLengths(Hd) = calcSegmentSize(Hd,Segments.skull,  [], DefBody,MxNS,P);

% Arms: clavicle, humerus, radius
% Scapula
RS = startsWith(SegmentNames,'right_scapula_');
LS = startsWith(SegmentNames,'left_scapula_');
if Present.right_shoulder_back_marker
    MeanLengths(RS) = calcSegmentSize(RS,Segments.rHumerus, Markers.right_shoulder_back_marker, DefBody,MxNS,P);
end
if Present.left_shoulder_back_marker
    MeanLengths(LS) = calcSegmentSize(LS,Segments.lHumerus, Markers.left_shoulder_back_marker, DefBody,MxNS,P);
end

% Clavicle
RClavNo = startsWith(SegmentNames,'right_clavicle_');
LClavNo = startsWith(SegmentNames,'left_clavicle_');
MeanLengths(RClavNo) = calcSegmentSize(RClavNo,Segments.rClavicle,Segments.rHumerus,DefBody,MxNS,P);
MeanLengths(LClavNo) = calcSegmentSize(LClavNo,Segments.lClavicle,Segments.lHumerus,DefBody,MxNS,P);

% Humerus
RH = startsWith(SegmentNames,'right_humerus_');
LH = startsWith(SegmentNames,'left_humerus_');
MeanLengths(RH) = calcSegmentSize(RH,Segments.rHumerus, Segments.rRadius, DefBody,MxNS,P);
MeanLengths(LH) = calcSegmentSize(LH,Segments.lHumerus, Segments.lRadius, DefBody,MxNS,P);

% Radius
RR = startsWith(SegmentNames,'right_radius_');
LR = startsWith(SegmentNames,'left_radius_');
MeanLengths(RR) = calcSegmentSize(RR,Segments.rRadius,  Segments.rHand,   DefBody,MxNS,P);
MeanLengths(LR) = calcSegmentSize(LR,Segments.lRadius,  Segments.lHand,   DefBody,MxNS,P);

% Hand
RM = startsWith(SegmentNames,'right_hand_');
LM = startsWith(SegmentNames,'left_hand_');
MeanLengths(RM) = calcSegmentSize(RM,Segments.rHand, [], DefBody,MxNS,P);
MeanLengths(LM) = calcSegmentSize(LM,Segments.lHand, [], DefBody,MxNS,P);

% Legs
% femur
RF = startsWith(SegmentNames,'right_femur_');
LF = startsWith(SegmentNames,'left_femur_');
MeanLengths(RF) = calcSegmentSize(RF,Segments.rFemur,  Segments.rTibia, DefBody,MxNS,P);
MeanLengths(LF) = calcSegmentSize(LF,Segments.lFemur,  Segments.lTibia, DefBody,MxNS,P);
% tibia
RT = startsWith(SegmentNames,'right_tibia_');
LT = startsWith(SegmentNames,'left_tibia_');
MeanLengths(RT) = calcSegmentSize(RT,Segments.rTibia,  Segments.rFoot,  DefBody,MxNS,P);
MeanLengths(LT) = calcSegmentSize(LT,Segments.lTibia,  Segments.lFoot,  DefBody,MxNS,P);
% foot
RP = startsWith(SegmentNames,'right_foot_');
LP = startsWith(SegmentNames,'left_foot_');
MeanLengths(RP) = calcSegmentSize(RP,Segments.rFoot, Segments.rToes, DefBody,MxNS,P);
MeanLengths(LP) = calcSegmentSize(LP,Segments.lFoot, Segments.lToes, DefBody,MxNS,P);
% toes
RZ = startsWith(SegmentNames,'right_toes_');
LZ = startsWith(SegmentNames,'left_toes_');
MeanLengths(RZ) = calcSegmentSize(RZ,Segments.rToes, [], DefBody,MxNS,P);
MeanLengths(LZ) = calcSegmentSize(LZ,Segments.lToes, [], DefBody,MxNS,P);


%% Define the offsets (NextPos) of the next segments of the chain

% Pelvis joints (Pelvis_L5; l&r pelvis-thigh)
PelvNextNames = getNextNames(DefBody, PlvNo);
% Pelvis Origin is at mid ASI; x-ventral, y-left, z-caudal
% lumbar base: scale nextpos to segment length (=PelvisWidth)
L5No = startsWith(PelvNextNames,'L5_');
NextPos(PlvNo,L5No,:) = JointOffsets.L5 / PelvisWidth;
% Humerus base - hip joints: scale nextpos to segment length (=PelvisWidth)
LHipNo = startsWith(PelvNextNames,'left_femur_');
NextPos(PlvNo,LHipNo,:) = JointOffsets.LHip / PelvisWidth;
RHipNo = startsWith(PelvNextNames,'right_femur_');
NextPos(PlvNo,RHipNo,:) = JointOffsets.RHip / PelvisWidth;

% Ribcage / thorax joints (Thorax_C7; l&r pelvis-clavicle)
RibNextNames = getNextNames(DefBody, RibNo);
% TRhorax Origin is at Thorax Base (thoracal vertebra); x-ventral, y-left, z-caudal
% Clavicle joints: scale nextpos to segment length (=Thorax length)
NextNo = startsWith(RibNextNames,'left_clavicle_');
NextPos(RibNo,NextNo,:) = JointOffsets.lClavicle / ThoraxLength;
NextNo = startsWith(RibNextNames,'right_clavicle_');
NextPos(RibNo,NextNo,:) = JointOffsets.rClavicle / ThoraxLength;
NextNo = startsWith(RibNextNames,'C7_');
NextPos(RibNo,NextNo,:) = [0 0 1];

% Humerus
% Shoulder joint is caudal to clavicle endpoint
LClavNextNames = getNextNames(DefBody, LClavNo);
RClavNextNames = getNextNames(DefBody, RClavNo);
LClavNext = startsWith(LClavNextNames,'left_clavicle_');
RClavNext = startsWith(RClavNextNames,'right_clavicle_');
NextPos(LClavNo,LClavNext,:) = JointOffsets.lHumerus / MeanLengths(LClavNo);
NextPos(RClavNo,RClavNext,:) = JointOffsets.rHumerus / MeanLengths(RClavNo);

% Corrections for simple joints
NextPos(RBMNo,1,:) = [0 0 0]; % no translation
Nos = find(startsWith(SegmentNames,{'left_tibia_','right_tibia_','left_toes_','right_toes_'}));
for i=1:length(Nos)
    NextPos(Nos(i),1,:) = [0 0 1]; % translation to the end of the segment
end

% scapula
% radius
% hand
% patella
% thigh

% Scale the nextpos to compensate for the differences in the length of the vertebral column 
% (rather than the gyr and segment size)
NextPos(LumbarNos,1,3) = scaleNextPos2UserVertColumn(DefBody, LumbarNos, LumbarLength, SubjectHeight);
NextPos(NeckNos,  1,3) = scaleNextPos2UserVertColumn(DefBody, NeckNos,   NeckLength,   SubjectHeight);

% fprintf('==== extractBodyFromSegments:\n');
% nextpostest = reshape(permute(NextPos,[1 3 2]),[36,9])
% nextpos1ztest = squeeze(NextPos(:,1,3))

%% normalize sizes
MeanLengths = MeanLengths / SubjectHeight;
ScaledLen = (MeanLengths * SubjectHeight) ./ (DefBody.segments.size * DefBody.general.height);

%% Write estimated segment sizes
SubjectBody.segments.id      = DefBody.segments.id;
SubjectBody.segments.name    = DefBody.segments.name;
SubjectBody.segments.scale   = ScaledLen;
SubjectBody.segments.size    = MeanLengths;
SubjectBody.segments.nextPos = NextPos;
SubjectBody.general          = DefBody.general; % initialize substructure with default values
end

%% ========================================================================

function NextNames = getNextNames(DefBody, PlvNo)
%% get the names of the segments that follow on the current row number in the segments table
% version 221024 (MdL) first version

NextNos = zeros(3,1);
for i = 1:3
    NextNos(i) = find(DefBody.segments.id == DefBody.segments.next(PlvNo,i));
end
NextNames = DefBody.segments.name(NextNos);
end

%% ========================================================================

function Height = estimateBodyHeightFromMarkers(Markers,Segments,DefaultHeight,MaxNSamplesToUse,IsHorse)
%% Estimate the body height from the head markers
% 
% version 221015 (MdL) use the 5th dimension of the skull segment rather than the markers; implement horse
% version 221026 (MdL) implement selectStablePoseWindow
% version 221107 (MdL) select vertical component only

% stop if no markers are present
Height = DefaultHeight;
if IsHorse
    MinPlausibleHeight = 0.5;
    MaxPlausibleHeight = 2.5;
else
    MinPlausibleHeight = 1;
    MaxPlausibleHeight = 2.2;
end

% Get body height, if the necessary information is present
if IsHorse
    if ~all(isnan(Markers.thorax_05_marker),'all')
        HeadTop = Markers.thorax_05_marker;
    else
        HeadTop = Markers.thorax_16_marker;
    end
elseif all(isnan(Segments.skull(5,:,:)),'all')
    warning('Not enough head markers found: using default value of %.2f',Height);
    return;
else
    HeadTop = squeeze(Segments.skull(5,:,:));
end
% select the initial non-nan range of length MaxNSamplesToUse
HeightArray = selectStablePoseWindow(HeadTop,MaxNSamplesToUse);
HeightArray = HeightArray(3,:);
HeightMean = mean(HeightArray,'omitnan');
HeightStd = std(HeightArray,'omitnan');

IsPlausible = ...
    HeightMean < MaxPlausibleHeight && ...
    HeightMean > MinPlausibleHeight && ...
    HeightStd  < 0.05;
if IsPlausible 
    Height  = round(HeightMean,2); % round to cm
    fprintf('Estimated height %.0f +/- %.1fcm\n',100*Height,100*HeightStd);
else
    warning('The estimated height %.2f +/- %.2fm is not plausible -> using default value of %.2f\n',HeightMean,HeightStd,Height);
end
end

% ==================================================================================================

function NextPos_Z_User = scaleNextPos2UserVertColumn(DefBody, VertebralNos, SegmentLength, BodyHeight)
%% correct the radius of gyration of neck and lumbar vertebrae
% Note to radius of gyration (see example of 221105, MdL): 
% - Moment of inertia with respect to centre of mass: I = m * r_(x,y,z)^2
% - Radius of gyration g = √(I/m) with unit [m]
% - Since the radius r is normalized to the segment length (i.e. r_z), it follows :
%     Suppose ellipsoid with r_(x0,y0,z0) is scaled by factor F 
%     - g_x, g_y scale with 1/√F
%     - g_z scales with 1/F
%% Since this does not work out correctly, we scale the nextpos rather than the gyr and segment size

% the total length of the segment according to def_body
SegmentLengthNormDefault = sum(DefBody.segments.size(VertebralNos) .* DefBody.segments.nextPos(VertebralNos,1,3));
% the total length of the intervertebral spaces according to nextpos
NormIntervertebralSpace = DefBody.segments.nextPos(VertebralNos,1,3) - 1; % this is a factor so minus 1 is the difference
SegmentDisksNormDefault = sum(DefBody.segments.size(VertebralNos) .* NormIntervertebralSpace);
% the total length of the segment user body
SegmentLengthNormUser = (SegmentLength/BodyHeight);
% the difference in size that is to be compensated by scaling the nextpos
NextPosCorrection = SegmentLengthNormUser - SegmentLengthNormDefault;
% the scaling factor
ScaleFactorNextPos = (NextPosCorrection + SegmentDisksNormDefault) / SegmentDisksNormDefault;
% scale the intervertebral space and add one to obtain a scaling factor again
NextPos_Z_User = (NormIntervertebralSpace * ScaleFactorNextPos) + 1;
% Totallength_check = sum(DefBody.segments.size(VertebralNos) .* NextPos_Z_User)
end

