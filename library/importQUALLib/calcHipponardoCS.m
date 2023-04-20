function Segments = calcHipponardoCS(Markers, Present, Info)
%% calculate the time series of the segment coordinate systems from the markar data for hipponardo
%
% SYNTAX
%   [Segments, mid_asi] = calcHipponardoCS(Markers,Present,RefBody, Info)
%
% INPUT
%     Markers  - (struct of 3 x N) time series for each marker that is present
%     Present  - (struct of logical) for each marker: is it present
%     Info     - (stuct) information, e.g. NSintQTM = number of samples; Tail for interpolation
%
% OUTPUT
%    Segments (struct) Coordinate systems of the segments
%
% Local functions: 
%
% See also: importQUAL, createCS
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Version 210830 based on calcMainCS
% Version 221015 define end points of distal segments, define neck and lumbar spine

%% Init
NSamples          = Info.NSintQTM;
Segments.skull    = nan(5,3,NSamples); % 3x3 orientation & 3x1 base vector & 3x1 end vector
Segments.neck     = nan(4,3,NSamples);
Segments.thorax   = nan(4,3,NSamples);
Segments.lumbar   = nan(4,3,NSamples);
Segments.pelvis   = nan(4,3,NSamples);
Segments.rClavicle= nan(4,3,NSamples);
Segments.lClavicle= nan(4,3,NSamples);
Segments.rHumerus = nan(4,3,NSamples);
Segments.lHumerus = nan(4,3,NSamples);
Segments.rRadius  = nan(4,3,NSamples);
Segments.lRadius  = nan(4,3,NSamples);
Segments.rHand    = nan(4,3,NSamples);
Segments.lHand    = nan(4,3,NSamples);
Segments.rFemur   = nan(4,3,NSamples);
Segments.lFemur   = nan(4,3,NSamples);
Segments.rTibia   = nan(4,3,NSamples);
Segments.lTibia   = nan(4,3,NSamples);
Segments.rFoot    = nan(4,3,NSamples);
Segments.lFoot    = nan(4,3,NSamples);
Segments.rToes    = nan(4,3,NSamples);
Segments.lToes    = nan(4,3,NSamples);

%% Do not continue if no markers are present
% convert the Present structure to logical array
IsPresent = cell2mat(struct2cell(Present));
if ~any(IsPresent)
    return
end

%% Pelvis and lumbar spine

% Die Wirbelsäule hat keine lordose, und verläuft in etwa zwischen den
% hueftmarkern, und zwischen den scapula markern
PelvisFrontMarkersPresent  = Present.left_asi_marker && Present.right_asi_marker;
HipMarkersPresent          = Present.left_femur_lat_epi_marker && Present.right_femur_lat_epi_marker;
KneeMarkersPresent         = Present.left_lat_malleolus_marker && Present.right_lat_malleolus_marker;
ShoulderBackMarkersPresent = Present.left_shoulder_back_marker && Present.right_shoulder_back_marker;

%% Pelvis right to left vector and midpoint
if PelvisFrontMarkersPresent
    PelvisRToLVector = Markers.left_asi_marker - Markers.right_asi_marker;
else % head of femur markers (not so good)
    PelvisRToLVector = Markers.left_femur_trochanter_marker - Markers.right_femur_trochanter_marker;
end

%% midpoints of hip and leg markers
MidAsi =      (Markers.left_asi_marker              + Markers.right_asi_marker)/2; 
MidHip =      (Markers.left_femur_trochanter_marker + Markers.right_femur_trochanter_marker)/2;
MidKnees =    (Markers.left_femur_lat_epi_marker    + Markers.right_femur_lat_epi_marker) /2;
MidAnkles =   (Markers.left_lat_malleolus_marker    + Markers.right_lat_malleolus_marker) /2;
MidShoulder = (Markers.left_shoulder_back_marker    + Markers.right_shoulder_back_marker) /2;

%% Vektoren
% Rostral von vordere Huefte zum Schulter
SpineVector = MidShoulder - MidAsi;
% Schulter Vektor richtung ventral in sagittal Ebene
if     Present.thorax_05_marker
    ThoraxVentralVector = MidShoulder - Markers.thorax_05_marker;
    SchulterRefMarker   = Markers.thorax_05_marker;
else
    ThoraxVentralVector = MidShoulder - Markers.thorax_16_marker;
    SchulterRefMarker   = Markers.thorax_16_marker;
end
% Pelvis Vektor richtung ventral in sagittal Ebene
if Present.sacrum_marker
    HipVentralVector = MidKnees - Markers.sacrum_marker;
else
    HipVentralVector = MidKnees - Markers.sacrum_3_marker;
end
PelvisYVector = normLength(Markers.left_asi_marker - Markers.right_asi_marker);
PelvisLength  = meanDistance(MidHip, MidAsi);
% x-ventral, y-left, z-rostral
if PelvisFrontMarkersPresent
    Segments.pelvis = createCS(MidAsi, PelvisRToLVector,    SpineVector, 'yxz');
else
    PelvisRostralVector  = Markers.sacrum_marker - Markers.sacrum_3_marker;
    Segments.pelvis = createCS(MidAsi, PelvisRostralVector, HipVentralVector, 'zyx');
end
% Redefine the rostral (Z) Vector from the calculated pelvis CS
PelvisRostralVector = squeeze(Segments.pelvis(3,:,:));
% Redefine PelvisBase as the caudal end of the sacrum, with length of Pelvis
PelvisBase = MidAsi - normLength(PelvisRostralVector) * PelvisLength;
Segments.pelvis(4,:,:) = PelvisBase;

%% Thorax
% thorax_05_marker; thorax_06_marker; back_marker; thorax_16_marker
% Richtungsvektor rostral entlang der markern (nehme T10 anstatt T16 denn
% der Ruecken ist sehr biegsam)
ThoraxZVector       = SchulterRefMarker - Markers.back_marker; 
ThoraxLeftVector    =-normLength(cross(ThoraxZVector,ThoraxVentralVector,1));
ThoraxVentralVector = normLength(cross(ThoraxLeftVector,SpineVector,1));
DistanceBack2Spine  = median(point2LineDistance(Markers.back_marker,MidAsi,SpineVector),'omitnan');
ThoraxBase          = Markers.back_marker + DistanceBack2Spine * ThoraxVentralVector;
Segments.thorax     = createCS(ThoraxBase, ThoraxLeftVector, SpineVector,'yxz');
Segments.thorax(5,:,:) = MidShoulder;

LumbarZVector = ThoraxBase - MidAsi;
% get rostral, ventral and left (Y) vectors

%% lumbar spine segment
Segments.lumbar = createCS(MidAsi, LumbarZVector, HipVentralVector,'zyx');
Segments.lumbar(5,:,:) = ThoraxBase;

fprintf('%25s %5.0f%5.0f%5.0f\n','HipRostralVector',100*mean(PelvisRostralVector,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','HipVentralVector',100*mean(HipVentralVector,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','ThoraxVentralVector',100*mean(ThoraxVentralVector,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','ThoraxZVector',100*mean(ThoraxZVector,2)')
fprintf('%25s %5.0f%5.0f%5.0f : %f\n','n.ThoraxLatVector',100*mean(ThoraxLeftVector,2)',meanDistance(ThoraxLeftVector))
fprintf('%25s %5.0f%5.0f%5.0f : %f\n','n.ThoraxXVector',100*mean(ThoraxVentralVector,2)',meanDistance(ThoraxVentralVector))
fprintf('%25s %5.0f%5.0f%5.0f\n','SpineVector',100*mean(SpineVector,2)')
fprintf('%25s %5.0f\n','DistanceBack2Spine',100*mean(DistanceBack2Spine,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','MidAsi',100*mean(MidAsi,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','ThoraxBase',100*mean(ThoraxBase,2)')
fprintf('%25s %5.0f%5.0f%5.0f\n','MidShoulder',100*mean(MidShoulder,2)');
fprintf('%25s %5.0f%5.0f%5.0f\n','BackMarker',100*mean(Markers.back_marker,2)');
fprintf('%25s %5.0f\n','LumbarLen',100*meanDistance(ThoraxBase,MidAsi))
fprintf('%25s %5.0f\n','ThoraxLen',100*meanDistance(MidShoulder,ThoraxBase))
fprintf('%25s %5.0f\n','BodyLen',100*meanDistance(MidShoulder,MidAsi))

mean(Segments.thorax,3)
mean(Segments.pelvis,3)

%% head
if Present.head_left_marker && Present.head_right_marker && Present.head_front_marker && ...
        Present.left_cervical_2_marker && Present.right_cervical_2_marker
    HeadLeftVector = (Markers.head_left_marker       - Markers.head_right_marker)/2;
    MidHeadLR      = (Markers.head_left_marker       + Markers.head_right_marker)/2;
    MidNeck        = (Markers.left_cervical_2_marker + Markers.right_cervical_2_marker)/2;
    SkullBase      = (2*MidNeck + MidHeadLR)/3;
    RostralVector  = Markers.head_front_marker - MidHeadLR;
    Segments.skull = createCS(SkullBase, RostralVector, HeadLeftVector,'xzy');
    Segments.skull(5,:,:) = Markers.head_front_marker;
end

%% neck segment
NeckZVector = SkullBase - MidShoulder;
Segments.neck = createCS(MidShoulder, NeckZVector, ThoraxVentralVector,'zyx');
Segments.neck(5,:,:) = SkullBase;

%% Clavicles (these do not exist in horses, but are virtual and massless)
% estimate centre of shoulder joint (down from shoulder marker)
RightClav0    = MidShoulder;
LeftClav0     = MidShoulder;
if     Present.left_shoulder_marker      && Present.right_shoulder_marker
    rClavLongAxis = Markers.right_shoulder_marker - RightClav0;
    lClavLongAxis = Markers.left_shoulder_marker  - LeftClav0;
elseif Present.left_shoulder_back_marker && Present.right_shoulder_back_marker
    rClavLongAxis = Markers.right_shoulder_back_marker - RightClav0;
    lClavLongAxis = Markers.left_shoulder_back_marker  - LeftClav0;
else
    rClavLongAxis = nan(3,NSamples);
    lClavLongAxis = nan(3,NSamples);
end
% Coordinate system
% left: x-ventral, y-caudal, z-left
% right: x-dorsal, y-caudal, z-right
Segments.rClavicle = createCS(RightClav0,rClavLongAxis,-ThoraxZVector,'zxy');
Segments.lClavicle = createCS(LeftClav0, lClavLongAxis,-ThoraxZVector,'zxy');

%% right and left fore leg
% Anahme: jedes Gelenk liegt unter ein Marker, entlang der ThoraxXvector
JointRadius = 0.05; % 5 cm
rScapulaJoint  = Markers.right_shoulder_back_marker + JointRadius * ThoraxLeftVector;
rShoulderJoint = Markers.right_shoulder_marker      + JointRadius * ThoraxLeftVector;
rElbowJoint    = Markers.right_lat_epicondyle_marker+ JointRadius * ThoraxLeftVector; 
rWristJoint    = Markers.right_ulnar_wrist_marker   + JointRadius * ThoraxLeftVector;
rFingerJoint   = Markers.right_finger_marker        + JointRadius * ThoraxLeftVector;
lScapulaJoint  = Markers.left_shoulder_back_marker - JointRadius * ThoraxLeftVector;
lShoulderJoint = Markers.left_shoulder_marker      - JointRadius * ThoraxLeftVector;
lElbowJoint    = Markers.left_lat_epicondyle_marker- JointRadius * ThoraxLeftVector; 
lWristJoint    = Markers.left_ulnar_wrist_marker   - JointRadius * ThoraxLeftVector;
lFingerJoint   = Markers.left_finger_marker        - JointRadius * ThoraxLeftVector;
% Huf hat drei markern % left_front_hoof_dist_marker; left_front_hoof_b_marker
rFrHufZVector  = Markers.right_front_hoof_dist_marker - Markers.right_front_hoof_prox_marker;
rFrHufXVector  = normLength(cross(ThoraxLeftVector,rFrHufZVector,1));
rFrHufJoint    = Markers.right_front_hoof_prox_marker + JointRadius * rFrHufXVector;
lFrHufZVector  = Markers.left_front_hoof_dist_marker  - Markers.left_front_hoof_prox_marker;
lFrHufXVector  = normLength(cross(ThoraxLeftVector,lFrHufZVector,1));
lFrHufJoint    = Markers.left_front_hoof_prox_marker  + JointRadius * lFrHufXVector;

%% right and left Legs
% mean distance between hip joint markers
WidthHip    = meanDistance(Markers.right_femur_trochanter_marker,Markers.left_femur_trochanter_marker);
rHipJoint   = Markers.right_femur_trochanter_marker + PelvisYVector * WidthHip/4 + SpineVector * JointRadius;
lHipJoint   = Markers.left_femur_trochanter_marker  - PelvisYVector * WidthHip/4 + SpineVector * JointRadius;

rKneeJoint  = Markers.right_femur_lat_epi_marker + JointRadius * PelvisYVector;
rAnkleJoint = Markers.right_lat_malleolus_marker + JointRadius * PelvisYVector;
rToeJoint   = Markers.right_DM2_marker           + JointRadius * PelvisYVector;
lKneeJoint  = Markers.left_femur_lat_epi_marker  - JointRadius * PelvisYVector;
lAnkleJoint = Markers.left_lat_malleolus_marker  - JointRadius * PelvisYVector;
lToeJoint   = Markers.left_DM2_marker            - JointRadius * PelvisYVector;
% Huf hat drei markern 
rHiHufZVector = Markers.right_hind_hoof_dist_marker - Markers.right_hind_hoof_prox_marker;
rHiHufXVector = normLength(cross(ThoraxLeftVector,rHiHufZVector,1));
rHiHufJoint   = Markers.right_hind_hoof_prox_marker + JointRadius * rHiHufXVector;
lHiHufZVector = Markers.left_hind_hoof_dist_marker  - Markers.left_hind_hoof_prox_marker;
lHiHufXVector = normLength(cross(ThoraxLeftVector,lHiHufZVector,1));
lHiHufJoint   = Markers.left_hind_hoof_prox_marker  + JointRadius * lHiHufXVector;

% Koordinatensysteme
Segments.rScapula = createCS(rScapulaJoint, rShoulderJoint - rScapulaJoint,ThoraxLeftVector,'zyx');
Segments.lScapula = createCS(lScapulaJoint, lShoulderJoint - lScapulaJoint,ThoraxLeftVector,'zyx');
Segments.rHumerus = createCS(rShoulderJoint,rElbowJoint  - rShoulderJoint, ThoraxLeftVector,'zyx');
Segments.lHumerus = createCS(lShoulderJoint,lElbowJoint  - lShoulderJoint, ThoraxLeftVector,'zyx');
Segments.rRadius  = createCS(rElbowJoint,   rWristJoint  - rElbowJoint,    ThoraxLeftVector,'zyx');
Segments.lRadius  = createCS(lElbowJoint,   lWristJoint  - lElbowJoint,    ThoraxLeftVector,'zyx');
Segments.rHand    = createCS(rWristJoint,   rFingerJoint - rWristJoint,    ThoraxLeftVector,'zyx');
Segments.lHand    = createCS(lWristJoint,   lFingerJoint - lWristJoint,    ThoraxLeftVector,'zyx');
Segments.rFinger  = createCS(rFingerJoint,  rFrHufJoint  - rFingerJoint,   ThoraxLeftVector,'zyx');
Segments.lFinger  = createCS(lFingerJoint,  lFrHufJoint  - lFingerJoint,   ThoraxLeftVector,'zyx');
Segments.rFrHuf   = createCS(rFrHufJoint,   rFrHufZVector,       ThoraxLeftVector,'zyx');
Segments.lFrHuf   = createCS(lFrHufJoint,   lFrHufZVector,       ThoraxLeftVector,'zyx');
Segments.rFemur = createCS(rHipJoint,   rKneeJoint  -rHipJoint,  PelvisYVector,'zxy');
Segments.lFemur = createCS(lHipJoint,   lKneeJoint  -lHipJoint,  PelvisYVector,'zxy');
Segments.rTibia = createCS(rKneeJoint,  rAnkleJoint -rKneeJoint, PelvisYVector,'zxy');
Segments.lTibia = createCS(lKneeJoint,  lAnkleJoint -lKneeJoint, PelvisYVector,'zxy');
Segments.rFoot  = createCS(rAnkleJoint, rToeJoint - rAnkleJoint, PelvisYVector,'zxy');
Segments.lFoot  = createCS(lAnkleJoint, lToeJoint - lAnkleJoint, PelvisYVector,'zxy');
Segments.rToes  = createCS(rToeJoint,   rHiHufJoint - rToeJoint, PelvisYVector,'zxy');
Segments.lToes  = createCS(lToeJoint,   lHiHufJoint - lToeJoint, PelvisYVector,'zxy');
% Segments.rHiHuf = cs_create(rHiHufJoint, rFootVector, HipYVector,'zxy');
% Segments.lHiHuf = cs_create(lHiHufJoint, lFootVector, HipYVector,'zxy');

%% define endpoints of distal fore- and hind leg segments
if Present.right_finger_marker
    Segments.rFinger(5,:,:) = Markers.right_finger_marker;
end
if Present.left_finger_marker
    Segments.lFinger(5,:,:) = Markers.left_finger_marker;
end
if Present.right_DM2_marker
    Segments.rToes(5,:,:) = Markers.right_DM2_marker;
    %Segments.rHiHuf(5,:,:) = Markers.right_hind_hoof_b_marker;
end
if Present.left_DM2_marker
    Segments.lToes(5,:,:) = Markers.left_DM2_marker;
    %Segments.lHiHuf(5,:,:) = Markers.left_hind_hoof_b_marker;
end
end
