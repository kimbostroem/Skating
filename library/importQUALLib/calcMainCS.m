function Segments = calcMainCS(Markers, Present, SubjectBody, Info)
%% calculate the time series of the segment coordinate systems from the marker data
% A coordinate system (CS) has for each time sample
% - a rotation matrix for the global orientation of the segment
% - a vector to the base of the segment (global coordinates)
% - optionally a vector to the end of the segment (global coordinates)
% - size: 4 x 3 x NSamples or 5 x 3 x NSamples
% Coordinates: 
% - z always distal, starting from the hips
% - y to the left; arms in t-pose: up
% - x ventral in pelvis, trunk, head; dorsal in legs; left arm dorsal; right arm ventral
% Usage: 
% - calcMainCS is called by importQUAL and importBodyFromQUAL
% - importBodyFromQUAL : subject body is empty; joint positions of pelvis, thorax and 
%     clavicle/shoulder are based on marker positions and relative scaling of def_body positions.
% - importQUAL : positions are based entirely on the user_body
% 
% Position of hip joints in Pelvis: 
% Fiorentino et al. (2016) Accuracy of functional and predictive methods to calculate the hip joint center in young non-pathologic asymptomatic adults with dual fluoroscopy as a reference standard. Ann. Biomed. Eng. doi: 10.1007/s10439-015-1522-1.
% Harrington et al. (2007) Prediction of the hip joint centre in adults, children, and patients with cerebral palsy based on magnetic resonance imaging. J. Biomech. doi: 10.1016/j.jbiomech.2006.02.003.
%
% SYNTAX
%   Segments = calcMainCS(Markers, Present, SubjectBody, Info)
%
% INPUT
%     Markers     (struct of 3 x N) time series for each marker that is present
%     Present     (struct of logical) for each marker: is it present
%     SubjectBody (struct) def_body tables including user-specific parameters replacing general and segments tables of def_body
%     Info        (stuct) information, e.g. NSintQTM = number of samples; Tail for interpolation
%
% OUTPUT
%    Segments     (struct) Coordinate systems of the segments
% 
% EXAMPLE
% Using the default body, when processing the Body-file
%     Segments = calcMainCS(Markers,Present,[], Info);
% Using the Subject-body, which contains the correct segment lengths and the positions 
% of the translated joints:
%     Segments = calcMainCS(Markers,Present,SubjectBody, Info);
%
% Local functions: calcKneeElbAxis, vector_interpol
%
% See also: importQUAL, createCS
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet and Lena Kloock
% version 230215 (MdL & MvdH) FIXED thorax translation and thorax top position; various further corrections
% after review
% version 230217 (MdL) fix in vector_interpol
% version 230217 (MdL) added estimate of hip position from rotation of knee marker: at present not
% yet ready (so commented out); define BodyHeight at Init Section
% version 230221 (MdL) calculate subject height if it is not yet defined; update comments; comment-out currently non-used code
% version 230227 (MdL) some fine tuning of head scaling
% version 230306 (MdL) FIX check for exist of variable SkullBase

%% Init

% Do not continue if no markers are present
Segments = [];
if isempty(Present)
    return;
end
% convert the Present structure to a logical array
IsPresent = cell2mat(struct2cell(Present));
if ~any(IsPresent)
    return
end

% Parameters
X=1;Y=2;Z=3; %#ok<NASGU> 


%% DEFAULT VALUES POSSIBLY IN PREFS FILE?
%% (->try to read this value from the T-pose and store in user body)

% Default lengths are only used if median marker is missing
Kneewidth         = 0.08; % (m)
Anklewidth        = 0.08; % (m)
MarkerRadius      = 0.01; % (m)
% These angles are used for handling extended knee, elbow and ankle angles
ArmTransit        = [160 170];  % transition range for switch between arm plane
ArmExtendThreshold= 160;  % vector_interpol threshold for interpolating the joint vector when extended
%  Tail [samples] for vector_interpol (for getting elbow and knee axis)
Tail_Sec          = 0.15;
Tail              = floor(Tail_Sec * Info.QTMFreqSRint/2);   %
NSamples          = Info.NSintQTM;

% Check if the SubjectBody has been defined
if isempty(SubjectBody)
    % If the subject body is not yet defined, the leg length is needed to estimate the hips and the
    % trunk length to estimate the clavicle positions
    BodyTables          = loadBody('def_body');
    IsSubjectBodyLoaded = false;
    % Estimate the stable pose period for estimating the length measures
    StablePeriod        = getStablePeriod(Markers,Present,Info,false);
else
    BodyTables          = SubjectBody;
    IsSubjectBodyLoaded = true;
end

% Extract the default local offsets from the def_body tables
PelvisID    = 20;
L5_ID       = 25;
LumbarIDs   = 21:25;
ThoraxID    = 2;
LFemurID    = 27;
RFemurID    = 26;
LClavicleID = 13;
RClavicleID = 12;
LHumerusID  = 15;
%SkullID     =  1;
PelvisNo    = BodyTables.segments.id == PelvisID;
ThoraxNo    = BodyTables.segments.id == ThoraxID;
LClavicleNo = BodyTables.segments.id == LClavicleID;
%SkullNo     = BodyTables.segments.id == SkullID;
LumbarNos   = logical(sum(BodyTables.segments.id == LumbarIDs, 2));
L5_NextNo       = BodyTables.segments.next(PelvisNo,:) == L5_ID;    % which entry of nextpos gives the L5 joint
LHipNextNo      = BodyTables.segments.next(PelvisNo,:) == LFemurID; % which entry of nextpos gives the left hip
RHipNextNo      = BodyTables.segments.next(PelvisNo,:) == RFemurID; % which entry of nextpos gives the right hip
LShoulderNextNo = BodyTables.segments.next(LClavicleNo,:) == LHumerusID;
ClavNextNo      = logical(sum(BodyTables.segments.next(ThoraxNo,:) == [LClavicleID; RClavicleID]));
BodyPelvisSize   = BodyTables.segments.size(PelvisNo);
BodyThoraxLength = BodyTables.segments.size(ThoraxNo);
% NOTE (2023-02-25): The value of 0.117 in def_body is clearly too much. 
BodySkullSize    = .0889; %BodyTables.segments.size(SkullNo);
BodyLShoulder    = squeeze(BodyTables.segments.nextPos(LClavicleNo,LShoulderNextNo,:));
BodyDefaultSize  = BodyTables.general.height;

% Init the Segments structure
Segments.skull    = nan(4,3,NSamples); % 3x3 orientation & 3x1 base vector
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
ThoraxTopMarker   = nan(3,NSamples);
ThoraxLowMarker   = nan(3,NSamples);


%% Define head
% We need the head first, to estimate the body height

DoCalcHeadFront = ~Present.head_front_marker && Present.head_front_right_marker && Present.head_front_left_marker;
IsHeadLeftRight = Present.head_left_marker && Present.head_right_marker;
IsHeadOblique   = ~IsHeadLeftRight && ...
    Present.head_front_left_marker && Present.head_front_right_marker && ...
    Present.head_back_left_marker  && Present.head_back_right_marker;
if DoCalcHeadFront
    Markers.head_front_marker = (Markers.head_front_right_marker + Markers.head_front_left_marker)/2;
    Present.head_front_marker = true;
end
% vector in Y direction with length 1 radius
if IsHeadLeftRight
    HeadLeftVector = (Markers.head_left_marker - Markers.head_right_marker)/2;
    HeadCentre     = (Markers.head_left_marker + Markers.head_right_marker)/2;
elseif IsHeadOblique
    % correct the length of the left vector assuming that the markers sit at 45 degrees
    HeadLeftVector = (...
        Markers.head_front_left_marker - Markers.head_front_right_marker + ...
        Markers.head_back_left_marker  - Markers.head_back_right_marker   )/(4 * sind(45));
    HeadCentre     = (...
        Markers.head_front_left_marker + Markers.head_front_right_marker + ...
        Markers.head_back_left_marker  + Markers.head_back_right_marker   )/4;    
elseif Present.head_front_left_marker && Present.head_front_right_marker
    % correct the length of the left vector assuming that the markers sit at 45 degrees
    HeadLeftVector = (Markers.head_front_left_marker - Markers.head_front_right_marker)/(2 * sind(45));
    HeadCentre     = (Markers.head_front_left_marker + Markers.head_front_right_marker)/2;
end
% vector in X direction (front) if possible
if IsHeadOblique
    HeadXVector    = normLength( ...
        Markers.head_front_left_marker + Markers.head_front_right_marker - ...
        Markers.head_back_left_marker  - Markers.head_back_right_marker   );
elseif Present.head_front_marker && exist('HeadCentre','var')
    HeadXVector    = normLength(Markers.head_front_marker - HeadCentre);
end
% The vertical skull vector depending on which variables have been computed
if exist('HeadLeftVector','var') && exist('HeadXVector','var')
    HeadZVector    = cross(HeadXVector, normLength(HeadLeftVector));
elseif Present.head_top_marker && exist('HeadCentre','var')
    HeadZVector    = normLength(Markers.head_top_marker   - HeadCentre);
end
% The skull Coordinate System (CS)
if     exist('HeadCentre','var') && exist('HeadXVector','var') && exist('HeadLeftVector','var')
    Segments.skull = createCS(HeadCentre, HeadXVector, HeadLeftVector,'xzy');
    HeadZVector    = squeeze(Segments.skull(3,:,:));
elseif exist('HeadCentre','var') && exist('HeadZVector','var') && exist('HeadLeftVector','var')
    Segments.skull = createCS(HeadCentre, HeadZVector,-HeadLeftVector,'zxy');
end
% Skull height
if exist('HeadLeftVector','var')
    SkullHeight = mean(vecnorm(HeadLeftVector, 2), 'omitnan') * 1.8;
else % use the default value for the skull height
    SkullHeight = BodySkullSize * BodyDefaultSize;
end
% The top and base of the head
if Present.head_top_marker
    SkullTop  = Markers.head_top_marker;
    SkullBase = SkullTop - HeadZVector * SkullHeight;
elseif exist('HeadCentre','var') && exist('HeadZVector','var')
    % If the top marker is absent, then use the marker-based center:
    % this is significantly higher than the base of the skull (probably 1/3 - 1/2 times the height)
    SkullTop  = HeadCentre + HeadZVector * SkullHeight*2/3;
    SkullBase = HeadCentre - HeadZVector * SkullHeight/3;
end
% Re-define base and end of skull segment
if exist('SkullBase','var') && exist('SkullTop','var')
    Segments.skull(4,:,:) = SkullBase;
    Segments.skull(5,:,:) = SkullTop;
end


%% Define body height
% The body height is central for scaling the segment sizes
if IsSubjectBodyLoaded
    SubjectHeight       = BodyTables.general.height;
else
    Segments.bottom_position(4,:,:) = Markers.bottom_position;
    DefaultHeight       = BodyTables.general.height;
    SubjectHeight       = estimateBodyHeightFromMarkers(Segments, DefaultHeight, StablePeriod);
end


%% Define pelvis
% The pelvis (at mid-asi) and RBM (at mid-hip) are the base from which the body is built

% select the reference markers
PelvisPsiMarkersPresent    = Present.right_psi_marker & Present.left_psi_marker; 
PelvisBackMarkersPresent   = PelvisPsiMarkersPresent  | Present.sacrum_marker;
PelvisFrontMarkersPresent  = Present.left_asi_marker  & Present.right_asi_marker;
PelvisTopMarkersPresent    = Present.left_mid_crista_iliaca_marker && Present.right_mid_crista_iliaca_marker;
UsePelvisTopInsteadOfBack  = PelvisTopMarkersPresent  &&  PelvisFrontMarkersPresent && ~PelvisBackMarkersPresent;
UsePelvisTopInsteadOfAsi   = PelvisTopMarkersPresent  && ~PelvisFrontMarkersPresent &&  PelvisBackMarkersPresent;
% Select reference markers, depending on which ones are present
if PelvisPsiMarkersPresent
    % midpoint hip on the posterior side: the sacrum
    PelvisBackMarker = (Markers.right_psi_marker + Markers.left_psi_marker)/2;
elseif Present.sacrum_marker
    % the sacrum
    PelvisBackMarker = Markers.sacrum_marker;
elseif Present.lumbar_5_marker
    % the lumbar marker sits about 5 cm above the back marker, though
    PelvisBackMarker = Markers.lumbar_5_marker;
elseif UsePelvisTopInsteadOfBack
    % in case nothing else works
    PelvisBackMarker = (Markers.right_mid_crista_iliaca_marker + Markers.left_mid_crista_iliaca_marker)/2;
else
    PelvisBackMarker = nan(3,NSamples); 
end
% When the trunk is bent (e.g. sitting, biking) the ASI markers are invisible and pelvis top markers
% must be used instead. The pelvis origin should be shifted down (-Z) accordingly. 
PelvisZOffset = 0.08; % 8 cm: default value
if UsePelvisTopInsteadOfAsi && Present.lumbar_5_marker
    % if lumbar-5 and other back markers are present, the vertical shift can be calculated:
    PelvisZOffset = mean(vecnorm(PelvisBackMarker - Markers.lumbar_5_marker, 2), 'omitnan');
    % The top markers are higher: use L5 marker to get the pelvis orientation right
    PelvisBackMarker = Markers.lumbar_5_marker;
end

% midpoint on anterior side and right to left vector
if PelvisFrontMarkersPresent
    LFrontPelvis     = Markers.left_asi_marker;
    RFrontPelvis     = Markers.right_asi_marker;
    Mid_asi          = (LFrontPelvis + RFrontPelvis)/2; % midpoint hip on the anterior side
    PelvisLeftVector =  LFrontPelvis - RFrontPelvis;
elseif PelvisTopMarkersPresent
    LFrontPelvis     = Markers.left_mid_crista_iliaca_marker;
    RFrontPelvis     = Markers.right_mid_crista_iliaca_marker;
    Mid_asi          = (LFrontPelvis + RFrontPelvis)/2; 
    PelvisLeftVector = (LFrontPelvis - RFrontPelvis) * 0.9; % *0.9: top markers are slightly further apart than Asi
else % head of femur markers (not so good)
    LFrontPelvis     = Markers.left_femur_trochanter_marker;
    RFrontPelvis     = Markers.right_femur_trochanter_marker;
    Mid_asi          = (LFrontPelvis + RFrontPelvis)/2;
    PelvisLeftVector = (LFrontPelvis - RFrontPelvis) * 0.85; % *0.85: PSI markers are even further apart
end
PelvisYVector = normLength(PelvisLeftVector);

% back to front vector
PelvisVentralVector = Mid_asi - PelvisBackMarker;
% correct the back marker if it was estimated from the pelvis top markers
if UsePelvisTopInsteadOfBack
    PelvisBackMarker = Mid_asi - PelvisVentralVector;
    PelvisVentralVector = 2 * PelvisVentralVector;
elseif UsePelvisTopInsteadOfAsi
    PelvisVentralVector = 2 * PelvisVentralVector;
end

% Coordinate system and offsets of the joints
Segments.pelvis = createCS(Mid_asi, PelvisLeftVector, -PelvisVentralVector,'yzx');
% subtract the vertical offset in case the asi markers could not be used
if UsePelvisTopInsteadOfAsi
    Segments.pelvis(4,:,:) = Segments.pelvis(4,:,:) - Segments.pelvis(3,:,:) * PelvisZOffset;
end

% % If the subject body is not yet loaded, the position of the Hips must be estimated here
% % Note, that, at present, only the pelvis width can be scaled, whereas the height (z) and depth (x) 
% % of the L5 and hip joints are coupled to the estimated hip distance
% if ~IsSubjectBodyLoaded
% 
%     %% Hip position in Pelvis (Harrington et al 2007, see header; note, that Y and Z axes are swapped)
%     % % THIS DOES NOT RESULT IN A BETTER ESTIMATE OF THE HIP POSITION
%     % % calculate lengths and distances to pelvic joints
%     % PelvisWidth = mean(vecnorm(PelvisLeftVector(   :,StablePeriod), 2), 'omitnan');
%     % PelvisDepth = mean(vecnorm(PelvisVentralVector(:,StablePeriod), 2), 'omitnan') - 2 * MarkerRadius;
%     % 
%     % % estimate leg length: 
%     % % either the complete length ...
%     % LeftHip2AnkleLength  = vecnorm(LFrontPelvis - Markers.left_lat_malleolus_marker,2);
%     % RightHip2AnkleLength = vecnorm(RFrontPelvis - Markers.right_lat_malleolus_marker,2);
%     % Hip2AnkleLength = mean((LeftHip2AnkleLength(StablePeriod) + RightHip2AnkleLength(StablePeriod)) / 2, 'omitnan');
%     % % ... or the sum of upper and lower leg segments
%     % LeftFemurLength  = vecnorm(Markers.left_femur_lat_epi_marker  - LFrontPelvis,                      2);
%     % LeftTibiaLength  = vecnorm(Markers.left_femur_lat_epi_marker  - Markers.left_lat_malleolus_marker, 2);
%     % RightFemurLength = vecnorm(Markers.right_femur_lat_epi_marker - RFrontPelvis,                      2);
%     % RightTibiaLength = vecnorm(Markers.right_femur_lat_epi_marker - Markers.right_lat_malleolus_marker,2);
%     % FemurPlusTibiaLength = (LeftFemurLength + LeftTibiaLength + RightFemurLength + RightTibiaLength) / 2;
%     % FemurPlusTibiaLength = mean(FemurPlusTibiaLength(StablePeriod), 'omitnan');
%     % % maximum
%     % LegLength = max(Hip2AnkleLength, FemurPlusTibiaLength,'omitnan');
%     % if isnan(LegLength)
%     %     LegLength = 1;
%     % end
%     % 
%     % HipX =-0.24 * PelvisDepth - 0.0099;
%     % HipY = 0.28 * PelvisDepth + 0.16 * PelvisWidth + 0.0079;
%     % HipZ =-0.16 * PelvisWidth - 0.04 * LegLength - 0.0071;
%     % LHip_Offset = [HipX  HipY HipZ];
%     % RHip_Offset = [HipX -HipY HipZ];
%     % % The L5 joint is not based on the literature, since the pelvis is scaled only to the hip-distance
%     % %L5X  =-0.65 * PelvisDepth;
%     % %L5Z  = 0.05;
%     % % X= at about 2/3 toward the backmarker (Harrington et al); Y= middle; Z= about 3 cm above the pelvis plane
%     % %L5_Offset = [L5X 0 L5Z];
%     % 
%     % S=100;
%     % fprintf('L %5.1f %5.1f %5.1f\n',LHip_Offset(1)*S,LHip_Offset(2)*S,LHip_Offset(3)*S);
%     % fprintf('R %5.1f %5.1f %5.1f\n',RHip_Offset(1)*S,RHip_Offset(2)*S,RHip_Offset(3)*S);
%     
%     %% Estimate hip joint from the rotation of the knee marker about the joint
%     % IsLKneeMarkers = Present.left_femur_med_epi_marker  && Present.left_femur_lat_epi_marker;
%     % IsRKneeMarkers = Present.right_femur_med_epi_marker && Present.right_femur_lat_epi_marker;
%     % 
%     % if IsLKneeMarkers
%     %     LKneeMarker  = (Markers.left_femur_med_epi_marker +  Markers.left_femur_lat_epi_marker )/2;
%     % else
%     %     LKneeMarker  = Markers.left_femur_lat_epi_marker;
%     % end
%     % if IsRKneeMarkers
%     %     RKneeMarker  = (Markers.right_femur_med_epi_marker + Markers.right_femur_lat_epi_marker )/2;
%     % else
%     %     RKneeMarker  = Markers.right_femur_lat_epi_marker;
%     % end
%     % LKneeLocal = getLocalPosOnSegment(Segments.pelvis,LKneeMarker);
%     % RKneeLocal = getLocalPosOnSegment(Segments.pelvis,RKneeMarker);
%     % [LHip_Offset,LRadius,LStd,LRadius_t]=sphereFit(LKneeLocal');
%     % [RHip_Offset,RRadius,RStd,RRadius_t]=sphereFit(RKneeLocal');
%     % LRng = (max(LKneeMarker,[],2) - min(LKneeMarker,[],2)) ./ LRadius;
%     % RRng = (max(RKneeMarker,[],2) - min(RKneeMarker,[],2)) ./ RRadius;
%     % % figure;
%     % % subplot(1,2,1);hold on;plot(LKneeLocal');plot(LRadius_t);title('calcMainCS :','left hip');legend({'X' 'Y' 'Z' 'R'})
%     % % subplot(1,2,2);hold on;plot(RKneeLocal');plot(RRadius_t);title('knee marker movement in hip coordinates','right hip');
%     % S=100;
%     % fprintf('L %5.1f %5.1f %5.1f (R=%5.1f Sdev=%5.1f) %5.2f %5.2f %5.2f\n',LHip_Offset(1)*S,LHip_Offset(2)*S,LHip_Offset(3)*S,LRadius*S,LStd*S,LRng(1),LRng(2),LRng(3));
%     % fprintf('R %5.1f %5.1f %5.1f (R=%5.1f Sdev=%5.1f) %5.2f %5.2f %5.2f\n',RHip_Offset(1)*S,RHip_Offset(2)*S,RHip_Offset(3)*S,RRadius*S,RStd*S,RRng(1),RRng(2),RRng(3));
%     % 
%     % LHip_Offset(2) = nan;
%     % RHip_Offset(2) = nan;
%     % % lateral offset is not always good: to use the default:
%     % LHip_Offset(2) = squeeze(BodyTables.segments.nextPos(PelvisNo,LHipNextNo,2)) * BodyPelvisSize * SubjectHeight;
%     % RHip_Offset(2) = squeeze(BodyTables.segments.nextPos(PelvisNo,RHipNextNo,2)) * BodyPelvisSize * SubjectHeight;
% 
%     %% use the default position
%     LHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,LHipNextNo,:)) * BodyPelvisSize * SubjectHeight;
%     RHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,RHipNextNo,:)) * BodyPelvisSize * SubjectHeight;
% 
% else 
%     % the local offsets according to def_body
%     LHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,LHipNextNo,:)) * BodyPelvisSize * SubjectHeight;
%     RHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,RHipNextNo,:)) * BodyPelvisSize * SubjectHeight;
% end
% the local offsets according to def_body
LHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,LHipNextNo,:)) * BodyPelvisSize * SubjectHeight;
RHip_Offset = squeeze(BodyTables.segments.nextPos(PelvisNo,RHipNextNo,:)) * BodyPelvisSize * SubjectHeight;

% Joint positions in time
LHipJoint   = getGlobalPosFromSegment(Segments.pelvis, LHip_Offset);
RHipJoint   = getGlobalPosFromSegment(Segments.pelvis, RHip_Offset);
MidHip = (LHipJoint + RHipJoint) / 2;
fprintf('calcMainCS:  ===============\nBodyPelvisSize %.1fcm, SubjectHeight %.1fcm\n',100*BodyPelvisSize, 100*SubjectHeight);
fprintf('hip  dist = %.1fcm\n',100 * mean(vecnorm(LHipJoint - RHipJoint)));
KneeDist = vecnorm(Markers.left_asi_marker - Markers.right_asi_marker);
fprintf('knee dist = %.1fcm (%.1f -- %.1f)\n',100 * mean(KneeDist),100 * min(KneeDist),100 * max(KneeDist));


%% Define thorax
% The thorax defines the proximal lumbar spine, the base of the clavicles and the base of the
% cervical spine
% - check which markers are present
% - define the orientation : 
%     - Z (rostral) from the markers on the back
%     - ventral from back and ventral markers
% define the size from the default size (?)
% define the position such that the clavicle marker is at the correct position

%  Check which markers are present
UpperThoraxMarkers   = {'neck_marker','thorax_2_marker'}; % rostral to caudal order 
LowerThoraxMarkers   = {'thorax_12_marker','back_marker'}; % caudal to rostral order
MiddleThoraxMarkers  = {'thorax_05_marker','thorax_7_marker'};
BackMarkers          = [MiddleThoraxMarkers,LowerThoraxMarkers];
VentralThoraxMarkers = {'clavicula_marker','sternum_marker'}; % rostral to caudal order 
UpperThoraxMarkerNos   = matches(BodyTables.markers.name,UpperThoraxMarkers);
LowerThoraxMarkerNos   = matches(BodyTables.markers.name,LowerThoraxMarkers);
BackMarkerNos          = matches(BodyTables.markers.name,BackMarkers);
PresentThoraxTopMarker = any(IsPresent(UpperThoraxMarkerNos));
PresentThoraxLowMarker = any(IsPresent(LowerThoraxMarkerNos));
VentralThoraxMarkerNos = matches(BodyTables.markers.name,VentralThoraxMarkers);
CanCalcThoraxVentralVector    = any(IsPresent(BackMarkerNos)) && any(IsPresent(VentralThoraxMarkerNos));
PresentShoulderMarkers        = Present.left_shoulder_marker & Present.right_shoulder_marker;
PresentThoraxDorsalLRMarkers  = (Present.right_back_marker && Present.left_back_marker);
PresentThoraxVentralLRMarkers = (Present.right_chest_marker && Present.left_chest_marker);
PresentThoraxLRMarkers        = PresentThoraxVentralLRMarkers || PresentThoraxDorsalLRMarkers;

% Define the Z-axis
% Select thorax top marker (rostral) and lowest marker (caudal)
% ThoraxTop_Z_offset = 0; % [m] estimated missing length for thorax % LK: can be deleted?
% ThoraxLow_Z_offset = 0; % [m] estimated missing length for thorax % LK: can be deleted?
if PresentThoraxTopMarker
    % select the most rostral thorax marker
    for i=1:length(UpperThoraxMarkers)
        if Present.(UpperThoraxMarkers{i})
            ThoraxTopMarker = Markers.(UpperThoraxMarkers{i});
            % ThoraxTopMarkerName = UpperThoraxMarkers{i};
            % if matches(ThoraxTopMarkerName,'thorax_2_marker')
            %     ThoraxTop_Z_offset = 0.05; % [m] estimated missing length for thorax
            % end
            break;
        end
    end
end
if PresentThoraxLowMarker
    % select the lowest (most caudal) thorax marker
    for i=1:length(LowerThoraxMarkers)
        if Present.(LowerThoraxMarkers{i})
            ThoraxLowMarker = Markers.(LowerThoraxMarkers{i});
            % ThoraxLowMarkerName = LowerThoraxMarkers{i};
            % if matches(ThoraxLowMarkerName,'thorax_12_marker')
            %     % T12 is the lowest possible thorax marker
            % elseif matches(ThoraxLowMarkerName,'back_marker')
            %     ThoraxLow_Z_offset = 0.05; % [m] add a "small" value if the lowest marker sits higher
            % else
            %     warning('non-defined case of ThoraxLowMarkerName = "%s"',ThoraxLowMarkerName);
            % end
            break;
        end
    end
else
    % if none of the back markers is present, estimate it as the mean between neck and sacrum
    ThoraxLowMarker = (ThoraxTopMarker + PelvisBackMarker)/2;
end
CanCalcThoraxVector = PresentThoraxTopMarker && ~all(isnan(ThoraxLowMarker),'all');
if CanCalcThoraxVector
    ThoraxRostralVector = ThoraxTopMarker - ThoraxLowMarker; % up
else
    ThoraxRostralVector = nan(3,NSamples);
end

% Define the ventral vector
if CanCalcThoraxVentralVector
    % select a dorsal thorax marker
    for i=1:length(BackMarkers)
        if Present.(BackMarkers{i})
            BackMarker = Markers.(BackMarkers{i});
            break;
        end
    end
    % select a ventral marker
    for i=1:length(VentralThoraxMarkers)
        if Present.(VentralThoraxMarkers{i})
            FrontMarker = Markers.(VentralThoraxMarkers{i});
            break;
        end
    end
    ThoraxVentralVector = FrontMarker - BackMarker; % ventral
else
    ThoraxVentralVector = nan(3,NSamples);
end
ThoraxZVector = normLength(ThoraxRostralVector);

% Define the shouder- mid-point and vertical offset from the clavicle
if ~IsSubjectBodyLoaded
    ShoulderMarkerDistance = vecnorm(Markers.left_shoulder_marker - Markers.right_shoulder_marker, 2);
    ShoulderMarkerDistance = mean(ShoulderMarkerDistance(StablePeriod), 'omitnan');
    ShoulderVertOff = BodyLShoulder(Y) * ShoulderMarkerDistance / 2;
else
    DefaultLClavSize = BodyTables.segments.size(LClavicleNo) * SubjectHeight;
    ShoulderVertOff = BodyLShoulder(Y) * DefaultLClavSize;
end

% Estimate the position of the shoulder joint from the thorax axis and the 
% vertical offset from the shoulder marker
LShoulderJoint = Markers.left_shoulder_marker  + ThoraxZVector * ShoulderVertOff;
RShoulderJoint = Markers.right_shoulder_marker + ThoraxZVector * ShoulderVertOff;
ShoulderVector = Markers.left_shoulder_marker - Markers.right_shoulder_marker;
% Mid-shoulder is needed in case thorax-top marker and/or clavicle marker are absent
MidShoulder = (LShoulderJoint + RShoulderJoint)/2;

% Translate the thorax axis to the vertebral column
ThoraxBase = ThoraxLowMarker;
ThoraxTop  = ThoraxTopMarker;
if CanCalcThoraxVentralVector
    % get the ventral (X) vector that is orthogonal
    ThoraxLeftVector = cross(ThoraxZVector,ThoraxVentralVector,1);
    ThoraxXVector = normLength(cross(ThoraxLeftVector,ThoraxZVector,1));
%     SpineDepth = mean(vecnorm(ThoraxVentralVector, 2), 'omitnan') / 3;
%     SpineVentralTranslation = SpineDepth * ThoraxXVector;
%     ThoraxBase = ThoraxBase + SpineVentralTranslation;
%     ThoraxTop  = ThoraxTop  + SpineVentralTranslation;
elseif ~PresentThoraxTopMarker && PresentShoulderMarkers
    ThoraxTop  = MidShoulder;
else
    ThoraxTop  = ThoraxTopMarker;
end
% hierarchical decision for how to estimate the thorax right to left vector
if CanCalcThoraxVentralVector
    ThoraxYVector = cross(ThoraxZVector,ThoraxXVector);
elseif PresentThoraxLRMarkers
    ThoraxDorsalYVector  = (Markers.left_back_marker  - Markers.right_back_marker);
    ThoraxVentralYVector = (Markers.left_chest_marker - Markers.right_chest_marker);
    if     ~ThoraxDorsalYVector,  ThoraxYVector = ThoraxVentralYVector;
    elseif ~ThoraxVentralYVector, ThoraxYVector = ThoraxDorsalYVector;
    else,                         ThoraxYVector = ThoraxDorsalYVector + ThoraxVentralYVector;
    end
elseif PresentShoulderMarkers
    ThoraxYVector = ShoulderVector;
else
    ThoraxYVector = PelvisYVector;
end
% normalize the length
ThoraxYVector = normLength(ThoraxYVector);

% We can now obtain the ventral vector if this was not yet done
if ~CanCalcThoraxVentralVector
    ThoraxXVector = cross(ThoraxYVector,ThoraxZVector);
end

% Compensate for the possible rostral and caudal offsets
% ThoraxTop   = ThoraxTop  + ThoraxTop_Z_offset * ThoraxZVector;
% ThoraxBase  = ThoraxBase - ThoraxLow_Z_offset * ThoraxZVector;

% calculate the thorax CS
Segments.thorax = createCS(ThoraxBase, ThoraxRostralVector, ThoraxVentralVector,'zyx');
Segments.thorax(5,:,:) = ThoraxTop;

% Define the thorax size
DefaultLumbarLength = sum(BodyTables.segments.size(LumbarNos) .* BodyTables.segments.nextPos(LumbarNos,1,3));
DefaultL5Local = BodyPelvisSize * squeeze(BodyTables.segments.nextPos(PelvisNo,L5_NextNo,:));
DefaultMidHipLocal = squeeze(BodyTables.segments.nextPos(PelvisNo,LHipNextNo,:));
DefaultMidHipLocal(Y) = 0; % we need the mid-hip position so y=0
Default_HipToL5 = DefaultL5Local - DefaultMidHipLocal;
DefaultTrunkLength = Default_HipToL5(Z) + DefaultLumbarLength + BodyThoraxLength;
TrunkLength  = max(vecnorm(MidHip - ThoraxTop, 2),[], 'omitnan');
if ~IsSubjectBodyLoaded
    ThoraxLength = BodyThoraxLength * TrunkLength / DefaultTrunkLength;
else
    ThoraxLength = BodyThoraxLength * SubjectHeight;
end


%% Define clavicles 
% ... and shift the Thorax to match the clavicle base

% estimate the origin of the clavicula from the clavicula marker if present
if Present.clavicula_marker
    ClavicleBase = Markers.clavicula_marker - ThoraxXVector * 2 * MarkerRadius;
elseif Present.sternum_marker % the Qualisys Sportsskeleton requires sternum marker, which can be attached at the clavicula base
    ClavicleBase = Markers.sternum_marker - ThoraxXVector * 2 * MarkerRadius;
else % if no clavicle marker present, assume midpoint between the shoulders
    ClavicleBase = MidShoulder - ThoraxXVector * 2 * MarkerRadius; 
end

% Shift the thorax such that the clavicle base is at the correct position
% Clavicle base in thorax coordinates from user body (mean of left and right clavicles)
BodyClavBaseLocal = squeeze(mean(    BodyTables.segments.nextPos(ThoraxNo,ClavNextNo,:))) * ThoraxLength;
BodyClav_Y_Offset = squeeze(mean(abs(BodyTables.segments.nextPos(ThoraxNo,ClavNextNo,2)))) * ThoraxLength;
if Present.clavicula_marker
    ClavBaseLocal = mean(getLocalPosOnSegment(Segments.thorax,ClavicleBase), 2, 'omitnan');
    ThoraxShift = BodyClavBaseLocal - ClavBaseLocal;
    ThoraxBase = getGlobalPosFromSegment(Segments.thorax,-ThoraxShift);
    ThoraxTop  = getGlobalPosFromSegment(Segments.thorax,-ThoraxShift,'end');
    %ThoraxTop = ThoraxTop + ThoraxBase - squeeze(Segments.thorax(4,:,:));
    Segments.thorax(4,:,:) = ThoraxBase;
    Segments.thorax(5,:,:) = ThoraxTop;
end


LClavicleBase = ClavicleBase + BodyClav_Y_Offset * ThoraxYVector; 
RClavicleBase = ClavicleBase - BodyClav_Y_Offset * ThoraxYVector;
LClavicleEnd  = Markers.left_shoulder_marker  - ThoraxZVector * 2 * MarkerRadius;
RClavicleEnd  = Markers.right_shoulder_marker - ThoraxZVector * 2 * MarkerRadius;
% figure;plot(sqrt(sum((RClavicleEnd-RClavicleBase).^2)))
% figure;plot(sqrt(sum((Markers.clavicula_marker-Markers.left_shoulder_marker).^2)))
% figure;hold on; plot(Markers.clavicula_marker(3,:)); plot(Markers.sternum_marker(3,:))

% estimate long axis of clavicle (Z)
LClavLongAxis  = LClavicleEnd - LClavicleBase;
RClavLongAxis  = RClavicleEnd - RClavicleBase;
% Coordinate system
% left: x-front, y-down, z-left
% right: x-back, y-down, z-right
Segments.lClavicle = createCS(LClavicleBase,LClavLongAxis,-ThoraxZVector,'zxy');
Segments.rClavicle = createCS(RClavicleBase,RClavLongAxis,-ThoraxZVector,'zxy');
Segments.lClavicle(5,:,:) = LClavicleEnd;
Segments.rClavicle(5,:,:) = RClavicleEnd;


%% lumbar segment

if ~IsSubjectBodyLoaded
    L5_Offset = DefaultL5Local * TrunkLength / DefaultTrunkLength;
else
    L5_Offset = DefaultL5Local * SubjectHeight;
end
LumbarJoint = getGlobalPosFromSegment(Segments.pelvis, L5_Offset);
LumbarRostralVector = ThoraxBase - LumbarJoint;
Segments.lumbar = createCS(LumbarJoint, LumbarRostralVector, PelvisVentralVector,'zyx');
Segments.lumbar(5,:,:) = ThoraxBase;


%% neck segment

if PresentThoraxTopMarker && exist('SkullBase','var')
    NeckZVector = SkullBase - ThoraxTop;
    Segments.neck = createCS(ThoraxTop, NeckZVector, ThoraxVentralVector,'zyx');
    Segments.neck(5,:,:) = SkullBase;
end


%% Define right and left arm

IsLeftElbowMarkers  = Present.left_med_epicondyle_marker  && Present.left_lat_epicondyle_marker;
IsRightElbowMarkers = Present.right_med_epicondyle_marker && Present.right_lat_epicondyle_marker;

% estimate wrist joint as the midpoint of the wrist markers
rWristJoint  = (Markers.right_radial_wrist_marker + Markers.right_ulnar_wrist_marker)/2;
lWristJoint  = (Markers.left_radial_wrist_marker  + Markers.left_ulnar_wrist_marker) /2;

% Vector as difference between wrist markers  
rWristVector = Markers.right_ulnar_wrist_marker   - Markers.right_radial_wrist_marker;
lWristVector =  Markers.left_ulnar_wrist_marker   - Markers.left_radial_wrist_marker;

% estimate position of the elbow joints and vectors
if IsLeftElbowMarkers % if median marker is present, take the midpoint of elbow markers
    lElbowJoint   = (Markers.left_med_epicondyle_marker +  Markers.left_lat_epicondyle_marker)/2;
    lElbowLatMed  =  Markers.left_med_epicondyle_marker -  Markers.left_lat_epicondyle_marker;
else
    lElbowJoint   = Markers.left_lat_epicondyle_marker;
    lElbowLatMed  = nan(3,NSamples);
end
if IsRightElbowMarkers
    rElbowJoint   = (Markers.right_med_epicondyle_marker + Markers.right_lat_epicondyle_marker)/2;
    rElbowLatMed  =  Markers.right_med_epicondyle_marker - Markers.right_lat_epicondyle_marker;
else
    rElbowJoint   = Markers.right_lat_epicondyle_marker;
    rElbowLatMed  = nan(3,NSamples);
end
% normvector to the arm plane (shoulder - elbow - wrist)
[lNormArmplane,lAng] = normvector_angle(lWristJoint, lElbowJoint, LShoulderJoint);
[rNormArmplane,rAng] = normvector_angle(rWristJoint, rElbowJoint, RShoulderJoint);
% interpolate at large angle to prevent singularities
lNormArmplane = vector_interpol(lNormArmplane, lAng, ArmExtendThreshold,Tail,false);
rNormArmplane = vector_interpol(rNormArmplane, rAng, ArmExtendThreshold,Tail,false);

% hand markers
DoUseAlternativeLHandMarker = ~Present.left_finger_marker  && Present.left_ulnar_hand_marker;
DoUseAlternativeRHandMarker = ~Present.right_finger_marker && Present.right_ulnar_hand_marker;
if DoUseAlternativeLHandMarker
    LHandMarker = Markers.left_ulnar_hand_marker;
else
    LHandMarker = Markers.left_finger_marker;
end
if DoUseAlternativeRHandMarker
    RHandMarker = Markers.right_ulnar_hand_marker;
else
    RHandMarker = Markers.right_finger_marker;
end
    
% estimate the joint axis vector of elbow: This is difficult when (almost) stretched.
% This function follows a strategy of interpolation and switching
rElbowVector = calcKneeElbowAxis(rNormArmplane,rAng,-rElbowLatMed,ArmTransit,false);
lElbowVector = calcKneeElbowAxis(lNormArmplane,lAng, lElbowLatMed,ArmTransit,false);
Segments.rHumerus = createCS(RShoulderJoint, rElbowJoint - RShoulderJoint, rElbowVector, 'zxy');
Segments.lHumerus = createCS(LShoulderJoint, lElbowJoint - LShoulderJoint, lElbowVector, 'zxy');
Segments.rRadius  = createCS(rElbowJoint,    rWristJoint - rElbowJoint,    rNormArmplane,'zxy');
Segments.lRadius  = createCS(lElbowJoint,    lWristJoint - lElbowJoint,    lNormArmplane,'zxy');
Segments.rHand    = createCS(rWristJoint,    RHandMarker - rWristJoint,    rWristVector, 'zyx');
Segments.lHand    = createCS(lWristJoint,    LHandMarker - lWristJoint,    lWristVector, 'zyx');
Segments.rHand(5,:,:) = RHandMarker;
Segments.lHand(5,:,:) = LHandMarker;


%% Define right and left Legs

% Define toes
% reconstruct second toe base marker, if needed 
% – as weighed mean of 2*inner and 1*outer marker
% - if only the inner marker is present, take that one
IsRDM2marker    = Present.right_DM2_marker;
IsRDM1_5markers = Present.right_DM1_marker && Present.right_DM5_marker;
IsLDM2marker    = Present.left_DM2_marker;
IsLDM1_5markers = Present.left_DM1_marker && Present.left_DM5_marker;
if IsRDM2marker
    RFootDistal = Markers.right_DM2_marker;
else
    if IsRDM1_5markers
        RFootDistal = (2*Markers.right_DM1_marker + Markers.right_DM5_marker)/3;
    else
        RFootDistal = Markers.right_DM1_marker;
    end
end
if IsLDM2marker
    LFootDistal  = Markers.left_DM2_marker;
else
    if IsLDM1_5markers
        LFootDistal = (2*Markers.left_DM1_marker + Markers.left_DM5_marker) /3;
    else
        LFootDistal  = Markers.left_DM1_marker;
    end
end
% tip of first toe
if Present.right_distal_hallux_marker
    RToesDistal = Markers.right_distal_hallux_marker;
end
if Present.left_distal_hallux_marker
    LToesDistal = Markers.left_distal_hallux_marker;
end


%% Define ankles

% Estimate position of the ankle joints and the lat to med vectors; 
% - if the med marker is absent, use the hip vector as a first estimate
IsLMedAnkleMarkers = Present.left_med_malleolus_marker  && Present.left_lat_malleolus_marker;
IsRMedAnkleMarkers = Present.right_med_malleolus_marker && Present.right_lat_malleolus_marker;
if IsLMedAnkleMarkers
    lAnkleLatMed =  Markers.left_med_malleolus_marker -  Markers.left_lat_malleolus_marker;
    lAnkleJoint  = (Markers.left_med_malleolus_marker +  Markers.left_lat_malleolus_marker )/2;
else
    % if no inner markers are present, use the hip vector 
    lAnkleLatMed =-PelvisYVector; 
    lAnkleJoint  = Markers.left_lat_malleolus_marker - PelvisYVector * Anklewidth/2;
end
if IsRMedAnkleMarkers
    rAnkleLatMed =  Markers.right_med_malleolus_marker - Markers.right_lat_malleolus_marker;
    rAnkleJoint  = (Markers.right_med_malleolus_marker + Markers.right_lat_malleolus_marker )/2;
else
    % if no inner markers are present, use the hip vector 
    rAnkleLatMed = PelvisYVector; 
    rAnkleJoint  = Markers.right_lat_malleolus_marker + PelvisYVector * Anklewidth/2;
end
% normalize
rAnkleLatMed = normLength(rAnkleLatMed);
lAnkleLatMed = normLength(lAnkleLatMed);

%% Define Toes
if IsRDM1_5markers
    rToesLatMed = Markers.right_DM1_marker - Markers.right_DM5_marker;
elseif IsRDM2marker && Present.right_DM5_marker
    rToesLatMed = Markers.right_DM2_marker - Markers.right_DM5_marker;
else
    rToesLatMed = rAnkleLatMed;
end
if IsLDM1_5markers
    lToesLatMed = Markers.left_DM1_marker - Markers.left_DM5_marker;
elseif IsLDM2marker && Present.left_DM5_marker
    lToesLatMed = Markers.left_DM2_marker - Markers.left_DM5_marker;
else
    lToesLatMed = lAnkleLatMed;
end


%% Define knees

% Estimate position of the knee joints and the lat to med vectors; 
% - if the med marker is absent, use the average of hip and ankle vector as a first estimate
% - if the med marker is absent, use the knee vector to estimate the joint position
IsLKneeMarkers = Present.left_femur_med_epi_marker  && Present.left_femur_lat_epi_marker;
IsRKneeMarkers = Present.right_femur_med_epi_marker && Present.right_femur_lat_epi_marker;
if IsLKneeMarkers 
    lKneeLatMed =  Markers.left_femur_med_epi_marker -  Markers.left_femur_lat_epi_marker;
    LKneeJoint  = (Markers.left_femur_med_epi_marker +  Markers.left_femur_lat_epi_marker )/2;
else
    lKneeLatMed = (lAnkleLatMed + (-PelvisYVector)) / 2;
    LKneeJoint  = Markers.left_femur_lat_epi_marker + lKneeLatMed * Kneewidth / 2;
end
if IsRKneeMarkers
    rKneeLatMed =  Markers.right_femur_med_epi_marker - Markers.right_femur_lat_epi_marker;
    RKneeJoint  = (Markers.right_femur_med_epi_marker + Markers.right_femur_lat_epi_marker )/2;
else
    rKneeLatMed = (rAnkleLatMed + PelvisYVector) / 2;
    RKneeJoint  = Markers.right_femur_lat_epi_marker + rKneeLatMed * Kneewidth / 2;
end
% min(vecnorm(Markers.left_femur_lat_epi_marker-Markers.right_femur_lat_epi_marker)) % LK: can be deleted?
% min(vecnorm(LKneeJoint - RKneeJoint)) % LK: can be deleted?

%% Norm-vectors of knee plane and ankle plane
% i.e. the plane that is spanned between toe, ankle and knee
% create vector that points to the right (i.e. -y direction)
lNormAnklePlane = normvector_angle(LFootDistal, lAnkleJoint, LKneeJoint);
rNormAnklePlane = normvector_angle(RFootDistal, rAnkleJoint, RKneeJoint);

%% if the medial ankle marker is absent, recalculate the ankle using the normvector
if ~IsLMedAnkleMarkers
    lAnkleLatMed = lNormAnklePlane;
    lAnkleJoint  = Markers.left_lat_malleolus_marker + lAnkleLatMed * Anklewidth/2;
end
if ~IsRMedAnkleMarkers
    rAnkleLatMed =-rNormAnklePlane;
    rAnkleJoint  = Markers.right_lat_malleolus_marker + rAnkleLatMed * Anklewidth/2;
end

%% if the medial knee- and ankle markers are absent, recalculate the knee also
UseAlternativeLKneeLatMed = ~IsLKneeMarkers && ~IsLMedAnkleMarkers;
UseAlternativeRKneeLatMed = ~IsRKneeMarkers && ~IsRMedAnkleMarkers;
if UseAlternativeLKneeLatMed
    lKneeLatMed = (lAnkleLatMed + (-PelvisYVector)) / 2;
    LKneeJoint  = Markers.left_femur_lat_epi_marker + lKneeLatMed * Kneewidth;
end
if UseAlternativeRKneeLatMed
    rKneeLatMed = (rAnkleLatMed + PelvisYVector) / 2;
    RKneeJoint  = Markers.right_femur_lat_epi_marker + rKneeLatMed * Kneewidth;
end

% femur long axis
rFemurVector = RKneeJoint  - RHipJoint;
lFemurVector = LKneeJoint  - LHipJoint;
% tibia long axis
rTibiaVector = rAnkleJoint - RKneeJoint;
lTibiaVector = lAnkleJoint - LKneeJoint;
% foot long axis
rFootVector  = RFootDistal - rAnkleJoint;
lFootVector  = LFootDistal - lAnkleJoint;
% toes long axis
if Present.right_distal_hallux_marker
    rToesVector = RToesDistal - RFootDistal;
else
    rToesVector = rFootVector;
end
if Present.left_distal_hallux_marker
    lToesVector = LToesDistal - LFootDistal;
else
    lToesVector = lFootVector;
end

% Femur coordinate system; use knee as reference
Segments.rFemur = createCS(RHipJoint,   rFemurVector, -rKneeLatMed, 'zxy');
Segments.lFemur = createCS(LHipJoint,   lFemurVector,  lKneeLatMed, 'zxy');
% Tibia coordinate system; use knee as reference
Segments.rTibia = createCS(RKneeJoint,  rTibiaVector, -rKneeLatMed, 'zxy');
Segments.lTibia = createCS(LKneeJoint,  lTibiaVector,  lKneeLatMed, 'zxy');
% Foot Coordinate system:  use ankle as reference
Segments.rFoot  = createCS(rAnkleJoint, rFootVector,  -rAnkleLatMed,'zxy');
Segments.lFoot  = createCS(lAnkleJoint, lFootVector,   lAnkleLatMed,'zxy');
% Toes
Segments.rToes  = createCS(RFootDistal, rToesVector,  -rToesLatMed, 'zxy');
Segments.lToes  = createCS(LFootDistal, lToesVector,   lToesLatMed, 'zxy');

% LfemurLen = vecnorm(LHipJoint-LKneeJoint);
% RfemurLen = vecnorm(RHipJoint-RKneeJoint);
% figure;hold on;plot(LfemurLen);plot(RfemurLen)

%% define endpoints of distal foot segments
Segments.rFoot(5,:,:) = RFootDistal;
Segments.lFoot(5,:,:) = LFootDistal;
if Present.right_distal_hallux_marker
    Segments.rToes(5,:,:) = RToesDistal;
end
if Present.left_distal_hallux_marker
    Segments.lToes(5,:,:) = LToesDistal;
end
end


%% ========================== LOCAL FUNCTIONS =======================================
%% ==================================================================================


function OutVect = calcKneeElbowAxis(Vect,Angle,AltVect,TransitRng,isPlot)
% Handle elbow and knee joint vector when in near extension
% Use the AltVect with a smooth transition phase for the angular range TransitRng
% Vect         vector of jont axis as time series [3,nsamples]
% Angle        extension angle (180 deg is extended) [1,nsamples]
% AltVect      alternative vector (deg) [3,nsamples]
% TransitRng   transition range (deg) [min max]
% isPlot       (optional)
% 
% Author: Marc de Lussanet
% version 220915 : better checking for missing vectors


% Check if any or both Vect, AltVect are missing
if all(isnan(Vect),'all') && (isempty(AltVect) || all(isnan(AltVect),'all'))
    OutVect = nan(size(Vect));
    return;
elseif all(isnan(Vect),'all')
    OutVect = normLength(AltVect);
    return;
elseif (isempty(AltVect) || all(isnan(AltVect),'all'))
    OutVect = normLength(Vect);
    return;
else % both vectors are present: interpolate across the range
    Vect    = normLength(Vect);
    AltVect = normLength(AltVect);
    
    % theoretically, there can still be gaps at the beginning or end of the time series, because
    % gapfilling does not extrapolate. In that case, overwrite the nans with the other vector
    if any(isnan(AltVect),'all')
        AltVect(isnan(AltVect)) = Vect(isnan(AltVect));
    end
    if any(isnan(Vect),'all')
        Vect(isnan(Vect)) = AltVect(isnan(Vect));
    end  

    % Switch the regions of near extension with the AltVect across the TransitRange
    Interp = Vect+(Angle-TransitRng(1))/(TransitRng(2)-TransitRng(1)).*(AltVect-Vect);
    Small  = Angle<TransitRng(1);
    Large  = Angle>TransitRng(2);
    OutVect          = Interp;
    OutVect(:,Small) = Vect(:,Small);
    OutVect(:,Large) = AltVect(:,Large);
    OutVect = normLength(OutVect);
end


% plotting if desired
if isPlot
    figure; title('interpolate norm-vector in vicinity of full extension')
    subplot(2,2,1); hold on; plot(OutVect(1,:),'linewidth',1); plot(Vect(1,:));plot(AltVect(1,:));
    subplot(2,2,2); hold on; plot(OutVect(2,:),'linewidth',1); plot(Vect(2,:));plot(AltVect(2,:));
    subplot(2,2,3); hold on; plot(OutVect(3,:),'linewidth',1); plot(Vect(3,:));plot(AltVect(3,:));
    subplot(2,2,4); hold on; plot(Angle); yline(TransitRng(1)); yline(TransitRng(2)); title('angle');
end
end

%% ========================================================================

function [OutVect,Gaps] = vector_interpol(Vector,Angle,ExtThreshold,Tail,DoPlot,MaxGap)
%% Interpolate Vector when angle>threshold.
% The idea: if a joint vector (e.g. elbow, knee) is estimated from the normal vector to the flanking 
% segments, the vector becomes really unreliable near extension. Thus, interpolate phases near extension 
% 
% SYNTAX
%   [OutVect,Gaps] = vector_interpol(Vector,Angle,ExtThreshold,Tail,isPlot,MaxGap)
%   OutVect = vector_interpol(Vector,Angle,ExtThreshold,Tail)
%
% INPUT
%     Vector       - (Nsamp x 3 double) time series of 3D vector
%     Angle        - (Nsamp double) time series of joint angle
%     ExtThreshold - (double) angle [deg] threshold
%     Tail         - (int) number of samples in the tail for interpolation
%     isPlot       - (logical) optional flag; default false
%     MaxGap       - (int) maximal length of gaps to be filled [samples]
%
% OUTPUT
%     OutVect - (Nsamp x 3 double) time series of 3D vector
%     Gaps    - (N x 2 double) start and end sample of the filled periods
%
% See also: calcMainCS
% 
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Last revision: 220407 added header and padding for extrapolation; normalize the length
% Last revision: 230217 (MdL) corrected the dim of the  normlength calculation

% handle optional input arguments
narginchk(4,6);
if nargin<5 || isempty(DoPlot), DoPlot = false; end
if nargin<6 || isempty(MaxGap), MaxGap = 0;     end

if all(isnan(Vector),'all') || all(isnan(Angle))
    OutVect = Vector;
    Gaps    = [1 length(Vector)];
    return;
end
Deg     = 2;
Out1    = Vector';
Out0    = Vector;
Out0(:,Angle>ExtThreshold)= nan;
Out1(:,1)        = polyfillgaps(Out0(1,:)',Tail,Deg);
Out1(:,2)        = polyfillgaps(Out0(2,:)',Tail,Deg);
[Out1(:,3),Gaps] = polyfillgaps(Out0(3,:)',Tail,Deg);

if ~isempty(Gaps)
    % omit too long periods of extension
    for i=1:size(Gaps,1)
        if MaxGap && Gaps(i,2) > MaxGap  % if MaxGap==0, fill everything
            Out1(Gaps(i,1) : Gaps(i,1)+Gaps(i,2)-1, :) = nan;
        elseif i==1 && Gaps(1,1)==1  % gap at beginning
        elseif i==size(Gaps,1) && Gaps(end,1)+Gaps(end,2)-1 == length(Vector) % gap at end
        else
            Gaps(i,:) = nan; % mark this gap as filled
        end
    end
    % omit extrapolations
    if Gaps(1,1)==1
        Out1(1:Gaps(1,2), :)=nan;
    end
    if Gaps(end,1)+Gaps(end,2)-1 == length(Vector)
        Out1(Gaps(end,1):end, :) = nan;
    end
    % remove the nan values
    Gaps(isnan(Gaps(:,1)), :) = [];
end
% normalize length
Out1 = normLength(Out1,2); % scaled to unit length

% the interpolation process is now ready
OutVect = Out1';

% Extrapolation does not work well, so assume constant axis
if ~isempty(Gaps)
    if Gaps(1,1)==1  % beginning of trajectory
        OutVect(1:Gaps(1,1), :)  = OutVect(Gaps(1,  1)+1, :);
    end
    if Gaps(end,1)+Gaps(end,2)-1 == length(OutVect) && Gaps(end,1)>1  % end of trajectory
        OutVect(Gaps(end,1):end, :) = OutVect(Gaps(end,1)-1, :);
    end
end

% figure if desired
if DoPlot %&& any(any(isnan(Out0)))
    Gap  = max(Gaps(:,2)); if isempty(Gap), Gap=0; end
    Max  = max(max(abs(Out1)));
    if MaxGap && Gap > MaxGap || Max > 1
        warning('Long periods of extended limb (%d samp, max %d); or too large %f.',Gap,MaxGap,Max);
    end
    figure; title('interpolate norm-vector in vicinity of full extension')
    subplot(2,1,1); hold on; plot(Vector','color','r');plot(Out1); title('vector\_interpol: 3D vector');
    subplot(2,1,2); hold on; plot(Angle); yline(ExtThreshold); title('angle and threshold');
end
end

%% ========================================================================



