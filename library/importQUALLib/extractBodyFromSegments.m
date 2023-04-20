function [SubjectBody,StablePeriod] = extractBodyFromSegments(Segments, DefBody, Markers, Present, Info)
%% Estimate body segment lengths from marker-defined segments to fill the segments table
%% fill the segments sheet for the user_body table
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
% Local functions: 
%
% See also: extractBodyFromQUALSkel, extractBodyFromTheia, extractHorseFromSegments
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230215 (MdL) Correction in clavicle length; FIX in body height estimation; workaround fro higher floor level
% in drop jump
% version 230216 (MdL) FIXED right humerus length calculation; 
% version 230217 (MdL) commented out adjustment for the RBM nextpos (not yet ready)

%% Params
DoPlotVar = false; % do plot on error in calcSegmentSize
DoPlotStable = false; % plot the selected period of minimal movement (i.e. the "stable pose")

%% Init
X = 1; Y = 2; Z = 3;
% NextPos       = nan(size(DefBody.segments.nextPos));
SizeDefault   = DefBody.segments.size; % default segment sizes
NextPosDefault= DefBody.segments.nextPos;
MeanLengths   = nan(size(SizeDefault));
SegmentNames  = DefBody.segments.name;

%% estimate the stable pose period
StablePeriod  = getStablePeriod(Markers,Present,Info,DoPlotStable);


%% Estimate body height
DefaultHeight = DefBody.general.height;
SubjectHeight = estimateBodyHeightFromMarkers(Segments, DefaultHeight, StablePeriod);


%% Get the segment lengths from Segment CSs

%% Pelvis (at mid-asi) and RBM (at mid-hip) have the same orientation and size
% the size is the "default asi distance" but scaled to the actual hip distance
PlvNo = startsWith(SegmentNames,'pelvis_');
RBMNo = startsWith(SegmentNames,'RBM_');
LFemurID = 27;
RFemurID = 26;
LHipNextNo = DefBody.segments.next(PlvNo,:) == LFemurID; % which entry of nextpos gives the left hip
RHipNextNo = DefBody.segments.next(PlvNo,:) == RFemurID; % which entry of nextpos gives the right hip
% pelvis default width scaled to the body height
PelvisWidthDefault = SizeDefault(PlvNo) * SubjectHeight; % [m]
% The position of the hip joints is given by nextPos, expressed in the pelvis width
NormalizedLHipPosDefault = squeeze(NextPosDefault(PlvNo, LHipNextNo, :));
NormalizedRHipPosDefault = squeeze(NextPosDefault(PlvNo, RHipNextNo, :));
NormalizedRBMPosDefault  = (NormalizedLHipPosDefault + NormalizedRHipPosDefault) / 2;
NormalizedHipWidthDefault = abs(NormalizedLHipPosDefault(Y) - NormalizedRHipPosDefault(Y));
HipDistanceDefault = PelvisWidthDefault * NormalizedHipWidthDefault;
HipDistance = mean(vecnorm(squeeze(Segments.lFemur(4,:,:) - Segments.rFemur(4,:,:))),'omitnan');

% Pelvis scaling factor in X
% PelvisScalingFactor_X = HipDistance / HipDistanceDefault;
% Pelvis scaling factor in Y
PelvisScalingFactor_Y = HipDistance / HipDistanceDefault;
% Pelvis scaling factor in Z
% PelvisScalingFactor_Z = HipDistance / HipDistanceDefault;

PelvisScalingFactor = PelvisScalingFactor_Y;
PelvisWidth = PelvisWidthDefault * PelvisScalingFactor;
% set the sizes
MeanLengths(PlvNo) = PelvisWidth / SubjectHeight;
MeanLengths(RBMNo) = PelvisWidth / SubjectHeight;

% shift the X and Z positons of the RBM, such that it lands at mid-hips
RBMPosDefault = PelvisWidthDefault * NormalizedRBMPosDefault;
LHipLocal = getLocalPosOnSegment(Segments.pelvis,squeeze(Segments.lFemur(4,:,:)));
RHipLocal = getLocalPosOnSegment(Segments.pelvis,squeeze(Segments.rFemur(4,:,:)));
RBMPosLocal = (mean(LHipLocal + RHipLocal, 2, 'omitnan')) / 2;
RBMnextpos = - RBMPosLocal / (PelvisWidth * SubjectHeight);
RBMnextposDefault = squeeze(NextPosDefault(RBMNo, 1, :));
NextPosRBM = RBMnextpos - RBMnextposDefault;
% NextPosRBM = RBMnextposDefault - RBMnextpos;

% report
fprintf('Estimated shift of  %22s: %6.1f%6.1f%6.1f cm %6.2f%6.2f%6.2f\n','hips', ...
    RBMPosLocal(X)*100,RBMPosLocal(Y)*100,RBMPosLocal(Z)*100, NextPosRBM(1),NextPosRBM(2),NextPosRBM(3));
fprintf('Estimated distance between %15s: %4.1f\n','hips',HipDistance*100);

%% Trunk
% Lumbar Back (this is set using the "trans" field in calcJointStream)

% Ribcage
ThoraxNo = startsWith(SegmentNames,'ribcage_');
MeanLengths(ThoraxNo) = calcSegmentSize(ThoraxNo,Segments.thorax, [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% Neck (this is set using the "trans" field in calcJointStream)

% Head
Hd = startsWith(SegmentNames,'skull_');
MeanLengths(Hd) = calcSegmentSize(Hd,Segments.skull,  [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% Arms: clavicle, humerus, radius
% Scapula
% RS = startsWith(SegmentNames,'right_scapula_');
% LS = startsWith(SegmentNames,'left_scapula_');
% if Present.right_shoulder_back_marker
%     MeanLengths(RS) = calcSegmentSize(RS,Segments.rHumerus, Markers.right_shoulder_back_marker, DefBody,SubjectHeight,StablePeriod,DoPlotVar);
% end
% if Present.left_shoulder_back_marker
%     MeanLengths(LS) = calcSegmentSize(LS,Segments.lHumerus, Markers.left_shoulder_back_marker, DefBody,SubjectHeight,StablePeriod,DoPlotVar);
% end

% Clavicle
RClavNo = startsWith(SegmentNames,'right_clavicle_');
LClavNo = startsWith(SegmentNames,'left_clavicle_');
MeanLengths(RClavNo) = calcSegmentSize(RClavNo,Segments.rClavicle,[],DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LClavNo) = calcSegmentSize(LClavNo,Segments.lClavicle,[],DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% Humerus
RH = startsWith(SegmentNames,'right_humerus_');
LH = startsWith(SegmentNames,'left_humerus_');
MeanLengths(RH) = calcSegmentSize(RH,Segments.rHumerus, Segments.rRadius, DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LH) = calcSegmentSize(LH,Segments.lHumerus, Segments.lRadius, DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% Radius
RR = startsWith(SegmentNames,'right_radius_');
LR = startsWith(SegmentNames,'left_radius_');
MeanLengths(RR) = calcSegmentSize(RR,Segments.rRadius,  Segments.rHand,   DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LR) = calcSegmentSize(LR,Segments.lRadius,  Segments.lHand,   DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% Hand
RM = startsWith(SegmentNames,'right_hand_');
LM = startsWith(SegmentNames,'left_hand_');
MeanLengths(RM) = calcSegmentSize(RM,Segments.rHand, [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LM) = calcSegmentSize(LM,Segments.lHand, [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);

%% Legs
% femur
RF = startsWith(SegmentNames,'right_femur_');
LF = startsWith(SegmentNames,'left_femur_');
MeanLengths(RF) = calcSegmentSize(RF,Segments.rFemur,  Segments.rTibia, DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LF) = calcSegmentSize(LF,Segments.lFemur,  Segments.lTibia, DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% tibia
RT = startsWith(SegmentNames,'right_tibia_');
LT = startsWith(SegmentNames,'left_tibia_');
MeanLengths(RT) = calcSegmentSize(RT,Segments.rTibia,  Segments.rFoot,  DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LT) = calcSegmentSize(LT,Segments.lTibia,  Segments.lFoot,  DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% foot
RP = startsWith(SegmentNames,'right_foot_');
LP = startsWith(SegmentNames,'left_foot_');
MeanLengths(RP) = calcSegmentSize(RP,Segments.rFoot, Segments.rToes, DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LP) = calcSegmentSize(LP,Segments.lFoot, Segments.lToes, DefBody,SubjectHeight,StablePeriod,DoPlotVar);

% toes
RZ = startsWith(SegmentNames,'right_toes_');
LZ = startsWith(SegmentNames,'left_toes_');
MeanLengths(RZ) = calcSegmentSize(RZ,Segments.rToes, [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);
MeanLengths(LZ) = calcSegmentSize(LZ,Segments.lToes, [], DefBody,SubjectHeight,StablePeriod,DoPlotVar);


%% Get the factor by which the segment was scaled with respect to the default
ScaleToDefault = MeanLengths ./ DefBody.segments.size;
% NextPosAll     = nan(size(DefBody.segments.nextPos));
% NextPosAll(RBMNo,1,:) = NextPosRBM;

%% Write estimated segment sizes
SubjectBody.segments.id      = DefBody.segments.id;
SubjectBody.segments.name    = DefBody.segments.name;
SubjectBody.segments.scale   = ScaleToDefault;
SubjectBody.segments.size    = MeanLengths;
% SubjectBody.segments.nextPos = NextPosAll;
SubjectBody.general.height   = SubjectHeight;
end


% ==================================================================================================

% function NextPos_Z_User = scaleNextPos2UserVertColumn(DefBody, VertebralNos, SegmentLength, BodyHeight)
% %% correct the radius of gyration of neck and lumbar vertebrae
% % Note to radius of gyration (see example of 221105, MdL): 
% % - Moment of inertia with respect to centre of mass: I = m * r_(x,y,z)^2
% % - Radius of gyration g = √(I/m) with unit [m]
% % - Since the radius r is normalized to the segment length (i.e. r_z), it follows :
% %     Suppose ellipsoid with r_(x0,y0,z0) is scaled by factor F 
% %     - g_x, g_y scale with 1/√F
% %     - g_z scales with 1/F
% %% Since this does not work out correctly, we scale the nextpos rather than the gyr and segment size
% 
% %% segment lengths are normalized to body size
% % the total length of the user segment normalized to user body height
% NormTotalSegmentLengthUser = SegmentLength / BodyHeight;
% % the total length of the segment according to def_body
% NormTotalSegmentLengthDefault = sum(DefBody.segments.size(VertebralNos) .* DefBody.segments.nextPos(VertebralNos,1,3));
% % the difference in size that is to be compensated by scaling the nextpos
% NextPosCorrection = NormTotalSegmentLengthUser - NormTotalSegmentLengthDefault;
% 
% %% nextpos are normalized to segment length
% % -> get the intervertebral space (i.e. if nextpos_z = 1.4, this is 0.4 time segment length)
% % the relative lengths of the intervertebral spaces according to nextpos
% IntervertebralSpacesDefault = DefBody.segments.nextPos(VertebralNos,1,3) - 1; % this is a factor so minus 1 is the difference
% % the sizes of the intervertebral spaces scaled to body size
% NormIntervertebralSpaceDefault = sum(DefBody.segments.size(VertebralNos) .* IntervertebralSpacesDefault);
% % the scaling factor
% ScaleFactorNextPos = (NextPosCorrection + NormIntervertebralSpaceDefault) / NormIntervertebralSpaceDefault;
% % scale the intervertebral space and add one to obtain a scaling factor again
% NextPos_Z_User = (IntervertebralSpacesDefault * ScaleFactorNextPos) + 1;
% % Totallength_check = sum(DefBody.segments.size(VertebralNos) .* NextPos_Z_User)
% end

% ==================================================================================================
