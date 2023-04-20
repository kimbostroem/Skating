function UserBody = extractBodyFromTheia(SkelData, SkelParam)
%% Estimate body segment sizes from Theia
% UserBody.general.bodyHeight
% Theia documentation:
%   https://www.theiamarkerless.ca/docs/model.html#model-description 
% Position of hip joints in Pelvis: 
%   Harrington et al. (2007). Prediction of the hip joint centre in adults, children, and patients 
%   based on magnetic resonance imaging. J. Biomech. doi: 10.1016/j.jbiomech.2006.02.003.
%
% SYNTAX
% UserBody = extractBodyFromTheia(SkelData, SkelParam)
%
% INPUT
%     SkelData        (struct) Skeleton segment data
%     SkelParam       (struct) Skeleton parameters from loadc3d
%
% OUTPUT
%    UserBody         (struct) Structure containing body parameters
%
% Local functions: calcAnkleHeight, fprintfATableOfEstimatedLengths
%
% See also: calcMainCS, extractBodyFromQUALSkel, importTheiaSkel, calcSegmentLength
%
% (c) 2022 by Predimo GmbH, http://www.predimo.com 
% Author: Marc de Lussanet, Kim Bostroem 
% version 230214 (MdL) fixed unknown index "PelvisNo"

%% Preparations

% get the conversion table
SrcType = 'theia';
DefBody = loadBody('def_body');
% copy the Theia fields to the generic fields (for calcSegmentLength)
DefBody.conversion.srcName        = DefBody.conversion.([SrcType, '_name']);
DefBody.conversion.srcOrientation = DefBody.conversion.(SrcType);

% Intermediate parameters
Y = 2;
NextPosDefault  = DefBody.segments.nextPos;
SizeDefault     = DefBody.segments.size; % default segment sizes
SkelSegData     = SkelData.signals.joints;
SkelSegNames    = {SkelSegData.name};
SkelParamGroups = {SkelParam.group};
SkelParams      = {SkelParam.parameter};
TheiaParamNos   = matches(SkelParamGroups, 'THEIA3D'); % No for row number in table
ThoraxNo        = matches(SkelParams, 'THORAX_LENGTH') & TheiaParamNos;

% Initialize
UserBody     = struct;
SegmentSizes = nan(size(DefBody.segments.size)); % init segment size array with NaN
NextPos      = nan(size(DefBody.segments.nextPos)); % init nextpos with NaN


%% Calc body measures

AnkleHeight       = calcAnkleHeight(SkelSegNames, SkelSegData);
[FemurLength,D1]  = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [27 31; 26 30]); % hip knee
[TibiaLength,D2]  = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [31 33; 30 32]); % knee ankle
LegLength         = TibiaLength + FemurLength;
HipWidth          = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [27 26]); % hips
% Torso between neck and hips: use the max rather than the mean, in case there is bending during the measurement
[TorsoLength,D3]  = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [2 20], [],'max'); % thorax pelvis
% Cervical vertebral column between neck and head: use the max (see above @ TorsoLength)
[CervicalLength,D4]=calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [ 1  2], [],'max'); % neck head
ClavicleLength    = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [13 15; 12 14]); % clav upparm
% Theia does not export the head length, so assume the default
DefaultHeadLength = DefBody.segments.size(1) * DefBody.general.height;

% Calculate body height
UseDefaultHeight = D1 & D2 & D3 & D4;
if UseDefaultHeight
    % if none of the segments can be estimated, use the default height.
    BodyHeight = DefBody.general.height;
    fprintf('Cannot estimate body height from Theia data: assuming default (%.2fm)\n', BodyHeight);
else
    % if at least one segment length can be calculated
    BodyHeight = AnkleHeight + LegLength + TorsoLength + CervicalLength + DefaultHeadLength;
    BodyHeight = round(BodyHeight, 2); % round to cm
    fprintf('Estimated height: %.2f m (estimated from Theia Skeleton)\n', BodyHeight);
end

% Compute the scaling factor for the pelvis
PelvisID = 20;
PelvisNo = DefBody.segments.id == PelvisID;
LFemurID = 27;
RFemurID = 26;
LHipNextNo = DefBody.segments.next(PelvisNo,:) == LFemurID; % which entry of nextpos gives the left hip
RHipNextNo = DefBody.segments.next(PelvisNo,:) == RFemurID; % which entry of nextpos gives the right hip
% Pelvis width is the Asi-distance in m for the default body
PelvisWidthDefault = SizeDefault(PelvisNo) * BodyHeight; % [m]
% The position of the hip joints is given by nextPos, expressed in the pelvis width
NormalizedHipWidthDefault = abs(NextPosDefault(PelvisNo, LHipNextNo, Y) - NextPosDefault(PelvisNo, RHipNextNo, Y));
HipWidthDefault = PelvisWidthDefault * NormalizedHipWidthDefault;
PelvisScalingFactor = HipWidth / HipWidthDefault;

% Compute the length of the lumbar vertebra
ThoraxLength         = str2double(SkelParam(ThoraxNo).data{1}); % given as additional parameter in data

%% Calc segment sizes

for SegmentID = 1:length(DefBody.segments.id)
    SegmentNo = DefBody.segments.id == SegmentID;
    switch SegmentID

        case 36 % "RBM" is welded to the pelvis. 
            % Its origin is at the mid-asi and the orientation is aligned to the pelvis
            
            SegmentSizes(SegmentNo) = SizeDefault(SegmentNo) * PelvisScalingFactor;

        case 20 % Pelvis segment -> hip joints and to lumbar joint

            SegmentSizes(SegmentNo) = SizeDefault(SegmentNo) * PelvisScalingFactor;

        case {21,22,23,24,25} % lumbar vertebra
            % do not estimate: these are estimated later through the translation

        case 2 % thorax segment
            NormThoraxLength = ThoraxLength / BodyHeight;
            SegmentSizes(SegmentNo) = NormThoraxLength;

        case {3,4,5,6,7,8,9} % cervical vertebra
            % do not estimate: these are estimated later through the translation

        case {12,13} % clavicle segment
            % the Theia clavicle segment originates on the neck joint and inserts on the humerus
            % base. 
            % Currently, the clavicle is rotated backward by 20 deg in Myonardo.slx; the according
            % shift of the clavicle base is defined by nextPos
            SegmentSizes(SegmentNo) = ClavicleLength / BodyHeight;
            
        case 18 % right hand segment (with fingers)
            % the theia representation is without fingers, so multiply by two to get a plausible length 
            RHandNo = matches(SkelParams,'RHAND_LENGTH') & TheiaParamNos;
            SegmentSizes(SegmentNo) = 2 * str2double(SkelParam(RHandNo).data{1}) / BodyHeight;

        case 19 % left hand segment (with fingers)
            % the theia representation is without fingers, so multiply by two to get a plausible length 
            LHandNo = matches(SkelParams,'LHAND_LENGTH') & TheiaParamNos;
            SegmentSizes(SegmentNo) = 2 * str2double(SkelParam(LHandNo).data{1}) / BodyHeight;

        case 34 % right toes segment
            RToeNo = matches(SkelParams,'RTOE_LENGTH') & TheiaParamNos;
            SegmentSizes(SegmentNo) = str2double(SkelParam(RToeNo).data{1}) / BodyHeight;

        case 35 % left toes segment
            LToeNo = matches(SkelParams,'LTOE_LENGTH') & TheiaParamNos;
            SegmentSizes(SegmentNo) = str2double(SkelParam(LToeNo).data{1}) / BodyHeight;

        otherwise % all other segments
            % 1     head : midpoint between the ears (not defined in theia)
            % 10,11 scapula  -> otherwise
            % 14,15 humerus  -> otherwise
            % 16,17 forearm  -> otherwise
            % 26,27 femur  -> otherwise
            % 30,31 tibia  -> otherwise
            % 32,33 foot   -> otherwise
            
            Size = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, SegmentID);
            % size relative to nextPos and body height
            ScaleFactor = 1 / (abs(NextPosDefault(SegmentNo, 1, 3)) * BodyHeight);
            SegmentSizes(SegmentNo) = Size * ScaleFactor;
    end
    % if the size is not defined, it's zero and should be set to NaN
    if SegmentSizes(SegmentNo) == 0
        SegmentSizes(SegmentNo) = NaN;
    end
end


%% Store parameters

% indices of defined segment sizes
IsDefined = ~isnan(SegmentSizes) | ~isnan(NextPos(:,1,1));
% Segment scaling factors
Scales = SegmentSizes ./ DefBody.segments.size;
% fill the user body
UserBody.general.height = BodyHeight;
UserBody.segments.id = DefBody.segments.id(IsDefined);
UserBody.segments.name = DefBody.segments.name(IsDefined);
UserBody.segments.size = SegmentSizes(IsDefined);
UserBody.segments.scale = Scales(IsDefined);

% fprintfATableOfEstimatedLengths(UserBody);
end

%% =================================================================================================

function AnkleHeight = calcAnkleHeight(SkelSegNames, SkelSegData)
%% calculate ankle height
% version 221222 (MdL) first version

MinPlausibleHeight = 0.04; % m
PlausibleHeight = 0.06; % m
LAnkleNo = matches(SkelSegNames, 'l_foot_4X4');
RAnkleNo = matches(SkelSegNames, 'r_foot_4X4');
LAnkleHeight = min(SkelSegData(LAnkleNo).dynamic.pos.data(:,3));
RAnkleHeight = min(SkelSegData(RAnkleNo).dynamic.pos.data(:,3));
AnkleHeight = min([LAnkleHeight RAnkleHeight]);
if AnkleHeight < MinPlausibleHeight
    AnkleHeight = PlausibleHeight;
end
end

%% =================================================================================================

% function LumbarVertebraLength = calcLumbarLength(ThoraxLength,MidhipToL5Length,TorsoLength,DefBody)
% %% calculate length of the lumbar spine: Torso - Thorax - Sacrum
% % Plausibility check: if the lumbar length is too short, the default length is returned
% % version 230125 (MdL) comments and explicit min plausible length
% 
% MinPlausibleLength = 0.05; % Lumbar size should be more than 5 cm
% 
% LumbarVertebraLength = TorsoLength - ThoraxLength - MidhipToL5Length;
% 
% if LumbarVertebraLength < MinPlausibleLength
%     warning('extractBodyFromTheia: Length of Lumbar column (%.3fm) is not plausible.\n\tPlease use a better body file with an upright pose', LumbarVertebraLength);
%     TorsoIds = [20 25 24 23 22 21 2]; % pelvis, L5-L1, sacrum
%     SegmentIds = [DefBody.segments.id];
%     TorsoNos = sum(SegmentIds == TorsoIds, 2,'native');
%     LumbarVertebraLength = sum(DefBody.segments.size(TorsoNos)) * DefBody.general.height;
% end
% end

%% =================================================================================================

function fprintfATableOfEstimatedLengths(UserBody) %#ok<DEFNU> 
BodyHeight = UserBody.general.height;
fprintf('Estimated segment lengths (cm) and next positions (np) for body height %.2fm:\n\n',BodyHeight);
fprintf('%5s %25s %6s %5s |%6s %6s %6s |%6s %6s %6s |%6s %6s %6s\n', ...
    'id','name','len_cm','scale','np1-x','np1-y','np1-z','np2-x','np2-y','np2-z','np3-x','np3-y','np3-z')
for No = 1:length(UserBody.segments.id)
    if ~isempty(No)
        fprintf('%5d %25s %6.1f %5.2f |%6.2f %6.2f %6.2f |%6.2f %6.2f %6.2f |%6.2f %6.2f %6.2f\n', ...
            UserBody.segments.id(No), UserBody.segments.name{No}, ...
            UserBody.segments.size(No)*BodyHeight*100, UserBody.segments.scale(No), ...
            UserBody.segments.nextPos(No,1,1), UserBody.segments.nextPos(No,1,2), UserBody.segments.nextPos(No,1,3), ...
            UserBody.segments.nextPos(No,2,1), UserBody.segments.nextPos(No,2,2), UserBody.segments.nextPos(No,2,3), ...
            UserBody.segments.nextPos(No,3,1), UserBody.segments.nextPos(No,3,2), UserBody.segments.nextPos(No,3,3) ...
            );
    end
end
end
