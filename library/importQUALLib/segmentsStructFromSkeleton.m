function Segments = segmentsStructFromSkeleton(GlobalMN)
%% Calculate segments structure for markerless data
% The Segments structure is necessary for placing force plate data on the force contact points.
% 
% SYNTAX
%   Segments = segmentsStructFromSkeleton(GlobalMN)
%
% INPUT
%     GlobalRotMN  - (struct MNData) spatial orientation of all segments in global coordinates
% 
% OUTPUT
%     Segments     - (struct) structure of the Segments
%
% EXAMPLE 
%   Segments = segmentsStructFromSkeleton(GlobalMN)
%
% Local functions:
%
% See also: importQUAL, extractForcePlates
% 
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% (c) 2023 by Predimo GmbH
% version 230204 (MdL) outsourced from placeVirtualMarkersOnJoints

% Init
DefBody  = loadBody('def_body');
Segments = [];
if isempty(GlobalMN)
    return;
end

%% Get the global orientation of desired segments

% loop the segments as defined in the segments table
JointNames = {GlobalMN.signals.joints.name};
for SegmentIdx = 1:length(DefBody.segments.id)

    % The label of the segments that is written in the segmentLabel column of the segments table
    SegmentLabel = DefBody.segments.segmentLabel{SegmentIdx};
    % If a segmentLabel is provided, save the global position and orientation to the segments structure
    if ~isempty(SegmentLabel)

        % The rownumber of the segment (if present) in the conversion table
        ConversionIdx = matches(DefBody.conversion.segmentName,DefBody.segments.name{SegmentIdx});
        IsSegentDefined = any(ConversionIdx);
        if IsSegentDefined

            % Name of the joint in the GlobalRotMN structure
            JointName = DefBody.conversion.name{ConversionIdx};
            % item number of the joint in the GlobalRotMN structure
            GlobalMN_No = matches(JointNames,JointName);
            if any(GlobalMN_No)
                % write global position and rotation to the segments structure as CS (coordinate system)
                GlobalR = GlobalMN.signals.joints(GlobalMN_No).dynamic.globalR.data;
                GlobalR = permute(GlobalR,[2 3 1]); % NSamples x3x3 -> 3x3xNsamples
                Segments.(SegmentLabel)(1:3,:,:) = GlobalR;
                GlobalPos = GlobalMN.signals.joints(GlobalMN_No).dynamic.pos.data;
                GlobalPos = GlobalPos';
                Segments.(SegmentLabel)(4,:,:) = GlobalPos;
            end
        end
    end
end
