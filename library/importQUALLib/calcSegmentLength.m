function [SegmentLength,IsDefault] = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, SegmentIds, Samples, Measure)
%% Calculate segment length on the basis of the listed SegmentIds
% - SegmentIds represent current segment and the follower one.
% - If SegmentIds lists more than one segment (e.g. left and right shin), the average is returned
%   (e.g. for estimating the body length).
% - The range of samples over which the length is calculated can be defined (e.g. the first sample
%   in skeleton data)
% - If no length can be calculated, the default length according to DefBody is returned
% 
% SYNTAX
% [SegmentLength] = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, SegmentIds, Samples, Measure)
%
% INPUT
%     SkelSegNames (cell array) names of the segments in the input data
%     SkelSegData  (struct) Input data
%     DefBody      (struct) Note, that the DefBody.conversion.srcName should be defined
%     SegmentIds   (double) 2 x NSegments matrix of Id of the segments (see explanation above and
%                  example below)
%     Samples      (double) [default all). Samples across which the length is calculated
%     Measure      (char) [default 'mean']. Other: 'max', 'min', 'median'
%
% OUTPUT
%     SegmentLength (double) length (m)
%     IsDefault     (logical) true if the default value is returned (no value could be calculated)
%
% Examples:
% average lower leg length: Ids give l,r shin and foot segments. From first sample
%     FemurLength = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [27 31; 26 30], 1);
% maximal cervical length : Ids give neck and head segments. From ALL samples
%     NeckLength = calcSegmentLength(SkelSegNames, SkelSegData, DefBody, [1 2], [],'max');
% 
% See also: extractBodyFromTheia
%
% (c) 2022 by Predimo GmbH, http://www.predimo.com 
% Author: Marc de Lussanet
% Version 230106 (MdL) commenting and improved header for review

% Init
NSegments  = size(SegmentIds,1);
narginchk(4,6);
NSamples = length(SkelSegData(1).dynamic.pos.data(:,1));
if nargin<5 || isempty(Samples), Samples = 1:NSamples; end
if nargin<6 || isempty(Measure), Measure = 'mean'; end
IsDefault = false;

% if the distal id is not provided, retrieve it from the segments table
if size(SegmentIds,2)==1
    % add a row for the distal position
    SegmentIds(:,2) = nan(length(SegmentIds),1);
    for iSeg = 1:NSegments
        SegmentNo_segments = DefBody.segments.id == SegmentIds(iSeg);
        NextSegId = DefBody.segments.next(SegmentNo_segments,1);
        SegmentIds(iSeg,2) = NextSegId;
    end
end

% loop the segments
SegLengths = nan(1,NSegments);
for iSeg = 1:NSegments
    % get the row numbers in the conversion table
    SegmentNo_Conv = sum(DefBody.conversion.segmentID == SegmentIds(iSeg,:), 2,'native');
    SegmentNames = DefBody.conversion.srcName(SegmentNo_Conv);
    SkelSegNos = find(matches(SkelSegNames, SegmentNames));
    % only fill a value if the segment is in the data
    IsSegmentPresent = ~isempty(SkelSegNos) && length(SkelSegNos)==2;
    if IsSegmentPresent
        % row numbers as output arguments
        % get the position of base and and distal side of the segment for the selected range of samples
        BasePos = SkelSegData(SkelSegNos(1)).dynamic.pos.data(Samples,:);
        DistPos = SkelSegData(SkelSegNos(2)).dynamic.pos.data(Samples,:);
        % calculate segment length over time
        Length_t = vecnorm(BasePos-DistPos, 2,2); % ...,2,2) -> euclidean length on 2nd dim
        if length(Samples)>1
            % get single value according to Measure
            switch Measure
                case 'min'
                    SegLengths(iSeg) = min(Length_t,[],'omitnan');
                case 'max'
                    SegLengths(iSeg) = max(Length_t,[],'omitnan');
                case 'median'
                    SegLengths(iSeg) = median(Length_t,'omitnan');
                otherwise % 'mean'
                    SegLengths(iSeg) = mean(Length_t,'omitnan'); 
            end
        else
            SegLengths(iSeg) = Length_t;
        end
    end
end
% average across the segments that have been calculated
SegmentLength = mean(SegLengths,'omitnan');

% if the segment length was not calculated, use the default length
if isnan(SegmentLength)
    % get the first entries of the segments (base)
    SegmentNos = sum(DefBody.segments.id == SegmentIds(:,1)', 2,'native');
    % get the average
    SegmentLength = mean(DefBody.segments.size(SegmentNos)) * DefBody.general.height;
    IsDefault = true;
end
end
