function [Prefix,LabelsClean] = getQTMPrefix(LabelNames,ValidLabelNames,DesiredPrefix,GiveFeedback)
%% Get the Qualisys Skeleton-Prefix or Vicon-Prefix from the markers
% The QTM/Vicon prefix contains of all characters of a label before the first underscore or colon ({'_', ':'})
%
% input:
%     LabelNames      [cell array] labels of the marker data
%     ValidLabelNames [cell array] all valid marker names
%     DesiredPrefix   [cell string] optional if set, fill only labels of this prefix into LabelsClean
%     GiveFeedback    [logical] 
%     
% output:
%     Prefix      [cell array] list of the removed prefixes
%     LabelsClean [cell array] list of all valid labels without the prefixes
% 
%  Local functions: prefixStatistics
% 
% (c) 2021 by Predimo GmbH
% version 221110 (MdL) added plausibility check and some further cleaning (orthosportslab conflict)
% version 221128 (MdL) cherry picked from Marc.combine-Ergobike; 
% remove detected prefix also from non-valid labels

if nargin<3 || contains(DesiredPrefix, '*') || isempty(DesiredPrefix)
    DesiredPrefix = {''};
end
if nargin<4
    GiveFeedback = false;
end

% If a specific MarkerPrefix is required, then only remove this one
if ~cellfun('isempty',DesiredPrefix)
    AllowOnlyDesiredPrefix = true;
else
    AllowOnlyDesiredPrefix = false;
end

% init the list of prefixes and the selection criteria
PotentialPrefixes    = cell(size(LabelNames));
WithPrefixIsValid    = false(size(LabelNames));
WithoutPrefixIsValid = false(size(LabelNames));
% The prefix is set before the label and separated by an underscore '_' or colon ':'
for LabelNo=1:length(LabelNames)
    PotentialPrefixes(LabelNo) = {''};
    if ~isempty(contains(LabelNames{LabelNo},{'_',':'})) % if there is at least one '_' or ':'
        [Pr,Label] = strtok(LabelNames{LabelNo},{'_',':'}); % separate by the '_' or ':'
        % Is the label with potential prefix a valid label name?
        WithPrefixIsValid(LabelNo) = any(strcmp(ValidLabelNames,LabelNames{LabelNo}));
        % Is the rest (without the {'_',':'}) a valid label name?
        WithoutPrefixIsValid(LabelNo) = any(strcmp(ValidLabelNames,Label(2:end)));
        if WithoutPrefixIsValid(LabelNo)
            PotentialPrefix = any(contains(DesiredPrefix,Pr));
            if AllowOnlyDesiredPrefix
                % If a specific MarkerPrefix is required, then only remove this one
                if PotentialPrefix
                    PotentialPrefixes(LabelNo) = {Pr};
                end
            else
                % add the prefix to the list
                PotentialPrefixes(LabelNo) = {Pr};
            end
        end
    end
end
% Check potential prefixes have been detected and how many of each
[DetectedPotentialPrefixes,NDetected,NEmptyPrefixes] = prefixStatistics(PotentialPrefixes);

%% Resolve conflicts
% If a label is valid both with and without its potential prefix there is a conflict
Conflicts = WithPrefixIsValid & WithoutPrefixIsValid;
for LabelNo=1:length(LabelNames)
    % In case a conflict exists, retain the most plausible option
    if Conflicts(LabelNo)
        % the number of Labels that have the same potential prefix as the current one (if there are
        % no others, it is unlikely that this is a true prefix)
        N_InThisPrefix = NDetected(matches(DetectedPotentialPrefixes,PotentialPrefixes(LabelNo)));
        % is the prefix plausible?
        ThePrefixIsPlausible = N_InThisPrefix>2 && N_InThisPrefix > NEmptyPrefixes;
        % if not: reject it
        if ThePrefixIsPlausible == false
            PotentialPrefixes(LabelNo) = {''};
        end
    end
end

% Get the prefix and check that there is just one prefix (multiple subjects in a single recording are currently not supported)
if GiveFeedback
    [DetectedPotentialPrefixes,NDetected] = prefixStatistics(PotentialPrefixes);
    NonEmptyPrefixes = DetectedPotentialPrefixes(~cellfun('isempty',DetectedPotentialPrefixes));
    if length(NonEmptyPrefixes)>1
        fprintf('More than one subject was detected according to the prefixes: %s\n',strjoin(NonEmptyPrefixes,', '));
    end
    if all(cellfun('isempty',DetectedPotentialPrefixes))
        %fprintf('No marker prefixes were detected\n')
    else
        for PrefixNo=1:length(DetectedPotentialPrefixes)
            fprintf('%2d Labels have prefix "%s"\n',NDetected(PrefixNo),DetectedPotentialPrefixes{PrefixNo});
        end
    end
end

% Create the list of clean labels
LabelsClean = LabelNames; % default = do not remove a prefix
RemovedPrefixes = PotentialPrefixes;
Prefix = {}; Nprx = 0;
for LabelNo=1:length(LabelNames)
    % A prefix may have been recognized on markers with a valid label, i.e. ones that exist in the markers table
    if ~isempty(RemovedPrefixes{LabelNo})
        if isempty(Prefix)
            % first marker on which a prefix has been detected: add the prefix to the list
            Nprx      =  Nprx+1;
            Prefix(1) = RemovedPrefixes(LabelNo);
        elseif ~any(strcmp(Prefix, RemovedPrefixes{LabelNo}))
            % on further markers: if the prefix is novel, add it to the list
            Nprx      =  Nprx+1;
            Prefix(Nprx) = RemovedPrefixes(LabelNo); %#ok<AGROW> 
        end
        % get the labelname without the prefix
        [~,LabelsClean{LabelNo}] = strtok(LabelNames{LabelNo},{'_',':'});
        LabelsClean{LabelNo}(1) = '';
    end
end
% Remove prefixes also from "non-valid" labels, i.e. ones that exist in the markers table
for LabelNo=1:length(LabelNames)
    if isempty(RemovedPrefixes{LabelNo})
        PotentialPrefix = strtok(LabelNames{LabelNo},{'_',':'});
        HasPrefix = matches(PotentialPrefix,Prefix);
        if HasPrefix
            % remove the prefix
            [~,LabelsClean{LabelNo}] = strtok(LabelNames{LabelNo},{'_',':'});
            LabelsClean{LabelNo}(1) = '';
        end            
    end
end
end

function [DetectedPotentialPrefixes,NDetected,NEmptyPrefixes] = prefixStatistics(PotentialPrefixes)
%% Check potential prefixes have been detected and how many of each

DetectedPotentialPrefixes = unique(PotentialPrefixes);
NDetected = zeros(length(DetectedPotentialPrefixes),1);
for i=1:length(DetectedPotentialPrefixes)
    NDetected(i) = sum(matches(PotentialPrefixes,DetectedPotentialPrefixes{i}));
end
EmptyPrefixNo = cellfun('isempty',DetectedPotentialPrefixes);
if any(EmptyPrefixNo)
    NEmptyPrefixes = NDetected(EmptyPrefixNo);
else
    NEmptyPrefixes =0;
end
end

