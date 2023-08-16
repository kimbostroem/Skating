function Stream = extractErgoForces(Stream, Data, RefBody, Info, DoInsert_IsContact)
%% fill forces from ergometer contact points
% The fully instrumented ergometer measures contact forces on steer (each hand), saddle and pedals,
% as well as the orientation and angular speed of the pedals. The data are processed by the script
% "Auswertung_Ergomessungen.m" and exported as QTM-type Matlab structure.
%
% SYNTAX
% Stream = extractErgoForces(Stream, Data, RefBody, Info, DoInsert_IsContact);
%
% INPUT
%     Stream  (MNData struct) Extracted Data stream to which ergometer forces are added
%     Data    (QTM struct) Imported Data possibly containing ergometer forces
%     RefBody (Def_Body struct) struct of tables imported from Def_Body.xlsx
%     Info    (struct) containing data information such as time
%     DoInsert_IsContact (logical)
%
% OUTPUT
%     Stream  (MNData struct) Imported Data stream to which ergometer forces are added
%
% See also: extractForcePlates
%
% (c) 2020 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% Version 221219 (MdL) restructured and organized

% Init
narginchk(4,5);
if nargin<5 || isempty(DoInsert_IsContact)
    DoInsert_IsContact = false;
end

% check for force elements : everything that is not a force plate
if isfield(Data,'Force') && ~isempty(Data.Force)
    ForceElementNames = {Data.Force.ForcePlateName};
    %  skip force plates
    IsForceElement = contains(ForceElementNames,Info.ForcePlateLabels,'IgnoreCase',true);
    IsErgoForce = ~IsForceElement & ~cellfun(@isempty,ForceElementNames);
else
    IsErgoForce = false;
end

if ~any(IsErgoForce)
    return;
end

% copy and resample the force elements
QTMTimeSRint = Info.QTMTimeSRint;
NSinternal   = length(QTMTimeSRint);
NForceElements = size(Data.Force,2);
ErgoNames = cell(NForceElements,1);
% Number of samples of the contact force points
NErgoSamp = min([Data.Force(IsErgoForce).NrOfSamples]); % (min: the number can vary slightly)
ErgoForce = NaN(sum(IsErgoForce),3,length(QTMTimeSRint));
ErgoElementNo = 0;
for ForceElementNo=1:NForceElements
    if IsErgoForce(ForceElementNo)
        ErgoElementNo = ErgoElementNo+1;
        ErgoNames(ErgoElementNo) = {Data.Force(ForceElementNo).ForcePlateName};
        TimeOffset = 0;
        TimeForce = (0:(NErgoSamp-1)) / Data.Force(ForceElementNo).Frequency + TimeOffset;
        for cc=1:3
            ErgoForce(ErgoElementNo,cc,:) = ...
                interp1(TimeForce,squeeze(Data.Force(ForceElementNo).Force(cc,1:NErgoSamp)),QTMTimeSRint);
        end
    end
end
ErgoNames(cellfun(@isempty,ErgoNames)) = [];

% assign the forces to the contact points
SubjectBody = RefBody;
Scale = SubjectBody.general.height;
RFootLen = SubjectBody.segments.size(matches(SubjectBody.segments.name,'right_foot_segment')) * Scale;
LFootLen = SubjectBody.segments.size(matches(SubjectBody.segments.name, 'left_foot_segment')) * Scale;
RHandLen = SubjectBody.segments.size(matches(SubjectBody.segments.name,'right_hand_segment')) * Scale;
LHandLen = SubjectBody.segments.size(matches(SubjectBody.segments.name, 'left_hand_segment')) * Scale;

% initialize the forces substructure, if it is not yet present (one of the other force
% extraction functions may already have initialised it)
ForceItemNames  = SubjectBody.forces.name;
if isfield(Stream,'forces')
    CumItemNo = length(Stream.forces);
else
    Stream.forces = [];
    CumItemNo = 0;
end

% loop the ergo-force items and add them to stream
for ItemNo = 1:length(ErgoNames)
    forceItem = ErgoNames{ItemNo}; % the current contact segment
    % Warn if the stream already contains this item
    for it=1:length(Stream.forces)
        if strcmp(forceItem,Stream.forces(it))
            warning('extractErgoForces: force item "%s" already exists',forceItem);
        end
    end
    % ErgoNames lists the contact points that are present in ErgoForce
    ErgoMatchNo = matches(ForceItemNames,ErgoNames{ItemNo});
    if ~any(ErgoMatchNo)
        %skip if there is no match
        warning('extractErgoForces: no matching force item found to "%s"',forceItem)
        continue;
    end

    % initiate
    CumItemNo = CumItemNo + 1;
    ID = SubjectBody.forces.id(matches(SubjectBody.forces.name,forceItem));
    Stream.forces(CumItemNo).id = ID;
    Stream.forces(CumItemNo).name = forceItem;
    Stream.forces(CumItemNo).dynamic.force.units = 'N';
    Stream.forces(CumItemNo).dynamic.force.data  = squeeze(ErgoForce(ItemNo,1:3,:))';
    Stream.forces(CumItemNo).dynamic.COP.units = 'm';
    Stream.forces(CumItemNo).dynamic.COP.data  = zeros(NSinternal,3);
    if  DoInsert_IsContact
        Stream.forces(CumItemNo).dynamic.isContact.units = 'bool';
        Stream.forces(CumItemNo).dynamic.isContact.data  = true(NSinternal,1);
    end

    % define the correct contact location on the local coordinate system
    switch forceItem
        case 'right_foot_force'
            Stream.forces(CumItemNo).dynamic.COP.data(:,1) = 0.05;
            Stream.forces(CumItemNo).dynamic.COP.data(:,3) = 0.8 * RFootLen;
        case 'left_foot_force'
            Stream.forces(CumItemNo).dynamic.COP.data(:,1) = 0.05;
            Stream.forces(CumItemNo).dynamic.COP.data(:,3) = 0.8 * LFootLen;
        case 'right_hand_force'
            Stream.forces(CumItemNo).dynamic.COP.data(:,2) = 0.03;
            Stream.forces(CumItemNo).dynamic.COP.data(:,3) = 0.8 * RHandLen;
        case 'left_hand_force'
            Stream.forces(CumItemNo).dynamic.COP.data(:,2) = -0.03;
            Stream.forces(CumItemNo).dynamic.COP.data(:,3) =  0.8 * LHandLen;
        case 'right_knee_force' % leave force at joint position
        case 'left_knee_force'  % leave force at joint position
        case 'belt_force'       % leave force at joint position
        case 'pelvis_force'     % leave force at joint position
        otherwise
            error('Unknown force %s',forceItem);
    end
end
    
% if only foot elements are present, calculate a hip force to compensate those forces and the
% subject's weight
if all(matches(ErgoNames,{'left_foot_force','right_foot_force'}))
    fprintf('Adding saddle forces');
    G = 9.81;
    Weight = RefBody.general.weight;
    CumItemNo = CumItemNo + 1;
    Stream.forces(CumItemNo).id = 8;
    Stream.forces(CumItemNo).name = 'pelvis_force';
    Stream.forces(CumItemNo).dynamic.force.units = 'N';
    Stream.forces(CumItemNo).dynamic.COP.units  = 'm';
    Stream.forces(CumItemNo).dynamic.COP.data   = zeros(NSinternal,3);
    Stream.forces(CumItemNo).dynamic.force.data = ...
        - Stream.forces(CumItemNo-1).dynamic.force.data ...
        - Stream.forces(CumItemNo-2).dynamic.force.data;
    Stream.forces(CumItemNo).dynamic.force.data(:,3) = ...
        Stream.forces(CumItemNo).dynamic.force.data(:,3) + G * Weight;
end
end

