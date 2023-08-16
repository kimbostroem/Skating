function [ Force, Frequency, COP, Moment, Location, Analog, BoardNames ] = kraftAusQualisys( Data, Scale, ForcePlateLabels )
%% Extract the force signals from a QTM mat structure
% Returns ground ACTION forces, i.e. the vertical component is in negative-Z when loaded
% 
% SYNTAX
% Force = kraftAusQualisys(Data);
% [Force, Frequency, COP, Moment, Location, Analog, ForceBoard] = kraftAusQualisys(Data,Scale,ForcePlateLabels);
%
% INPUT
%     Data          	(QTM struct) imported data structure according to Qualisys QTM format
%     Scale         	(double) scale parameter for COP values (default 1; e.g. 0.001 to scale from mm to m)
%     ForcePlateLabels	(char cell array) optional list of force plate labels to try import
%            
% OUTPUT
%     Force        	(NPlates x 3D x NSamples double) imported force data
%     Frequency 	(double) measurement frequency of the force plates
%     COP       	(NPlates x 3D x NSamples double) imported COP data
%     Moment    	(NPlates x 3D x NSamples double) imported force moment data
%     Location   	(NPlates x 3D double x 4) imported locations of the force plate corners
%     Analog     	(NChannels x NSamples) imported analog signals of the force plates
%     BoardNames	(cell array) BoardNames in Data.Analog structure
%
% EXAMPLE 
% see syntax above; note that just one input parameter is mandatory
%
% Local functions: none
%
% See also: datenAusQualisys, EMGAusQualisys, rigidbodyAusQualisys, loadc3d, extractForcePlates
% 
% (c) 2019 by Movement Science, WWU Muenster
% Author: Marc de Lussanet, WWU Muenster
% Version 11, 220908 (DK, MdL) FIXED : overwritten Variable; init IsForcePlate as logical; replace
% i,j with meaningfull variable names
% Version 230214 (MdL) more robust check for the sign of the force
% Version 230221 (MH) removed isMedianPositive Check
% Version 230419 (MdL) check if field ForcePlateOrientation exists (in old files it may be absent)
% Version 230421 (MdL) ignore extra analog channels
% Version 230516 (MdL) give BoardNames a value
% Version 230522 (MdL) fixed unknown variable Forces -> Force

narginchk(1,3);
if nargin < 2 || ~Scale || isempty(Scale),  Scale=1; end
if nargin < 3 || isempty(ForcePlateLabels)
    ForcePlateLabels = {'Force-plate','Force plate','Force platform',''};
end
% check if force data are present
IsForce   = 0;
Force     = [];
Frequency = nan;
COP       = [];
Moment    = [];
Location  = [];
Analog    = [];
BoardNames= {''};
Y=2;Z=3;XY=1:2; %#ok<NASGU>

% parameters for force data
if isfield(Data,'Force') && ~isempty(Data.Force)
    IsForce  = 1;
    NSamples = cell2mat({Data.Force.NrOfSamples});
    Frequency= cell2mat({Data.Force.Frequency});
    if ~all(Frequency == Frequency(1)) || ~all(NSamples == NSamples(1))
        warning('The force plates are measured at different frequencies or no of samples!');
    end
    Frequency = Frequency(1);
    NSamples  = NSamples(1);
    % some items may not be a force plate (e.g. ergo bike forces)
    IsForcePlate = false(size(Data.Force));
    for ForceElement = 1:length(IsForcePlate)
        if contains(Data.Force(ForceElement).ForcePlateName,ForcePlateLabels,'IgnoreCase',true)
            IsForcePlate(ForceElement) = true;
        end
    end
    NPlates  = sum(IsForcePlate);
    Force    = NaN(NPlates,Z,NSamples);
    COP      = NaN(NPlates,Z,NSamples);
    Moment   = NaN(NPlates,Z,NSamples);
    Location = NaN(NPlates,4,Z); % locations of the four corners
    % Sign of the force: when the flag is 1, we have Ground Action forces, if 0, we have Ground
    % Reaction forces (in the latter case, multiply by -1). In the wording by QTM: "Coordinate system 
    % in which force data is expressed: 0 (local force plate coordinates), 1 (global coordinate system)"
    if isfield(Data.Force,'ForcePlateOrientation')
        Sign = [Data.Force.ForcePlateOrientation] * 2 - 1;
    else % in old versions of QTM, this filed is not yet defined
        Sign = ones(size(Data.Force));
    end
    ForcePlateCounter = 0;
    for ForceElement = 1:length(IsForcePlate)
        if IsForcePlate(ForceElement) == true
            ForcePlateCounter = ForcePlateCounter + 1;
            Force(   ForcePlateCounter,:,:) = Data.Force(ForceElement).Force * Sign(ForceElement);
            if isfield(Data.Force(ForceElement),'COP')
                COP( ForcePlateCounter,:,:) = Data.Force(ForceElement).COP * Scale;
            end
            if isfield(Data.Force(ForceElement),'Moment')
                Moment(ForcePlateCounter,:,:) = Data.Force(ForceElement).Moment;
            end
            Location(ForcePlateCounter,:,:) = Data.Force(ForceElement).ForcePlateLocation * Scale;
        end
    end
end
% Plausibility check
IsReactionForces = isReactionForce(Force);
if IsReactionForces
    % in c3d files, the extraction of the sign of force is not reliable, so 
    % change the sign automatically
    if endsWith(Data.File,'.c3d')
        Force = -Force;
        warning('kraftAusQualisys : the sign of the forces is unplausible: reversing the sign');
    else
        warning('kraftAusQualisys : the sign of the forces seems unplausible');
    end
end

% parameters of raw analog channel data
if isfield(Data,'Analog') && isfield(Data.Analog,'BoardName')
    BoardNames = {Data.Analog.BoardName};
    for BoardNo = 1:size(Data.Analog,2)
        BoardName = BoardNames{BoardNo};
        ForcePlateItems = find(contains(BoardNames,{'Kistler','AMTI'},'IgnoreCase',true));
        NrOfChannels = 0;
        for ForceElement=1:length(ForcePlateItems)
            NrOfChannels = NrOfChannels + Data.Analog(ForcePlateItems(ForceElement)).NrOfChannels;
        end
        if contains(BoardName,{'USB-2533','Kistler','AMTI'})
            if     contains(BoardName,'Kistler'),  NChannelsPerPlate = 8;
            elseif contains(BoardName,'USB-2533'), NChannelsPerPlate = 6;
            elseif contains(BoardName,'AMTI'),     NChannelsPerPlate = 6;
            else,  continue;
            end
            if ~IsForce
                NSamples  = Data.Analog(BoardNo).NrOfSamples;
                warning('Analog channels but no computed force data found. Check the data!');
            end
        else
            continue;
        end
        % the number of channels is usually all information that is provided, but in some cases,
        % extra information is saved on the same analog board. Ignore these additional channels by
        % checking that the board cannot contain more plates than are present in the Forces section
        NPlatesInThisBoard = min(NPlates, Data.Analog(BoardNo).NrOfChannels / NChannelsPerPlate);
        AnalogBrd = NaN(NPlatesInThisBoard,NChannelsPerPlate,NSamples);
        for ForceElement = 1:NPlatesInThisBoard
            AnalogBrd(ForceElement,:,:)  = Data.Analog(BoardNo).Data((ForceElement-1)*NPlatesInThisBoard+(1:NChannelsPerPlate),:);
        end
        Analog = [Analog; AnalogBrd]; %#ok<AGROW> 
    end
end
end