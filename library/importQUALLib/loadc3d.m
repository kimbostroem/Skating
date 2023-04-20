function [Data, Header] = loadc3d(FullFileName)
%% loadc3d - Retrieve 3D coordinate/analog data from a C3D file. 
% Reads: markers, analog data, and segments (markerless).
% 
% A detailed c3d documentation can be found in the C3D_User_Guide: 
%     https://www.c3d.org/c3ddocs.html
%     https://www.c3d.org/docs/C3D_User_Guide.pdf
% 
% The binary c3d file consists of 512 byte chunks. The header, containing rudimentary information
% is read; the main information is contained by the parameter section and is stored in the
% - Data.ParameterInfo field. This maintains original stucture and copies all information from the
% input file.
% - All data is stored in QUALISYS-compatible mat data structure (Data). 
% The current script is really fast.
% 
% Tested for: Captury, Qualisys QTM, Theia Markerless, Vicon Nexus, EZC3D, OptiTrack Motive files.
%
% INPUT
%     FullFileName   - (char) file (including path) to be read
%
% OUTPUTS
%     Data           - (struct) All data, according to QTM matlab structure, 
%                               with additional fields Meta and ParameterInfo
%     Header         - (struct) C3D header (this is very minimalistic information)
% 
% local functions:
%     getC3DFileHandle, readParameterSection, defineMeta, defineTrajectories, defineAnalog, 
%     defineSegments, defineEvents, readC3DdataBlock, checkDataSize, restructureAnalog, 
%     restructureCalculated, defineForce, getTara, getCOP
% 
% See also: kraftAusQualisys, getCOPfromAnalog
%
% AUTHOR AND VERSION-HISTORY
% Ver. 1.0 Alan Morris, Toronto, October 1998 [originally named "getc3d.m"]
% Ver. 2.0 Jaap Harlaar, Amsterdam, april 2002 [readc3d.m]
% Ver. 3.0 Adrian Honc, Hagen, april 2017 [readc3d.m]
% Ver. 4.4 Marc de Lussanet, University of Muenster, 2022
%          Completely revised; for win and mac; automatic recognition of data content and generalizations.
%          Sped-up by reading larger chunks of data
%          Tested for Captury, Qualisys QTM, Theia markerless, Vicon Nexus, EZC3D files.
% Subversion 221004 (MdL) add EZC3D manufacturer and version
% Subversion 221027 (MdL) corrected Data.Segments.Units
% Subversion 221212 (MdL) Workaround for Theia-parameter section bug; copy units without trailing space;
% sign of force plate loading; handle gaps in Theia segments; check for the presence of segments
% Subversion 230104 (MdL) correct sign of the force (ForcePlateOrientation); handle aray of "units"
% Subversion 230227 (MdL) detect sign of forces automatically, since it is not reliable; FIXED
% return if file could not be opened
% Subversion 230305 (MdL) revised and tested getCOP
% Subversion 230308 (MdL) deal with missing DESCRIPTIONS in OptiTrack data; deal with empty BoardName field (theia c3d)
% Subversion 230320 (MdL) some fixes for dealing with units

%% Init
Data = struct;

%% Marker Dimension (marker data have an RMS value in the fourth dimension)
MarkerDimension = 4;

%% Open the file
[FId,Header] = getC3DFileHandle(FullFileName);
if FId == -1
    return
end

%% Read the file information block
ParameterInfo      = readParameterSection(FId,Header);
Data.ParameterInfo = ParameterInfo;
Data.File          = FullFileName;
if isempty(ParameterInfo(1).group)
    warning('no valid parameter section found','backtrace'); %#ok<CTPCT> 
    return;
end

%% Analyse the file info 
% analyse some of the important parameter groups
ManufacturerGroup = matches({ParameterInfo.group},'MANUFACTURER');
PointGroup        = matches({ParameterInfo.group},'POINT');
AnalogGroup       = matches({ParameterInfo.group},'ANALOG');
ForceGroup        = matches({ParameterInfo.group},'FORCE_PLATFORM');
RotGroup          = matches({ParameterInfo.group},'ROTATION');
EventGroup        = matches({ParameterInfo.group},'EVENT');

%% Meta (incl Manufacturer)
Data = defineMeta(Data,ParameterInfo,PointGroup,ManufacturerGroup,MarkerDimension);

%% Definition of marker data from the POINT group
[Data,Scale] = defineTrajectories(Data,ParameterInfo,Header,PointGroup);

%% Analog channels
Data = defineAnalog(Data,ParameterInfo,AnalogGroup);

%% Segments
[Data, NSegments, SegmentStart, NSegmentFrames] = defineSegments(Data,ParameterInfo,RotGroup);

%% Events
Data = defineEvents(Data,ParameterInfo,EventGroup);

%% Read data block
Data = readC3DdataBlock(Data, FId,Scale,NSegmentFrames,NSegments,SegmentStart);

%% Force plates
Data = defineForce(Data,ParameterInfo,ForceGroup,AnalogGroup);

%% Restructure Data struct to handle calculated and analog parameters
Data = restructureAnalog(Data);
Data = restructureCalculated(Data,ParameterInfo,PointGroup);
end 


%% =================================================================================================
%% =================================================================================================

function [FId,Header] = getC3DFileHandle(FullFileName)
%%  Open the file and return the file handle and the header
% The header is a binary section of 512 bytes. It allows a very minimalistic representation of the
% structure of the binary data. 
% The most important information in the header is the kind of endian representation and the position
% of the Parameter SecTion
% The header is read here for completeness, but most of the information is
% also stred in the Parameter Section that follows the header
% 
% INPUT
%    FId     - (double) File identifier
%    Header  - (struct) Minimal header of c3d file
%
% OUTPUT
%    FullFileName - (char) path and name of the imported C3D file
% 
% (c) 2022 by Predimo
% 220903 (MdL) Added header

[~,FileName,Ext] = fileparts(FullFileName);

% Init
Header = struct;
FId = fopen(FullFileName,'r','n'); % native format
if FId == -1
    warning(['File: ', FileName,Ext,' could not be opened']);
    return
end

% Read the file header, as described in the c3d manual
H.PointerToParameterSection = fread(FId,1,'int8'); % Byte 1: A pointer to the first block of the parameter section
H.DataStorageFlag        = fread(FId,1,'int8');    % Byte 2: A flag defining the Data section storage format.
H.N3DPoints              = fread(FId,1,'int16');   % The number of 3D points (markers) stored in each 3D frame.
H.NAnalogSamplesPerFrame = fread(FId,1,'int16');   % The total number of analog samples stored in each 3D frame.
H.FirstFrameTransferred  = fread(FId,1,'int16');   % The first frame number of raw data transfered to the C3D file
H.LastFrameTransferred   = fread(FId,1,'int16');   % The last frame number of raw data transfered to the C3D file
H.MaxInterpolationGap    = fread(FId,1,'int16');   % The maximum 3D frame interpolation gap.
H.FloatingPointFactor    = fread(FId,1,'int32');   % The floating-point factor that scales all 3D data values into system measurement units.
H.PointerToDataSection   = fread(FId,1,'int16');   % A pointer to the first block of the data storage section.
H.AnalogSampleRate       = fread(FId,1,'int16');   % The analog sample rate per 3D frame.
H.FrameRate              = fread(FId,1,'single');  % The 3D frame rate in Hz (32-bit floating-point).
NotUsed1                 = fread(FId,137,'int16'); %#ok<NASGU> % Currently not used
H.Is4CharEventLabels     = fread(FId,1,'int16');   % A key value indicating the file supports 4 char event labels.
H.NEventTimes            = fread(FId,1,'int16');   % Number of defined time events (O to 18)
NotUsed2                 = fread(FId,1,'int16');   %#ok<NASGU> % Currently not used
H.EventTimes             = fread(FId,18,'single'); % Event times (floating-point) in seconds (up to 18 events).
H.EventDisplayFlags      = fread(FId,18,'int16');  % 18 bytes - event display flags Ox00 = ON, Ox01 = OFF.
NotUsed3                 = fread(FId,1,'int16');   %#ok<NASGU> % Currently not used
H.EventLabels4Char       = fread(FId,72,'*char');  % Event labels. Each label is 4 characters long
NotUsed4                 = fread(FId,22,'int16');  %#ok<NASGU> % Currently not used
Header = H;

% Check for C3D format (the second byte of C3D files is always (int) 80 / (char) ('P'))
if Header.DataStorageFlag ~= 80
    warning(['File: ',FileName,Ext,' does not comply to the C3D format']);
    FId = [];
    return
end

% Processor type — Order of reading or writing bytes: 
% 1(INTEL-PC); 2(DEC-VAX); 3(MIPS-SUN/SGI)
fseek(FId,512*(H.PointerToParameterSection-1)+3,'bof'); % jump to processortype - field
ProcessorType = fread(FId,1,'int8')-83;
if ProcessorType == 2
    fclose(FId);
    FId = fopen(FullFileName,'r','l'); % DEC VAX D floating point and VAX ordering
elseif ProcessorType == 3
    fclose(FId);
    FId = fopen(FullFileName,'r','b'); % MIPS-SUN/SGI
end
end

%% =================================================================================================

function ParameterInfo = readParameterSection(FId,Header)
%% Read the parameter information section
% The parameter section is a sort of extended header (probably this is an extension to the original
% header which is very restricted; the parameter section duplicates and overrules parts of the
% header information).
% Parameters are grouped, each parameter and group have a name and description parameters have
% further information such as strings or (arrays of) values
% 
% INPUT
%    FId     - (double) File identifier
%    Header  - (struct) Minimal header of c3d file
%
% OUTPUT
%    ParameterInfo - (struct) containing grouped information about all parameters (i.e. the "large header")
% 
% (c) 2022 by Predimo
% 221129 (MdL) Workaround Theia-files: apparently, the header is over when a group-or-parameter name
% is empty

% Init
ParamNumber        = 0;
ParameterInfo      = struct('group',[],'description',{},'parameter',[],'data',[],'param_description',{});

% set pointer to 3rd byte of the parameter section 
ParamSectionStart  = 512*(Header.PointerToParameterSection-1);
SkipBytes          = 2;  % see manual
fseek(FId, ParamSectionStart + SkipBytes, 'bof');
NParameterBlocks   = fread(FId,1,'int8');    % negative number in some Vicon files [MdL 220827: Why?]
ParamSectionLength = 512 * abs(NParameterBlocks); 
Machinefmt         = fread(FId,1,'int8')-83; %#ok<NASGU> 

% loop for the hierarchical information table
NextRecordStart = 0;
GroupName       = {};
NChar           = 1; % init as workarond for EZC3D
while NextRecordStart < ParamSectionStart + ParamSectionLength 
    NChars        = abs(fread(FId,1,'int8'));  % characters in group/parameter name
    GroupId       = fread(FId, 1, 'int8');     % negative: group; pos: parameter of group
    IsGroupItem   = GroupId<0;
    GroupId       = abs(GroupId);              % index of the group
    Name          = fread(FId, NChars, '*char')';  % name of parameter or group
    PointerOffset = fread(FId, 1, 'int16');    % offset to next group or parameter in bytes
    PointerPosition = ftell(FId)-2;            % current file position

    EndOfParameterSection = PointerOffset == 0; % see manual
    if EndOfParameterSection || NChars == 0
        break;

    elseif IsGroupItem
        GroupName{GroupId}  = Name; %#ok<AGROW> 
        NCharDescr          = fread(FId, 1, 'uint8');    % n char in group description 
        GroupDescr{GroupId} = fread(FId, NCharDescr, '*char')'; %#ok<AGROW> 

    else % parameter data
        ParamNumber   = ParamNumber+1;
        ParamName     = Name;
        PrecisionCode = fread(FId, 1 ,'int8');

        % Precision of data: -1=char / 1=byte / 2=integer*2 / 4=real*4
        switch PrecisionCode
            case -1, Precision = '*char';
            case 1,  Precision = 'int8';
            case 2,  Precision = 'int16';
            case 4,  Precision = 'float32';
            otherwise
                error('Non-defined value for Precision (%d)',PrecisionCode);
        end
        ParameterInfo(ParamNumber).group       = GroupName{GroupId};  
        ParameterInfo(ParamNumber).description = GroupDescr{GroupId}; 
        ParameterInfo(ParamNumber).parameter   = ParamName;

        % Dimensional structure of the current parameter (0-7 dimensions are allowed)
        NDimensions = fread(FId, 1 ,'int8'); 
        NumberOfElements = 1;
        Dimension   = 1;
        if NDimensions
            Dimension = zeros(1,NDimensions);
            for jj = 1:NDimensions
                % According to the Manual, the number of characters (of the parameter name) can be 
                % negative, to indicate that the length is locked and should not be changed. 
                % HOWEVER: as int8 has a range of just [-127 128], many files use unsigned
                % integers here ([0 255]). 
                Dimension(jj) = fread(FId,1,'uint8');
                NumberOfElements = NumberOfElements*Dimension(jj);
            end
            % If there is just one field per definition, Dimension(2) is not yet set
            if NDimensions == 1
                Dimension(2) = 1;
            end
        end

        % copy into the Contents structure
        if strcmp(Precision, '*char')
            % the number of characters in the field is always in the first location of Dimension
            % (EZC3D does not redefine the number of characters for subsequent params, in which case
            % dimension(1) is zero)
            if Dimension(1) > 0
                NChar = Dimension(1); %length of character word
            end
            Dimension(1) = 1;
            NumberOfWords = NumberOfElements / NChar;
            % Create a cell with the remaining dimensions (i.e. without the first which encodes the 
            % number of characters) and read from the file
            Words = squeeze(cell(Dimension));
            for Word = 1:NumberOfWords
                % Read NChar and trim the trailing white space
                Words(Word) = {strtrim(fread(FId, NChar, '*char')')}; 
            end
            ParameterInfo(ParamNumber).data = Words; %strtrim: Remove leading and trailing whitespace
        else % read numerical
            % Length of data record for multi-dimensional array
            Datalength  = abs(PrecisionCode*NumberOfElements); 
            NDataBlocks = Datalength / abs(PrecisionCode);
            Numbers     = fread(FId, NDataBlocks, Precision);
            if all(size(Dimension)== 1), Dimension = [1 Dimension]; end %#ok<AGROW> 
            ParameterInfo(ParamNumber).data = reshape(Numbers,Dimension);
        end
        
        % Parameter description
        NCharDescr = fread(FId, 1, 'uint8');  % n char in parameter description 
        ParameterInfo(ParamNumber).param_description = {fread(FId, NCharDescr, '*char')'};
    end

    % move ahead to next record
    NextRecordStart = PointerPosition+PointerOffset;   
    fseek(FId,NextRecordStart,'bof');
end
end

%% =================================================================================================

function Data = defineMeta(Data,ParameterInfo,PointGroup,ManufacturerGroup,MarkerDimension)
%% Select the meta information and add to Data structure
% 
% INPUT
%    Data              - (struct) Input Data according to QTM matlab structure
%    ParameterInfo     - (struct) containing grouped information about all parameters (i.e. the header)
%    PointGroup        - (logical 1xn) indices of ParameterInfo that belong to the group
%    ManufacturerGroup - (logical 1xn) indices of ParameterInfo that belong to the group
%    MarkerDimension   - (double) dimension of the marker data (usually 3D or 4D)
%
% OUTPUT
%    Data         - (struct) see above
% 
% (c) 2022 by Predimo
% 221004 (MdL) add EZC3D manufacturer and version

% Define the Manufacturer substructure, if present, and write them into the Data structure
Parameters    = {ParameterInfo.parameter};
CompanyEntry  = ManufacturerGroup & matches(Parameters,'COMPANY'); % this returns a boolean array
SoftwareEntry = ManufacturerGroup & matches(Parameters,'SOFTWARE');
VersionEntry  = ManufacturerGroup & matches(Parameters,'VERSION');
Data.Meta.Company       = [];
Data.Meta.SoftwareEntry = [];
Data.Meta.Version       = [];
if any(CompanyEntry)
    Data.Meta.Company  = ParameterInfo(CompanyEntry).data{1};
elseif any(contains({ParameterInfo.group},'EZC3D'))
    Data.Meta.Company  = 'EZC3D';
end
if any(SoftwareEntry)
    Data.Meta.Software = ParameterInfo(SoftwareEntry).data{1};
end
if any(VersionEntry)
    Data.Meta.Version  = ParameterInfo(VersionEntry).data;
elseif any(contains({ParameterInfo.group},'EZC3D'))
    Data.Meta.Version  = ParameterInfo(contains({ParameterInfo.group},'EZC3D') & strcmp(Parameters,'VERSION')).data{:};
end

% Further meta informations
XscreenEntry    = PointGroup & matches(Parameters,'X_SCREEN'); % this returns a boolean array
YscreenEntry    = PointGroup & matches(Parameters,'Y_SCREEN'); 
LongFrmsEntry   = PointGroup & matches(Parameters,'LONG_FRAMES'); 
OrigLstFrmEntry = PointGroup & matches(Parameters,'ORIGINAL_LAST_FRAME'); 
if any(XscreenEntry),     Data.Meta.X_screen      = ParameterInfo(XscreenEntry).data{1}; end
if any(YscreenEntry),     Data.Meta.Y_screen      = ParameterInfo(YscreenEntry).data{1}; end
if any(LongFrmsEntry),    Data.Meta.LongFrames    = ParameterInfo(LongFrmsEntry).data;   end
if any(OrigLstFrmEntry),  Data.Meta.OrigLastFrame = ParameterInfo(OrigLstFrmEntry).data; end

% Implicit parameters for importing
Data.Meta.PointSize = MarkerDimension;
end

%% =================================================================================================

function [Data,Scale] = defineTrajectories(Data,ParameterInfo,Header,PointGroup)
%% Define the Trajectories substructure, if present, and write them into the Data structure
% Trajectory Parameters are stored in the POINT Group of the Parameter Section
% 
% INPUT
%    Data           - (struct) Input Data according to QTM matlab structure
%    ParameterInfo  - (struct) Containing grouped information about all parameters
%    Header         - (struct) The minimal header of c3d data
%    PointGroup     - (logical 1xn) indices of ParameterInfo that belong to the group
%
% OUTPUT
%    Data      - (struct) see above
%    Scale     - (double) 
% 
% (c) 2022 by Predimo
% 230103 (MdL) deal with array of units
% 230310 (MdL) define default units as m

Parameters      = {ParameterInfo.parameter};
DataStartField  = PointGroup & matches(Parameters,'DATA_START');
NInPointGroup   = ParameterInfo(PointGroup & matches(Parameters,'USED')).data;
Scale           = ParameterInfo(PointGroup & matches(Parameters,'SCALE')).data;
MarkerRate      = ParameterInfo(PointGroup & matches(Parameters,'RATE')).data;
NMarkerFrames   = ParameterInfo(PointGroup & matches(Parameters,'FRAMES')).data;
% Determine the start of the data section of the file: this is usually defined by the DATA_START
% parameter in the POINT group; alternatively it is given by the PointerToDataSection of the header
if any(DataStartField)
    StartOfDataSection  = ParameterInfo(DataStartField).data;
else
    StartOfDataSection  = Header.PointerToDataSection;
end
Data.StartOfDataSection = StartOfDataSection;
Data.Time       = ((StartOfDataSection : StartOfDataSection+NMarkerFrames-1)-1) / MarkerRate;
Data.FrameRate  = MarkerRate;
Data.Frames     = NMarkerFrames;
Data.Meta.NInPointGroup = NInPointGroup;
Dimension       = Data.Meta.PointSize; % the dimension of point data depends on the software from which was exported
if NInPointGroup % e.g. Theia markerless does not contain marker data
    MarkerLabels= ParameterInfo(PointGroup & matches(Parameters,'LABELS')).data;
    Data.Trajectories = struct('Labeled',[],'Unidentified',[]);
    Data.Trajectories.Labeled.Labels = MarkerLabels;
    Data.Trajectories.Labeled.Data   = NaN(length(MarkerLabels), Dimension, NMarkerFrames);
end
if any(PointGroup & matches(Parameters,'UNITS'))
    SizeUnits = length(ParameterInfo(PointGroup & matches(Parameters,'UNITS')).data);
    % copy the units string without trailing whitespace and null characters
    Data.Units = deblank(ParameterInfo(PointGroup & matches(Parameters,'UNITS')).data{1});
    if SizeUnits > 1 % the case in EZC3D
        UniqueUnits = unique(deblank(ParameterInfo(PointGroup & matches(Parameters,'UNITS')).data));
        if length(UniqueUnits) > 1
            warning('loadc3d : different units for different labels are not supported; using "%s"',Data.Units);
        end
    end
else
    % default unit is 'm'
    Data.Units = 'm';
    warning('loadc3d : No length units have been defined. Assuming "%s"',Data.Units);
end
end

%% =================================================================================================

function Data = defineAnalog(Data,ParameterInfo,AnalogGroup)
%% Define the Analog substructure, if present, and write them into the Data structure
% The force plate data in C3D are imported with the analog data and thus represent channels of the
% Analog substructure. 
% Combine the channels into 3D Forces and COPs
% 
% INPUT
%    Data          - (struct) Input Data according to QTM matlab structure
%    ParameterInfo - (struct) containing grouped information about all parameters (i.e. the header)
%    AnalogGroup   - (logical 1xn) indices of ParameterInfo that belong to the group
%
% OUTPUT
%    Data        - (struct) see above
% 
% (c) 2022 by Predimo
% 220830 (MdL) revised and cleaned and shortened
% 230308 (MdL) deal with missing DESCRIPTIONS (theia c3d)

Parameters      = {ParameterInfo.parameter};
IsDefinedAnalog = any(AnalogGroup) && ParameterInfo(AnalogGroup & matches(Parameters,'USED')).data>0;
if IsDefinedAnalog
    % general parameters
    NMarkerFrames   = Data.Frames;
    MarkerRate      = Data.FrameRate;

    % specific parameters
    NAnalogChan     = ParameterInfo(AnalogGroup & matches(Parameters,'USED')).data;
    AnalogScale     = ParameterInfo(AnalogGroup & matches(Parameters,'SCALE')).data;
    AnalogRate      = ParameterInfo(AnalogGroup & matches(Parameters,'RATE')).data;
    AnalogLabels    = ParameterInfo(AnalogGroup & matches(Parameters,'LABELS')).data;
    AnalogUnits     = ParameterInfo(AnalogGroup & matches(Parameters,'UNITS')).data;
    AnalogGenScale  = ParameterInfo(AnalogGroup & matches(Parameters,'GEN_SCALE')).data;
    AnalogOffset    = ParameterInfo(AnalogGroup & matches(Parameters,'OFFSET')).data;
    IsDefinedDescr  = any((AnalogGroup & matches(Parameters,'DESCRIPTIONS')));
    IsDefinedGain   = any((AnalogGroup & matches(Parameters,'GAIN')));
    if IsDefinedDescr
        AnalogDescripts  = ParameterInfo(AnalogGroup & matches(Parameters,'DESCRIPTIONS')).data;
    else
        AnalogDescripts  = cell(size(AnalogLabels));
    end
    if IsDefinedGain
        AnalogGain  = ParameterInfo(AnalogGroup & matches(Parameters,'GAIN')).data;
    else
        AnalogGain  = ones(size(AnalogOffset));
    end

    %% parameters for each channel
    % the number of analog frames that come to one marker sample
    NAnalogFramesPerSample = AnalogRate/MarkerRate;
    Data.Analog     = struct;
    for ac = 1 : NAnalogChan
        Data.Analog(ac).BoardName   = AnalogDescripts{ac};
        Data.Analog(ac).Labels      = AnalogLabels{ac};
        Data.Analog(ac).Units       = AnalogUnits{ac};
        Data.Analog(ac).NrOfSamples = NMarkerFrames * NAnalogFramesPerSample;
        Data.Analog(ac).Frequency   = AnalogRate;
        Data.Analog(ac).GenScale    = AnalogGenScale; % I do not know what VICON means with GEN_SCALE
        Data.Analog(ac).Gain        = AnalogGain(ac);
        Data.Analog(ac).Data        = zeros(Data.Analog(ac).NrOfSamples,1);
        Data.Analog(ac).Scale       = AnalogScale(ac);
        Data.Analog(ac).Offset      = AnalogOffset(ac);
    end
end
end

%% =================================================================================================

function [Data, NSegments, SegmentStart, NSegmentFrames] = defineSegments(Data,ParameterInfo,RotGroup)
%% Define the Segments substructure, if present, and write them into the Data structure
% 
% INPUT
%    Data          - (struct) Input Data according to QTM matlab structure
%    ParameterInfo - (struct) containing grouped information about all parameters (i.e. the header)
%    RotGroup      - (logical 1xn) indices of ParameterInfo that belong to the group
%
% OUTPUT
%    Data           - (struct) see above
%    NSegments      - (double) number of segments
%    SegmentStart   - (double) start of segment section in the c3d file
%    NSegmentFrames - (double) number of recorded frames in segment section
% 
% (c) 2022 by Predimo
% 221027 (MdL) corrected units entry

% segments are provided with the markersless recording, e.g., Theia.
Parameters       = {ParameterInfo.parameter};
UsedEntry        = RotGroup & matches(Parameters,'USED'); % this returns a boolean array
LabelsEntry      = RotGroup & matches(Parameters,'LABELS'); 
RateEntry        = RotGroup & matches(Parameters,'RATE'); 
DataStartEntry   = RotGroup & matches(Parameters,'DATA_START');
% Some parameters need not be added to the Data.Segments structure
%DescrptionsEntry= RotGroup & matches(Parameters,'DESCRIPTIONS'); 
%RatioEntry      = RotGroup & matches(Parameters,'RATIO'); 

if any(RotGroup) && ParameterInfo(UsedEntry).data>0
    MarkerRate    = Data.FrameRate;
    NMarkerFrames = Data.Frames;
    NSegments     = ParameterInfo(UsedEntry).data;
    SegmentLabels = ParameterInfo(LabelsEntry).data;
    SegmentRate   = ParameterInfo(RateEntry).data;
    SegmentStart  = ParameterInfo(DataStartEntry).data;
    Data.Segments = struct('Labeled',[],'Unidentified',[]);
    if SegmentRate == MarkerRate
        NSegmentFrames = NMarkerFrames;
    else
        NSegmentFrames = round(NMarkerFrames * SegmentRate / MarkerRate);
    end
    Data.Segments.Labeled.Count  = NSegments;
    Data.Segments.Labeled.Labels = SegmentLabels;
    Data.Segments.Labeled.Time   = (0:(NSegmentFrames-1))' / SegmentRate;
    Data.Segments.Labeled.Timeunits = 's';
    if isfield(Data,'Units')
        Data.Segments.Labeled.Units  = Data.Units;
    else
        Data.Segments.Labeled.Units  = 'm';
    end
    Data.Segments.Labeled.Rate   = SegmentRate;
    Data.Segments.Labeled.Rotmat = zeros(NSegments,NSegmentFrames,3,3);
    Data.Segments.Labeled.Pos    = zeros(NSegments,NSegmentFrames,3);
else
    NSegments = 0;
    SegmentStart = 0;
    NSegmentFrames = 0;
end
end

%% =================================================================================================

function Data = defineEvents(Data,ParameterInfo,EventGroup)
%% define the Events substructure and write them into the Data structure
if any(EventGroup)
    % extract the event information
    Parameters  = {ParameterInfo.parameter};
    nEvents     = ParameterInfo(EventGroup & matches(Parameters,'USED')).data;
    IsAnyTimesDefined = any(EventGroup & matches(Parameters,'TIMES'));
    if nEvents>0 && IsAnyTimesDefined
        EventLabels = ParameterInfo(EventGroup & matches(Parameters,'LABELS')).data;
        EventFrames = ParameterInfo(EventGroup & matches(Parameters,'TIMES')).data(1,:);
        EventTimes  = ParameterInfo(EventGroup & matches(Parameters,'TIMES')).data(2,:);
        % workaround for bug in QTM Export (the frame values are all zero)
        MarkerRate  = Data.FrameRate;
        if any(EventTimes~=0) && all(EventFrames==0)
            EventFrames = round(EventTimes * MarkerRate + 1);
        end
        % fill the Data structure
        for ee=1:nEvents
            Data.Events(ee).Label = EventLabels{ee};
            Data.Events(ee).Frame = EventFrames(ee);
            Data.Events(ee).Time  = EventTimes(ee);
        end
    end
end
end

%% =================================================================================================

function [Data,Log] = readC3DdataBlock(Data, FId,Scale,NSegmentFrames,NSegments,SegmentStart)
%% Read the binary data from the data block of the input file
% 
% INPUT
%     Data           - (struct) Input Data according to QTM matlab structure
%     FId            - (fid) file pointer
%     Scale          - (double) Scaling factor
%     NSegmentFrames - (double) Number of frmes for segments
%     NSegments      - (double) Number of segments (e.g. QTM skeleton)
%     SegmentStart   - (double) Pointer to start of segment block
%
% OUTPUT
%     Data         - (struct) see above
%     Log          - (struct) Log of file and data sizes
% 
% Local functions : checkDataSize
% 
% (c) 2022 by Predimo
% Version 221212 (MdL) FIXED Check for the presence of segments
% Version 230320 (MdL) FIXED Fixed hidden bug which redefined the input parameter "Scale"

NMarkerFrames = Data.Frames;
StartOfDataSection = Data.StartOfDataSection;
MarkerRate    = Data.FrameRate;
NaNRange      = -999998970; % found in theia segments
FReadSize     = 'float32';
FloatSize     = 4; % how many bytes fit in a float number
Log           = struct;
NAnalogChan   = 0;
NAnalogFramesPerSample  = 0;
% Define the sizes of analog data if present
Parameters    = {Data.ParameterInfo.parameter};
AnalogGroup   = matches({Data.ParameterInfo.group},'ANALOG');
if any(AnalogGroup)
    NAnalogChan  = Data.ParameterInfo(AnalogGroup & matches(Parameters,'USED')).data;
end
if NAnalogChan > 0
    AnalogRate   = Data.ParameterInfo(AnalogGroup & matches(Parameters,'RATE')).data;
    NAnalogFramesPerSample = AnalogRate/MarkerRate;
    % Modulo of AnalogRate/MarkerRate must be zero
    if mod(AnalogRate,MarkerRate)
        error('The analog rate (%d) must be a multiple of frame rate (%d)',AnalogRate,MarkerRate);
    end
end

% Define the size of a sample of data
NMarkers             = Data.Meta.NInPointGroup;
SizeofPoint          = Data.Meta.PointSize; % 3 or 4 for 3D or 4D
PointGroupSampleSize = NMarkers * SizeofPoint;
AnalogSampleSize     = NAnalogChan * NAnalogFramesPerSample;
SampleSize           = PointGroupSampleSize + AnalogSampleSize;
% Segments are rotation matrix (3x3 + 3D-point in a 4x4 matrix)
SizeofSegment        = 17; % 4x4+1 (i.e. there is a redundancy of 5 floats)
SegmentSampleSize    = SizeofSegment * NSegments; % -> Theia markerless

%% Read and process the data section
% Scale defines whether the data contain float or int values
if Scale < 0 % Negative scale factor means: these are "float32"-data

    % Check that the size of the data matches the filesize
    Data.Meta.Log = checkDataSize(FId, StartOfDataSection, NMarkerFrames, SampleSize, ...
        NSegmentFrames, SegmentSampleSize, FloatSize);

    %% Read the POINT and ANALOG Group data
    % Set file pointer to beginning of binary data block
    fseek(FId,(StartOfDataSection-1)*512,'bof');
    % fread the complete data section at once: this is much faster than single values
    Raw = fread(FId, NMarkerFrames * SampleSize, FReadSize);

    % Extract the POINT group and ANALOG group data from Raw
    for ii = 1:NMarkerFrames 
        OffsetMarkers = (ii-1) * SampleSize;
        OffsetAnalogs = OffsetMarkers + PointGroupSampleSize;
        if NMarkers
            % Markers: 3D position (in many cases 4D: 3D + "camera info" + residual error)
            MarkersSample = reshape(Raw(OffsetMarkers + (1 : PointGroupSampleSize)),[SizeofPoint NMarkers])';
            Truncate      = fix(MarkersSample(:, SizeofPoint));
            highbyte      = fix(Truncate/256);
            lowbyte       = Truncate - highbyte*256;
            CameraInfo    = highbyte;
            ResidualError = lowbyte*abs(Scale);
            for jj = 1:NMarkers
                Data.Trajectories.Labeled.Data(jj,1:3,ii) = MarkersSample(jj,1:3);
                % Often, markers have a fourth dimension containing the RMS and camera info
                if SizeofPoint == 4
                    Data.Trajectories.Labeled.Data(jj, 4 ,ii) = ResidualError(jj);
                    Data.Trajectories.Labeled.Filltype(jj,ii) = CameraInfo(jj);
                end
            end
        end
        if NAnalogChan
            % Analog: This is always a multiple of the sample rate of markers, so a multiple of
            % analog frames can be stored per marker frame
            AnalogSample = reshape(Raw(OffsetAnalogs + (1:NAnalogChan*NAnalogFramesPerSample)),[NAnalogChan NAnalogFramesPerSample])';
            for ac = 1:NAnalogChan
                Rng = (ii-1)*NAnalogFramesPerSample + (1:NAnalogFramesPerSample)';
                Data.Analog(ac).Data(Rng) = ...
                    Data.Analog(ac).Offset + ...
                    AnalogSample(1:NAnalogFramesPerSample,ac) * abs(Data.Analog(ac).Scale);
            end
        end
    end

    %% Read the SEGMENT group data (if present e.g. Theia markerless)
    % Set file pointer to beginning of binary data block
    fseek(FId,(SegmentStart-1)*512,'bof');
    % fread the complete data section at once
    Raw = fread(FId, NSegmentFrames * SegmentSampleSize, FReadSize);

    % Extract the SEGMENT group data from Raw
    if NSegmentFrames
        for ii = 1:NSegmentFrames
            % Segments: Rotation matrix + position in a [(4x4+1) x NJoints] block per sample
            OffSet       = (ii-1) * SegmentSampleSize;
            JointsSample = Raw(OffSet + (1:SizeofSegment*NSegments));
            JointsSample = reshape(JointsSample,[SizeofSegment, NSegments]);
            JointsSample = reshape(JointsSample(1:16,:),[4,4,NSegments]);
            for ss = 1:NSegments
                % the value stored in position (4,4) can scale the size, but usually it is one
                SegScale = JointsSample(4,4,ss);
                Data.Segments.Labeled.Rotmat(ss,ii,:,:) = squeeze(JointsSample(1:3,1:3,ss));
                Data.Segments.Labeled.Pos(ss,ii,:)      = squeeze(JointsSample(1:3,4,ss)) * SegScale;
                % the last value is the scale parameter and the other values are zero
                %Data.Segments.Labeled.Quatro(ss,ii,:)   = squeeze(JointsSample(4,1:4,ss)) * SegScale;
            end
        end
        % detect and remove gaps
        Gaps = all(Data.Segments.Labeled.Pos<NaNRange,3);
        for ss=1:NSegments
            GapsInSegment = Gaps(ss,:);
            Data.Segments.Labeled.Pos(ss,GapsInSegment,:) = nan;
            Data.Segments.Labeled.Rotmat(ss,GapsInSegment,:,:) = nan;
        end
        % figure;hold on;plot(squeeze(Data.Segments.Labeled.Pos(14,:,3)));
    end

else

    %% Positive scale factor: int-data scaled by the factor
    warning('This c3d file is not binary but int-defined: this import is not yet implemented');
    NMarkers = length(Data.Trajectories.Labeled.Labels);
    fseek(FId,(StartOfDataSection-1)*512,'bof');
    for ii = 1:NMarkerFrames
        for jj = 1:NMarkers
            Markers(ii,jj,1:3) = fread(FId,3,'int16')'.*Scale; %#ok<NASGU,AGROW> 
            ResidualError(ii,jj) = fread(FId,1,'int8'); %#ok<NASGU,AGROW> 
            CameraInfo(ii,jj) = fread(FId,1,'int8'); %#ok<NASGU,AGROW> 
        end
        for jj = 1:NAnalogFramesPerSample
            AnalogSignals(jj+NAnalogFramesPerSample*(ii-1),1:NAnalogChan) = ...
            fread(FId,NAnalogChan,'int16')'; %#ok<NASGU,AGROW> 
        end
    end

end
fclose(FId);
% Test = squeeze(Data.Trajectories.Labeled.Data(1,:,:));
% figure; plot(Test');
end

%% =================================================================================================

function Log = checkDataSize(FId, StartOfDataSection, NMarkerFrames, SampleSize, ...
        NSegmentFrames, SegmentSampleSize, FloatSize)
%% Check that the size of the data matches the filesize
% 
% INPUT
%     FId                - (file id) file handle
%     StartOfDataSection - (double) start of the data section in the c3d file
%     NMarkerFrames      - (double) number of marker frames
%     SampleSize         - (double) number of floats in a marker frame (including analog)
%     NSegmentFrames     - (double) number of segment frames
%     SegmentSampleSize  - (double) number of floats in a segment frame
%     FloatSize          - (double) number of bytes in the floating point numbers
%
% OUTPUT
%    Log         - (struct) parameters telling how much of the sile size is actiually read
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: MdL
% Version: 220903 first version

% determine size of data block:
fseek(FId,(StartOfDataSection-1)*512,'bof'); % set pointer to beginning of data section
Start               = ftell(FId); % the byte number at which the data section starts
fseek(FId,0,'eof');
Log.SizeOfDataInFile= ftell(FId)-Start;
Log.SizeOfFRead     = (NMarkerFrames * SampleSize + NSegmentFrames * SegmentSampleSize) * FloatSize;
Log.PercentImported = 100 * Log.SizeOfFRead / Log.SizeOfDataInFile;
Log.FileOversize    = Log.SizeOfDataInFile - Log.SizeOfFRead;
if Log.SizeOfDataInFile < Log.SizeOfFRead
    error('Files size (%d bytes) is smaller than the amount of data to be read (%d bytes)!',...
        Log.SizeOfDataInFile,Log.SizeOfFRead);
else
    fprintf('Reading Datablock (%.1f MB): => %.0f%% of the file, %d bytes are not read\n',...
        Log.SizeOfFRead/(2^20),Log.PercentImported,Log.FileOversize)
end
end

%% =================================================================================================

function Data = restructureAnalog(Data)
%% Restructure Data.Analog struct to groups the analog boards
% Analog boards usually have a number of channels in a logical group. Restructure the list of analog
% channels into a board list
% 
% INPUT
%    Data    - (struct) Input Data according to QTM matlab structure
%
% OUTPUT
%    Data    - (struct) see above
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: MdL
% Version: 31.8.2022 Combine boards according to QTM convention

% The Analog field may be absent
if ~isfield(Data, 'Analog')
    return;
end

% remove Scale and Offset fields since these have already been used to update the data
Data.Analog = rmfield(Data.Analog,"Scale");
Data.Analog = rmfield(Data.Analog,"Offset");
% The analog section can store data from multiple devices ("boards"). Loop the boards to group the
% channels of each board.
BoardNames = unique({Data.Analog.BoardName});
Analog = struct;
for i=1:length(BoardNames)
    Channels = find(matches({Data.Analog.BoardName},BoardNames(i)));
    Analog(i).BoardName    = BoardNames{i};
    Analog(i).NrOfChannels = length(Channels);
    Analog(i).Labels       = {Data.Analog(Channels).Labels};
    Analog(i).Units        = {Data.Analog(Channels).Units};
    Analog(i).NrOfSamples  = Data.Analog(Channels(1)).NrOfSamples;
    Analog(i).Frequency    = Data.Analog(Channels(1)).Frequency;
    Analog(i).Data         = [Data.Analog(Channels).Data]';
    if isfield(Data.Analog,'GenScale') % I do not know what VICON means with GEN_SCALE
        Analog(i).GenScale = Data.Analog(Channels(1)).GenScale;
    end
    if isfield(Data.Analog,'Gain')
        Analog(i).Gain     = Data.Analog(Channels(1)).Gain;
    end
end
Data.Analog = Analog;
end

%% =================================================================================================

function Data = restructureCalculated(Data,ParameterInfo,PointGroup)
%% Restructure Data struct to handle calculated parameters
% Vicon (and others) stores various calculated variables into the Trajectories. 
% These are moved to a separate "Calculated" substructure.
% 
% INPUT
%    Data           - (struct) Input Data according to QTM matlab structure
%    ParameterInfo  - (struct) containing grouped information about all parameters (i.e. the header)
%    PointGroup     - (logical 1xn) indices of ParameterInfo that belong to the group
%
% OUTPUT
%    Data         - (struct) see above
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: MdL
% Version: 220921 FIXED error in selection of items

% The Trajectories field may be absent (e.g. Theia markerless)
if ~isfield(Data, 'Trajectories')
    return;
end

% Get the list of all parameters and select the indices for the labels and the calculated parameters
Parameters       = {ParameterInfo.parameter};
PointLabelsId    = PointGroup & matches(Parameters,'LABELS');
HasFilltypeField = isfield(Data.Trajectories.Labeled,'Filltype');
LabelList        = ParameterInfo(PointLabelsId).data;

% Copy the Labels, Data and Filltype into the Calculated substructure
DataNames     = {'Angles','Forces','Moments','Powers','ModeledMarkers'};
ParamNames    = {'ANGLES','FORCES','MOMENTS','POWERS','MODELED_MARKERS'};
UnitNames     = {'ANGLE_UNITS','FORCE_UNITS','MOMENT_UNITS','POWER_UNITS','MODELED_MARKER_UNITS'};

% Initialize all labels as "not to be moved"
MovedItems = false(1,length(LabelList)); 
% Loop the kinds of calculated parameters
for Nm=1:length(DataNames)
    ParamId = PointGroup & matches(Parameters,ParamNames{Nm});
    if any(ParamId)
        % Select the indices of the calculated parameter from the LABELS list
        CurrentItems = matches(LabelList,ParameterInfo(ParamId).data); % logical array
        MovedItems = MovedItems | CurrentItems; % logical array
        % Copy the data, labels and units to the "Calculated" substructure
        Data.Calculated.(DataNames{Nm}).Data   = Data.Trajectories.Labeled.Data(CurrentItems,:,:);
        Data.Calculated.(DataNames{Nm}).Labels = LabelList(CurrentItems);
        Data.Calculated.(DataNames{Nm}).Units  = ParameterInfo(PointGroup & matches(Parameters,UnitNames{Nm})).data{:};
        % Copy the Filltype field if present
        if HasFilltypeField
            Data.Calculated.(DataNames{Nm}).Filltype = Data.Trajectories.Labeled.Filltype(CurrentItems,:);
        end
    end
end

% Delete the copied items from the Trajectories substructure
Data.Trajectories.Labeled.Labels(MovedItems) = [];
Data.Trajectories.Labeled.Data(MovedItems,:,:) = [];
if HasFilltypeField
    Data.Trajectories.Labeled.Filltype(MovedItems,:) = [];
end
end

%% =================================================================================================

function Data = defineForce(Data,ParameterInfo,ForceGroup,AnalogGroup)
%% Define the Force substructure, if present, and write them into the Data structure
% The force plate data in C3D are imported with the analog data and thus represent channels of the
% Analog substructure. 
% Combine the channels into 3D Forces and COPs
% 
% INPUT
%    Data           - (struct) Input Data according to QTM matlab structure
%    ParameterInfo  - (struct) containing grouped information about all parameters (i.e. the header)
%    ForceGroup     - (logical 1xn) indices of ParameterInfo that belong to the group
%    AnalogGroup    - (logical 1xn) indices of ParameterInfo that belong to the group
%
% OUTPUT
%    Data         - (struct) see above
% 
% local functions: getCOP, getTara
% 
% (c) 2022 by Predimo GmbH
% 230104 (MdL) correct sign of the force (ForcePlateOrientation)
% 230308 (MdL) deal with empty BoardName field (theia); set force plate location by Units
% 230321 (MdL) removed scale force plate location by units (after testing with PUMA data)

Parameters     = {ParameterInfo.parameter};
IsDefinedForce = any(ForceGroup) && ParameterInfo(ForceGroup & matches(Parameters,'USED')).data>0;
if IsDefinedForce
    % general parameters
    MarkerRate    = Data.FrameRate;
    NMarkerFrames = Data.Frames;
    AnalogRate    = ParameterInfo(AnalogGroup & matches(Parameters,'RATE')).data;
    NAnalogFrPSm  = AnalogRate/MarkerRate;

    % specific parameters
    NForcePlates  = ParameterInfo(ForceGroup & matches(Parameters,'USED')).data;
    ForceType     = ParameterInfo(ForceGroup & matches(Parameters,'TYPE')).data;
    ForceCorners  = ParameterInfo(ForceGroup & matches(Parameters,'CORNERS')).data;
    ForceOrigin   = ParameterInfo(ForceGroup & matches(Parameters,'ORIGIN')).data;
    ForceChannels = ParameterInfo(ForceGroup & matches(Parameters,'CHANNEL')).data;
    ForceZero     = ParameterInfo(ForceGroup & matches(Parameters,'ZERO')).data;
    if any(matches(Parameters,'CAL_MATRIX'))
        Cal_Matrix= ParameterInfo(ForceGroup & matches(Parameters,'CAL_MATRIX')).data;
    end
    
    % parameters for each force plate
    for fp = 1 : NForcePlates
        PlateChannels = squeeze(ForceChannels(:,fp)); % the first analog channel for this plate
        BoardName     = Data.Analog(fp).BoardName;
        % Make sure that it is clear that we are dealing with force plates
        if isempty(BoardName) || ~contains(BoardName,{'force plate','force-plate'},'ignorecase',true)
            BoardName = sprintf('Force-plate %s', BoardName);
            for Channel = 1:length(PlateChannels)
                Data.Analog(Channel).BoardName = BoardName;
            end
        end
        % These are the standard fields of the QTM structure
        Data.Force(fp).ForcePlateName = BoardName;
        Data.Force(fp).NrOfFrames     = NMarkerFrames;
        Data.Force(fp).SamplingFactor = NAnalogFrPSm;
        Data.Force(fp).NrOfSamples    = NMarkerFrames*NAnalogFrPSm;
        Data.Force(fp).Frequency      = MarkerRate * NAnalogFrPSm;
        Location                      = squeeze(ForceCorners(:,:,fp)) + squeeze(ForceOrigin(:,fp));
        % check that the corners are in the rows
        if all(size(Location) == [3 4])
            Location = Location'; % transpose
        end
        Data.Force(fp).ForcePlateLocation = Location;
        % - Sign of the force: the order of the force plate corners defines the sign of the force
        % If the order is clockwise, the force is down, if anti-clockwise, it is up (see manual).
        % - According to Qualisys: when the flag is 1, we have Ground Action forces, if 0, we have 
        % Ground Reaction forces. In the wording by QTM: "Coordinate system in which force data 
        % is expressed: 0 (local force plate coordinates), 1 (global coordinate system)"
        XAxis = Location(1,:)-Location(2,:); % corner 1 - corner 2
        YAxis = Location(2,:)-Location(3,:);
        ZVector = cross(XAxis, YAxis);
        Data.Force(fp).ForcePlateOrientation = ZVector(3) < 0;
        if any(matches(Parameters,'CAL_MATRIX'))
            Data.Force(fp).CalibrationMatrix = squeeze(Cal_Matrix(:,:,fp));
        end
        % The following are the default fields of c3d and are thus redunant; copy them nevertheless
        % for completeness
        Data.Force(fp).ForceChannel   = PlateChannels;
        Data.Force(fp).ChannelLabels  = Data.Analog(PlateChannels).Labels;
        Data.Force(fp).ChannelUnits   = Data.Analog(PlateChannels).Units;
        Data.Force(fp).ForceCorners   = squeeze(ForceCorners(:,:,fp));
        Data.Force(fp).ForceOrigin    = squeeze(ForceOrigin(:,fp));
        Data.Force(fp).ForceType      = ForceType(fp); % 1,2,3,or 4 depending on manufacturer, see c3d manual
        Data.Force(fp).ForceZero      = ForceZero; % this is the range of samples over which to extract the offset
    end
else
    NForcePlates = 0;
    ForceType    = 0;
end

%% Fill the force data into the Force substructure
for fp = 1 : NForcePlates
    Channels = Data.Force(fp).ForceChannel;
    % Add the force data according to the type of force plate, as defined by ForceType
    switch ForceType(fp)
        case 1 % Channel 1-3 : 3d Force; Channels 4-5 : COP in x and y; Channel 6 Moment z
            if length(Channels) ~=6
                warning('The ForceType (1) is not consistent with the c3d definition (N channels (%d) should be 6)',length(Channels));
            else
                Force  = [Data.Analog(Channels(1:3)).Data]';
                Force  = Force - getTara(Force,fp,Data); % subtract the offset if needed (Tara or ForceZero)
                COP    = [Data.Analog(Channels(4:5)).Data]';
                Moment = [Data.Analog(Channels(6)).Data]';
                Data.Force(fp).Force  = Force;
                Data.Force(fp).Moment = Moment;
                Data.Force(fp).COP    = COP;
            end
        case 2 % Channel 1-3 : 3d Force; Channels 4-6: 3d Moment  (e.g. AMTI and Qualisys QTM)
            if length(Channels) ~=6
                warning('The ForceType (2) is not consistent with the c3d definition (N channels (%d) should be 6)',length(Channels));
            else
                Force  = [Data.Analog(Channels(1:3)).Data]';
                Force  = Force - getTara(Force,fp,Data); % subtract the offset if needed (Tara or ForceZero)
                Moment = [Data.Analog(Channels(4:6)).Data]';
                Data.Force(fp).Force = Force;
                Data.Force(fp).Moment = Moment;
                Data.Force(fp).COP = getCOP(Force,Moment,squeeze(ForceCorners(:,:,fp)),ForceOrigin(:,fp)');
            end
        case 3 % Cannel 1-2: x components; Cannel 3-4: y components; Cannel 5-8: z components
            if length(Channels) ~=8
                warning('The ForceType (3) is not consistent with the c3d definition (N channels (%d) should be 6)',length(Channels));
            else
                Force(1,:)  = sum([Data.Analog(Channels(1:2)).Data],2);
                Force(2,:)  = sum([Data.Analog(Channels(3:4)).Data],2);
                Force(3,:)  = sum([Data.Analog(Channels(4:8)).Data],2);
                Force       = Force - getTara(Force,fp,Data); % subtract the offset if needed (Tara or ForceZero)
                Data.Force(fp).Force = Force;
                warning('Moment and COP are not calculated! ADD THIS FEATURE IN FUTURE RELEASE');
            end
        case 4 % as Type 2, but with calibration matrix
            if length(Channels) ~=6
                warning('The ForceType (4) is not consistent with the c3d definition (N channels (%d) should be 6)',length(Channels));
            else
                Force  = [Data.Analog(Channels(1:3)).Data]';
                Force  = Force - getTara(Force,fp,Data); % subtract the offset if needed (Tara or ForceZero)
                Moment = [Data.Analog(Channels(4:6)).Data]';
                Data.Force(fp).Force = Force;
                Data.Force(fp).Moment = Moment;
                Data.Force(fp).COP = getCOP(Force,Moment,squeeze(ForceCorners(:,:,fp)),...
                    ForceOrigin(:,fp)',Data.Force(fp).CalibrationMatrix);
            end
        otherwise
            warning('This type of ForceType definition is not implemented; see analog substructure')
    end
end
end

%% =================================================================================================

function Tara = getTara(Force,Plate,Data)
%% Get the forceplate offsets to be subtracted (Tara or ForceZero)
% 
% INPUT
%    Force  - (3 x nsamples) measured force of current plate [N]
%    Plate  - (double) force plate index 
%    Data   - (struct) Input Data structure
%
% OUTPUT
%    Tara   - (3x1 double) tara value of the current force plate [N]
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: MdL
% Last revision: 13-Mar-2022 improved warning

% The tara-range ("ForceZero") is the range of samples over which the offset is to be calculated
TaraRange = Data.Force(Plate).ForceZero; 
if ~all(TaraRange==0) && TaraRange(1) < TaraRange(2)
    % error checks: indices must be positive int and smaller than the array length
    if TaraRange(1)<1, TaraRange(1)=1; end
    TaraRange(2) = min(TaraRange(2),size(Force,2));
    Tara = mean(Force(:,TaraRange(1) : TaraRange(2)),2,'omitnan');
else % default: do not calculate a tara offset
    Tara = [0; 0; 0];
end
end

%% =================================================================================================

function COP = getCOP(Force,Moment,ForceCorners,ForceOrigin,CalibrationMatrix)
%% Calculate the center of pressure from the force and the moment.
% - We do no know if we are dealing with action or reaction force, so detect that automatically.
% - The X and Y direction of the force plates is determined from the ForceCorners.
% - The force plate cordinate system may be rotated with respect to the global system by multiples
% of 90 degrees: this is handled automatically.
% 
% INPUT
%    Force        - (double 3 x nsamples) measured force of current plate [N]
%    Moment       - (double 3 x nsamples) measured force moment of current plate [Nm]
%    ForceCorners - (double 3 x 4) locations of the corners of the plate
%    ForceOrigin  - (double 3 x 1) origin of the plate coordinate system
%    CalibrationMatrix - (double) calibration matrix of the force plate
%
% OUTPUT
%    COP - (3 x nsamples double) center of pressure
% 
% See also: getCOPfromAnalog, defineForce
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: MdL
% version 230227 (MdL) detect sign of forces automatically; FIXED 90-deg rotation of the plates
% version 230304 (MdL) detect sign of X and Y, and thus COP direction for all cases; revised header
% and comments
% version 230308 (MdL) define IsAligned parameter also if the plate is oblique

%% Init
narginchk(4,5);
if nargin>4
    % if the calibration matrix is provided and unequal from the identity matrix, issue a warning
    if any(abs(CalibrationMatrix - eye(6)) > 1e-6, 'all')
        warning('The calibration matrix is not yet implemented');
    end
end
X=1;Y=2;Z=3;
AlignThreshold = 5; % mm
LoadThreshold  = 50; % N
DoPlot         = false;

% The radius is the half size in x and y.
% Assume that the axes of the force plates align with the global coordinates.
% Per default, corners are numbered anti-clockwise, with minimum X and Y in corner 1

% In case the plate coords are aligned (+/- 180 deg) with the global coord:
AlignmentXCorners = ForceCorners(X,4) - ForceCorners(X,1); % corner 1 and 4 should have same x value
AlignmentYCorners = ForceCorners(Y,2) - ForceCorners(Y,1); % corner 1 and 2 should have same y value
% ... alternatively they are rotated by +/- 90 degrees
AlignmentXCornersRotated = ForceCorners(X,2) - ForceCorners(X,1); % rotated : corner 1 and 2 should have same x value
AlignmentYCornersRotated = ForceCorners(Y,4) - ForceCorners(Y,1); % rotated : corner 1 and 4 should have same y value
PlateCenter = mean(ForceCorners,2);

% Is the plate aligned with the global coordinates (or 90-deg rotated)?
isAlignedWithGlobalCoordinates  = abs(AlignmentXCorners)  < AlignThreshold && abs(AlignmentYCorners)  < AlignThreshold;
isAlignedWithRotatedCoordinates = abs(AlignmentXCornersRotated) < AlignThreshold && abs(AlignmentYCornersRotated) < AlignThreshold;
% check if the plates are rotated by 90 degrees
if      isAlignedWithGlobalCoordinates &&  ~isAlignedWithRotatedCoordinates
    IsAligned = true;
elseif ~isAlignedWithGlobalCoordinates &&   isAlignedWithRotatedCoordinates
    IsAligned = false;
else
    warning('loadc3d: The force plate might have an oblique orientation; this is currently not implemented (Fehler x,y: %.1f, %.1fmm)', ...
        AlignmentXCorners,AlignmentYCorners);
    IsAligned = true;
end

% Get the loaded phases: The COP is only defined as long as there is a load, 
% so detect when the plate is loaded
LoadZ = Force(Z,:);
% Detect the sign of the force (action or reaction force) automatically
IsMeanNegative = mean(LoadZ,'omitnan') < 0;
IsMaxNegative = max(LoadZ,[],'omitnan') < -min(LoadZ,[],'omitnan');
if IsMeanNegative && IsMaxNegative
    Loaded = -LoadZ > LoadThreshold;
else
    Loaded =  LoadZ > LoadThreshold;
end

% Compute the COP
COP = zeros(size(Force));
COP(X:Y,:) = Moment([Y X],:) ./ Force(Z,:); 
COP(:,~Loaded) = nan;

% If the plate is rotated by 180 deg, or if the corners are numbered clockwise (cf. c3d manual),
% then the X, Y, or both axes must be swapped
SignX = 1;
SignY = 1;
if IsAligned
    % Swap x and y is not needed
else
    % Rotate 90 degrees anticlockwise
    COP = COP([Y X Z],:);
    COP(X,:) = -COP(X,:);
    % Swap x and y, if necessary 
    if AlignmentXCorners < 0, SignX = -1; end
    if AlignmentYCorners < 0, SignY = -1; end
end
COP(X,:) = SignX * COP(X,:);
COP(Y,:) = SignY * COP(Y,:);

% Correct the up-down orientation
Norm = cross(ForceCorners(:,3)-ForceCorners(:,4),ForceCorners(:,1)-ForceCorners(:,4));
if Norm(Z) < 0
    COP(Z,:) = -COP(Z,:);
end
% Correct for oblique force direction
ObliqueX = - SignX * abs(ForceOrigin(Z)) * Force(X,:) ./ Force(Z,:);
ObliqueY = - SignY * abs(ForceOrigin(Z)) * Force(Y,:) ./ Force(Z,:);
ObliqueX(~Loaded) = 0;
ObliqueY(~Loaded) = 0;

% Add the force plate location
COP(X,:) = COP(X,:) + ObliqueX + ForceOrigin(X);
COP(Y,:) = COP(Y,:) + ObliqueY + ForceOrigin(Y);
COP      = COP      + PlateCenter;

% Plot controll figure if desired
% if abs(ForceOrigin(Z))>1, DoPlot = true; end
if DoPlot
    Start = find(~isnan(COP(X,:)),10); %#ok<UNRCH>
    figure; hold on; title('COP trajectory on plate') 
    plot(ForceCorners(X,[1 2 3 4 1]),ForceCorners(Y,[1 2 3 4 1]),'linewidth',2); axis('equal')
    plot(COP(X,:),COP(Y,:));     plot(COP(X,Start),COP(Y,Start),'k','linewidth',2);
    if abs(ForceOrigin(Z))>1
        plot(COP(X,:)-ObliqueX, COP(Y,:)-ObliqueY);
        legend('Plate','COP','start','No-oblique');
        figure; plot(Force'); title('force in 3-d')
    end
    fprintf('%d_%d Ax Ay Rx Ry: %.0f %.0f; %.0f %.0f\n',IsAligned,IsMeanNegative && IsMaxNegative, .1*AlignmentXCorners,.1*AlignmentYCorners,.1*AlignmentXCornersRotated,.1*AlignmentYCornersRotated);
end
end
