function [Data,SimInMN,Flag,Msgs] = loadInputFile(FPath,Flag,Msgs)
%% load the Data file and extract header data
% Data can be c3d or QTM-mat
% 
% SYNTAX
%     [Data,SimInMN,Flag,Msgs] = loadInputFile(FPath,Flag,Msgs)
%
% INPUT
%     FPath        (string)  qualisys data file
%     Flag         (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs         (char)   notifications, warning and error messages
%
% OUTPUT
%     Data         (QTM struct) all available information from input file in Qualisys format
%     SimInMN      (MNData struct) all available information from input file
%     Flag         (double) Exit status of the function (1 = success, 0 = warning, -1 = error)
%     Msgs         (char)   notifications, warning and error messages
%
% EXAMPLES
%     [Data,SimInMN,Flag,Msgs] = loadInputFile(FPath,Flag,Msgs)
% 
% Local functions: omitEmptyTheiaEnd
%
% see also: importBodyFromQUAL
% 
% (c) 2023 Predimo GmbH
% Website: http://www.predimo.com
% authors: Marc de Lussanet
% version 230208 (MdL) omitEmptyTheiaEnd : check for force field
% version 230221 (MdL) check for file length
% version 230306 (MdL) check if file exists
% version 230310 (MdL) check units
% version 230327 (MdL) Units for Qualisys QTM is mm
% version 230329 (MvdH) Corrected bug

% Init
SimInMN = struct;  
if ~isfile(FPath)
    Flag = -1;
    Msgs = [Msgs, {sprintf('File %s not readable',FPath)}];
end

% Detect the filetype
[~,Fname,Ext] = fileparts(FPath);
fprintf('Loading "%s%s" ...\n',Fname,Ext);
if matches(Ext,'.c3d','IgnoreCase',true)
    Data = loadc3d(FPath); % [Data, InfoC3D]
    fprintf('... finished loading c3d file\n\n')
elseif matches(Ext,'.mat','IgnoreCase',true)
    tmp    = load(FPath);
    fnames = fieldnames(tmp);
    Data   = tmp.(fnames{1});
    if isempty(Data) || ~isstruct(Data)
        Flag = -1; % error code
        Msgs = [Msgs, {sprintf('Unexpected datastructure in file "%s" : Is this QTM data? -> Aborting import',FPath)}];
        return;
    end
    if ~isfield(Data,'Units')
        % Qualisys mat has default mm, but no Units field
        Data.Units = 'mm';
        fprintf('... finished loading QTM file\n\n')
    end
else
    Flag = 0;
    Msgs = [Msgs, {sprintf('loadInputFile: file "%s%s" not recognized',Fname,Ext)}];
end
if isempty(fieldnames(Data)) || ~isfield(Data,'File')
    Flag = -1;
    Msgs = [Msgs, {sprintf('loadInputFile: Could not open input file')}];
    return;
end

% convert QTM unit [mm] to SI [m]
if isfield(Data,'Units')
    if strcmp(Data.Units,'mm')
        Scale = 0.001;
    elseif strcmp(Data.Units,'m')
        Scale = 1;
    else
        Flag = -1;
        Msgs = [Msgs, {sprintf('loadInputFile: non-defined units %s',Data.Units)}];
        return;
    end
else % if there is no Units field at this point, it is an error. We throw a warning and a fix
    Data.Units = 'm';
    Scale = 1; %#ok<NASGU> 
    Flag = 0;
    Msgs = [Msgs, {sprintf('Units field is missing in file "%s" -> Assuming [m] for lengths',FPath)}];
    return;
end
if isfield(Data,'Trajectories') && isfield(Data.Trajectories,'Labeled') && ...
        isfield(Data.Trajectories.Labeled,'Data')
    Data.Trajectories.Labeled.Data(:,1:3,:) = Data.Trajectories.Labeled.Data(:,1:3,:) * Scale;
end

% build up the structure of SimInMN
if ~isfield(SimInMN,'meta')
    SimInMN.meta = struct;
    if ~isfield(SimInMN.meta,'static')
    SimInMN.meta.static = struct;
    end
end
SimInMN.meta.static.measurementName = FPath;
% copy header info from the data
if ~isfield(Data,'StartFrame')
    Data.StartFrame = 1;
end

% detect and omit empty samples at end of Theia file
if isfield(Data,'Segments') % Theia Data
    Data = omitEmptyTheiaEnd(Data);
end
SimInMN.meta.static.fileNameOrig   = Data.File;
SimInMN.meta.static.nFramesOrig    = Data.Frames;
SimInMN.meta.static.startFrameOrig = Data.StartFrame;
if Data.Frames < 4
    Flag = -1;
    Msgs = [Msgs, {sprintf('loadInputFile: File is too short (%d<4 samples)',Data.Frames)}];
end
end

%% ========================================================================

function Data = omitEmptyTheiaEnd(Data)
%% We noticed, that Theia data sometimes have empty samples at the end of the file. 
% Since extrapolation of data is generally a bad idea, these empty samples are useless. 
% Moreover they may lead to unexpected results of smoothing etc. 
% 
% Note: At present, we only remove such samples for Theia data, but should they be present in other
% formats, these should be truncated also. Given the complexity of the data, each case should be
% thoroughly thought-through and tested. (hence at present only Theia)
% 
% 
% SYNTAX
%     Data = omitEmptyTheiaEnd(Data)
%
% INPUT
%     Data         (QTM struct) all available information from input file in Qualisys format
%
% OUTPUT
%     Data         (QTM struct) all available information from input file in Qualisys format
% 
% see also: loadInputFile
% 
% (c) 2023 Predimo GmbH
% Website: http://www.predimo.com
% authors: Marc de Lussanet
% version 230208 (MdL) check for force field

% the segments that do have any non-zero values
NonZeroSegments = ~all(Data.Segments.Labeled.Pos==0,[2 3]);
% Samples that have only nan values
NanSamples = all(isnan(Data.Segments.Labeled.Pos(NonZeroSegments,:,:)),[1 3]);
NSamples = length(NanSamples);
EmptyAtEnd = NSamples - find(~NanSamples,1,'last');
% Force-samples that have only zero values
if isfield(Data,'Force')
    NSamplesF = Data.Force(1).NrOfSamples;
    NPlatesF = length(Data.Force);
    Forces = reshape([Data.Force.Force],[3 NSamplesF NPlatesF]);
    ZerosSamplesForce = all(Forces==0,[1 3]);
    EmptyEndForce = NSamplesF - find(~ZerosSamplesForce,1,'last');
    EmptyEndSamplesByForce = EmptyEndForce / Data.Force(1).SamplingFactor;
else
    NPlatesF  = 0;
    EmptyEndSamplesByForce = 0;
end
NEmptySamples = max([EmptyAtEnd EmptyEndSamplesByForce]);
if NEmptySamples > NSamples / 2
    warning('loadInputFile : apparently more than half the samples are empty');
else
    NSamples = NSamples - NEmptySamples;
    Data.Frames = NSamples;
    Data.Segments.Labeled.Pos(:,:, NSamples+1 : end) = [];
    if isfield(Data,'Trajectories')
        Data.Trajectories.Labeled.Data(:,:, NSamples+1 : end) = [];
    end
    for i=1:NPlatesF
        NSamplesF = NSamples * Data.Force(1).SamplingFactor;
        Data.Force(i).NrOfFrames = NSamples;
        Data.Force(i).NrOfSamples = NSamplesF;
        Data.Force(i).Force(:, NSamplesF+1:end, :) = [];
        Data.Force(i).Moment(:, NSamplesF+1:end, :) = [];
        Data.Force(i).COP(:, NSamplesF+1:end, :) = [];
    end
    fprintf('Omitted %d empty samples at the end of the file\n',NEmptySamples);
end
end
