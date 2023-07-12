function makeMotorData()

fprintf('\nMaking Motor Data...\n');

motorDir = evalin('base', 'motorDir');

% load current state
Measurements = loadState;

% load table containing subjects info
Subjects = Measurements.Subjects;

% list content of input folder into cell array
dirInfo = dir(fullfile(motorDir, '*.mat'));
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxExclude = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fnames(idxExclude) = [];
fdirs(idxExclude) = [];
fpaths = strcat(fdirs, filesep, fnames);

gapMargin = 0.1; % margin around gaps in seconds, because there the COP is distorted
LoadThresh = 20; % threshold below which forces are set to zero [N]
InitForce = 200; % force threshold to cross to start trial [N]
InitMarg = 0.5; % additional time after crossing InitForce to start trial [s]
HumFreq = 50; % humming frequency in Hz

fprintf('Extract measurement data...\n');
ticAll = tic;

if isfield(Measurements, 'MotorData') && isfield(Measurements.MotorData, 'fileName')
    loadedFiles = {Measurements.MotorData.fileName}';
else
    loadedFiles = {};
end
nFiles = length(fpaths);
fileNr = length(loadedFiles) + 1; % init file index to end of already loaded files
for iFile = 1:nFiles
    fpath = fpaths{iFile};
    [~, fname, ~] = fileparts(fpath);

    if any(strcmp(fname, loadedFiles))
        continue
    end
    tic

    % split file name into parts separated by underscore
    parts = strsplit(fname, '_');

    % subject identity and code
    subject = parts{1};
    subjectCode = [subject, '_', parts{2}];
    subjectCodes = [Subjects.Code_I, Subjects.Code_II, Subjects.Code_III];
    [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
    if isempty(subjectIdx)
        continue
    end

    % jump position markers
    Beidbein_start = string(Subjects.Beidbein_start(subjectIdx));
    Beidbein_stop = string(Subjects.Beidbein_stop(subjectIdx));
    Einbein_start = string(Subjects.Einbein_start(subjectIdx));
    Einbein_stop = string(Subjects.Einbein_stop(subjectIdx));

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        side = parts{5};
    else
        side = parts{4};
    end
   
    %% Import Force data

    tmp = load(fpath);
    fields = fieldnames(tmp);
    MotorData = tmp.(fields{1});
    lengthScale = 0.001; % mm -> m
    [Forces, sampleRate, COPs, ~, location, Analog] = kraftAusQualisys(MotorData, lengthScale);
    PowerLineHum = 50;
    nSamples = length(Forces);
    if nSamples < 2^nextpow2(PowerLineHum)/2
        fprintf('\t\tNot enough data (%d samples) -> skipping\n', nSamples);
        continue
    end
    Forces = -Forces; % ground reaction forces are inverse of plate forces
    nPlates = size(Forces, 1);
    dt = 1/sampleRate; % time step size [s]
    HumPeriod = sampleRate/HumFreq;
    for iPlate = 1:nPlates
        Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
    end
    [COPs, ~, Forces, ~] = getCOPfromAnalog(Analog, Forces, [], location, sampleRate, [], [], [], COPs);
    % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
    for iPlate = 1:nPlates
        idxPlateGaps = squeeze(abs(Forces(iPlate, 3, :)) < LoadThresh)';
        Forces(iPlate, :, idxPlateGaps) = 0;
        COPs(iPlate, :, idxPlateGaps) = NaN;
    end
    Force = squeeze(sum(Forces, 1, 'omitnan'));
    absForces = abs(Forces(:, 3, :)); % z-component of force
    weightedCOPs = (absForces .* COPs) ./ sum(absForces, 1);
    COPxyz = squeeze(sum(weightedCOPs, 1, 'omitnan'));
    COP = COPxyz(1:2, :);

    iStart = find(abs(Force(3,:))>InitForce, 1, 'first');
    iStart = iStart + floor(InitMarg/dt);
    Force = Force(:, iStart:end);
    iStop = find(abs(Force(3,:))>InitForce, 1, 'last');
    iStop = iStop - floor(InitMarg/dt);
    Force = Force(:, 1:iStop);
    COP = COP(:, iStart:iStop+iStart-1);
    nSamples = size(Force, 2);
    if nSamples == 0
        fprintf('\t\tToo few samples with nonzero force -> skipping\n');
        continue
    end
    Time = (0:nSamples-1)/sampleRate;

    %% Remove non-contact phases from COP

    idxGaps = ~any(COP, 1);
    COP(:, idxGaps) = NaN;

    %% Add gap margins

    gapMarginSmp = round(gapMargin * sampleRate); % gap margin in samples
    % find the beginnings of each gap
    gapStart = [-1, find(idxGaps)];
    dGapStart = [1, diff(gapStart)];
    gapStart(dGapStart == 1) = [];
    % find the endings of the gaps
    gapStop = [find(idxGaps), nSamples+2];
    dGapStop = [diff(gapStop), 1];
    gapStop(dGapStop == 1) = [];
    % enlarge the gaps by margin
    gapStart = gapStart - gapMarginSmp;
    gapStop = gapStop + gapMarginSmp;
    % keep within data range
    gapStart(gapStart < 1)  = 1;
    gapStop(gapStop > nSamples) = nSamples;
    % re-create bool array for non-gap indices
    for iSample = 1:length(gapStart)
        idxGaps(gapStart(iSample):gapStop(iSample)) = true;
    end

    %% set COP within gaps to NaN

    COP(:, idxGaps) = NaN;

    %% Import marker data
    % get jump start and stop positions

    startPos = nan(3, 1);
    stopPos = nan(3, 1);
    if isfield(MotorData, 'Trajectories') && isfield(MotorData.Trajectories, 'Labeled') && isfield(MotorData.Trajectories.Labeled, 'Labels')
        labels = MotorData.Trajectories.Labeled.Labels;
        % get start and stop marker labels
        switch side
            case 'B'
                startLabel = Beidbein_start;
                stopLabel = Beidbein_stop;
            case {'L', 'R'}
                startLabel = Einbein_start;
                stopLabel = Einbein_stop;
            otherwise
        end
        % get start and stop marker data
        idx = strcmp(labels, startLabel);
        if any(idx)
            pos_raw = squeeze(MotorData.Trajectories.Labeled.Data(idx, 1:3, :))/1000;
            startPos = mean(pos_raw, 2);
        end
        idx = strcmp(labels, stopLabel);
        if any(idx)
            pos_raw = squeeze(MotorData.Trajectories.Labeled.Data(idx, 1:3, :))/1000;
            stopPos = mean(pos_raw, 2);
        end
    end

    %% Store data

    Measurements.MotorData(fileNr).fileName = fname;
    Measurements.MotorData(fileNr).nSamples = nSamples;
    Measurements.MotorData(fileNr).sampleRate = sampleRate;
    Measurements.MotorData(fileNr).Time = Time;
    Measurements.MotorData(fileNr).Force = Force;
    Measurements.MotorData(fileNr).COP = COP;
    Measurements.MotorData(fileNr).idxContact = ~idxGaps;
    Measurements.MotorData(fileNr).startPos = startPos;
    Measurements.MotorData(fileNr).stopPos = stopPos;    

    % report progress
    fprintf('\t-> %s (%d/%d = %.1f%% in %.3fs)\n', fname, iFile, nFiles, iFile/nFiles*100, toc);

    % export Measurements structure to base workspace
    assignin('base', 'Measurements', Measurements);

    % increment number of processed files
    fileNr = fileNr+1;
end

% save current state
saveState;

fprintf('Finished extracting data from %d files in %.3f s\n', fileNr - (length(loadedFiles) + 1), toc(ticAll));

end