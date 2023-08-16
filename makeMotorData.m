function makeMotorData()

fprintf('\nMaking Motor Data...\n');

motorDir = evalin('base', 'motorDir');

% load current state
Measurements = loadState();

% load table containing subjects info
SubjectsTable = Measurements.Subjects;

% dirInfo = dir(motorDir);
% motorDirs = {dirInfo.name}';
% idxExclude = startsWith(motorDirs', {'.', '~'}) | ~[dirInfo.isdir];
% motorDirs(idxExclude) = [];
% motorDirs = {'PR'};
motorDirs = cellstr(motorDir);

for iDir = 1:length(motorDirs)

    myMotorDir = motorDirs{iDir};
    % list content of input folder into cell array
    dirInfo = dir(fullfile(myMotorDir, '*.mat'));
    fdirs = {dirInfo.folder}';
    fnames = {dirInfo.name}';
    idxExclude = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
    fnames(idxExclude) = [];
    fdirs(idxExclude) = [];

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

    [~, fnames, fexts] = fileparts(fnames);
    [~, iA] = setdiff(fnames, loadedFiles);
    fpaths = strcat(fdirs(iA), filesep, fnames(iA), fexts(iA));

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
        subjectCodes = [SubjectsTable.Code_I, SubjectsTable.Code_II, SubjectsTable.Code_III];
        [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
        if isempty(subjectIdx)
            % fprintf('\t%s:\tNo subject with code ''%s'' -> skipping\n', fname, subjectCode);
            continue
        end

        % jump position markers
        Beidbein_start = string(SubjectsTable.Beidbein_start(subjectIdx));
        Beidbein_stop = string(SubjectsTable.Beidbein_stop(subjectIdx));
        Einbein_start = string(SubjectsTable.Einbein_start(subjectIdx));
        Einbein_stop = string(SubjectsTable.Einbein_stop(subjectIdx));

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
        try
            [Forces, sampleRateForce, COPs, ~, location, Analog] = kraftAusQualisys(MotorData, lengthScale);
        catch ME
            fprintf('\t%s:\tError when extracting forces -> skipping\n', fname);
            fprintf('\t\t%s\n', ME.message);
            continue
        end
        PowerLineHum = 50;
        nSamplesForce = length(Forces);
        if nSamplesForce < 2^nextpow2(PowerLineHum)/2
            fprintf('\t%s:\tNot enough data (%d samples) -> skipping\n', fname, nSamplesForce);
            continue
        end
        Forces = -Forces; % ground reaction forces are inverse of plate forces
        nPlates = size(Forces, 1);
        dtForce = 1/sampleRateForce; % time step size [s]
        HumPeriod = sampleRateForce/HumFreq;
        for iPlate = 1:nPlates
            Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
        end
        [COPs, ~, Forces, ~] = getCOPfromAnalog(Analog, Forces, [], location, sampleRateForce, [], [], [], COPs);
        % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
        for iPlate = 1:nPlates
            idxPlateGaps = squeeze(abs(Forces(iPlate, 3, :)) < LoadThresh)';
            if all(idxPlateGaps)
                Forces(iPlate, :, :) = 0;
            end
            COPs(iPlate, :, idxPlateGaps) = NaN;
        end
        Force = squeeze(sum(Forces, 1, 'omitnan'));        
        absForces = abs(Forces(:, 3, :)); % z-component of force
        weightedCOPs = (absForces .* COPs) ./ sum(absForces, 1);
        COPxyz = squeeze(sum(weightedCOPs, 1, 'omitnan'));
        COP = COPxyz(1:2, :);       

        iStart = find(abs(Force(3,:))>InitForce, 1, 'first');
        iStart = iStart + floor(InitMarg/dtForce);
        Force = Force(:, iStart:end);
        iStop = find(abs(Force(3,:))>InitForce, 1, 'last');
        iStop = iStop - floor(InitMarg/dtForce);
        Force = Force(:, 1:iStop);
        idxForce =  iStart:iStop+iStart-1;
        COP = COP(:, idxForce);
        nSamplesForce = size(Force, 2);
        if nSamplesForce == 0
            fprintf('\t%s:\tToo few samples with nonzero force -> skipping\n', fname);
            continue
        end
        TimeForce = (0:nSamplesForce-1)/sampleRateForce;
        sampleRateKin = MotorData.FrameRate;        
        iStartKin = ceil(iStart/sampleRateForce*sampleRateKin);
        iStopKin = ceil(iStop/sampleRateForce*sampleRateKin);
        idxKin =  iStartKin:iStopKin+iStartKin-1;   
        nSamplesKin = length(idxKin);
        TimeKin = (0:nSamplesKin-1)/sampleRateKin;

        %% Remove non-contact phases from COP

        idxGaps = ~any(COP, 1);
        COP(:, idxGaps) = NaN;

        %% Add gap margins

        gapMarginSmp = round(gapMargin * sampleRateForce); % gap margin in samples
        % find the beginnings of each gap
        gapStart = [-1, find(idxGaps)];
        dGapStart = [1, diff(gapStart)];
        gapStart(dGapStart == 1) = [];
        % find the endings of the gaps
        gapStop = [find(idxGaps), nSamplesForce+2];
        dGapStop = [diff(gapStop), 1];
        gapStop(dGapStop == 1) = [];
        % enlarge the gaps by margin
        gapStart = gapStart - gapMarginSmp;
        gapStop = gapStop + gapMarginSmp;
        % keep within data range
        gapStart(gapStart < 1)  = 1;
        gapStop(gapStop > nSamplesForce) = nSamplesForce;
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
                case 'L'
                    startLabel = Einbein_start;
                    stopLabel = Einbein_stop;
                    footLabels = {'LFM1'};
                case 'R'
                    startLabel = Einbein_start;
                    stopLabel = Einbein_stop;
                    footLabels = {'RFM1'};
                case 'B'
                    startLabel = Beidbein_start;
                    stopLabel = Beidbein_stop;
                    footLabels = {'LFM1', 'RFM1'};
                otherwise
            end
            % get start and stop marker data
            idx = strcmp(labels, startLabel);
            if any(idx)
                pos_raw = squeeze(MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin))/1000;
                startPos = mean(pos_raw, 2);
            end
            idx = strcmp(labels, stopLabel);
            if any(idx)
                pos_raw = squeeze(MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin))/1000;
                stopPos = mean(pos_raw, 2);
            end
            % get foot marker data
            idx = ismember(labels, footLabels);
            footPos = nan(1, nSamplesForce);
            if any(idx)
                footPos = squeeze(mean(MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin), 1, 'omitnan'))/1000;
            end
        end

        %% Store data

        Measurements.MotorData(fileNr).fileName = fname;
        Measurements.MotorData(fileNr).nSamplesForce = nSamplesForce;
        Measurements.MotorData(fileNr).sampleRateForce = sampleRateForce;
        Measurements.MotorData(fileNr).TimeForce = TimeForce;
        Measurements.MotorData(fileNr).Force = Force;
        Measurements.MotorData(fileNr).COP = COP;
        Measurements.MotorData(fileNr).idxContact = ~idxGaps;
        Measurements.MotorData(fileNr).startPos = startPos;
        Measurements.MotorData(fileNr).stopPos = stopPos;
        Measurements.MotorData(fileNr).nSamplesKin = nSamplesKin;
        Measurements.MotorData(fileNr).sampleRateKin = sampleRateKin;
        Measurements.MotorData(fileNr).TimeKin = TimeKin;
        Measurements.MotorData(fileNr).footPos = footPos;

        % report progress
        fprintf('\t-> %s (%d/%d = %.1f%% in %.3fs)\n', fname, iFile, nFiles, iFile/nFiles*100, toc);

        % increment number of processed files
        fileNr = fileNr+1;
    end

    % export Measurements structure to base workspace
    assignin('base', 'Measurements', Measurements);

    fprintf('Finished extracting data from %d files in %.3f s\n', fileNr - (length(loadedFiles) + 1), toc(ticAll));
    fprintf('If necessary, save current state using ''saveState''\n');

end

end