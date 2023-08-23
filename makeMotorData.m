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
        stage = parts{2};
        subjectCode = [subject, '_', stage];
        subjectCodes = [SubjectsTable.Code_I, SubjectsTable.Code_II, SubjectsTable.Code_III];
        [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
        if isempty(subjectIdx)
            % fprintf('\t%s:\tNo subject with code ''%s'' -> skipping\n', fname, subjectCode);
            continue
        end

        % jump position markers
        Beidbein_start = string(SubjectsTable.Beidbein_start(subjectIdx));
        Beidbein_stop = string(SubjectsTable.Beidbein_stop(subjectIdx));
        switch stage
            case 'I'
                Einbein_start = string(SubjectsTable.Einbein_start_I(subjectIdx));
                Einbein_stop = string(SubjectsTable.Einbein_stop_I(subjectIdx));
            case {'II', 'III'}
                Einbein_start = string(SubjectsTable.Einbein_start_II_III(subjectIdx));
                Einbein_stop = string(SubjectsTable.Einbein_stop_II_III(subjectIdx));
        end

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
            [Forces, sampleRate, COPs, ~, location, Analog] = kraftAusQualisys(MotorData, lengthScale);
        catch ME
            fprintf('\t%s:\tError when extracting forces -> skipping\n', fname);
            fprintf('\t\t%s\n', ME.message);
            continue
        end
        PowerLineHum = 50;
        nSamples = length(Forces);
        if nSamples < 2^nextpow2(PowerLineHum)/2
            fprintf('\t%s:\tNot enough data (%d samples) -> skipping\n', fname, nSamples);
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
        iStart = iStart + floor(InitMarg/dt);
        Force = Force(:, iStart:end);
        iStop = find(abs(Force(3,:))>InitForce, 1, 'last');
        iStop = iStop - floor(InitMarg/dt);
        Force = Force(:, 1:iStop);
        idxForce =  iStart:iStop+iStart-1;
        COP = COP(:, idxForce);
        nSamples = size(Force, 2);
        if nSamples == 0
            fprintf('\t%s:\tToo few samples with nonzero force -> skipping\n', fname);
            continue
        end
        Time = (0:nSamples-1)/sampleRate;
        sampleRateKin = MotorData.FrameRate;        
        iStartKin = ceil(iStart/sampleRate*sampleRateKin);
        iStopKin = ceil(iStop/sampleRate*sampleRateKin);
        idxKin =  iStartKin:iStopKin+iStartKin-1;   
        nSamplesKin = length(idxKin);
        TimeKin = (0:nSamplesKin-1)/sampleRateKin;

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
            % get start and stop marker positions
            idx = strcmp(labels, startLabel);
            if any(idx)
                pos_raw = MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin)/1000;
                pos_raw = squeeze(mean(pos_raw, 1, 'omitnan')); % avg over markers
                startPos = mean(pos_raw, 2, 'omitnan');
            end
            idx = strcmp(labels, stopLabel);
            if any(idx)
                pos_raw = MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin)/1000;
                pos_raw = squeeze(mean(pos_raw, 1, 'omitnan')); % avg over markers
                stopPos = mean(pos_raw, 2, 'omitnan');
            end
            % get foot marker data
            idx = ismember(labels, footLabels);
            footPos = nan(3, nSamples);
            if any(idx)
                pos_raw = MotorData.Trajectories.Labeled.Data(idx, 1:3, idxKin)/1000;
                pos_raw = squeeze(mean(pos_raw, 1, 'omitnan')); % avg over markers
                idxNan = any(isnan(pos_raw),1);
                if ~all(idxNan)
                    if any(idxNan)
                        pos_raw(:, idxNan) = interp1(TimeKin(~idxNan)', pos_raw(:,~idxNan)', TimeKin(idxNan)', 'pchip')';
                    end
                    footPos = interp1(TimeKin', pos_raw', Time', 'pchip')';
                end
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