function makeData

% load current state
loadState;

gapMargin = 0.1; % margin around gaps in seconds, because there the COP is distorted
LoadThresh = 20; % threshold below which forces are set to zero [N]
InitForce = 200; % force threshold to cross to start trial [N]
InitMarg = 0.5; % additional time after crossing InitForce to start trial [s]
HumFreq = 50; % humming frequency in Hz

inDir = Measurements.inDir; %#ok<NODEF>
outDir = Measurements.outDir;
nMeas = length(Measurements.Observations);

fprintf('Extract measurement data...\n');
ticAll = tic;
nProc = 0; % init number of processed files
for iMeas = 1:nMeas
    ticItem = tic;
    fileName = Measurements.Observations(iMeas).fileName;
    fpath = fullfile(inDir, fileName);
    task = Measurements.Observations(iMeas).task;
    side = Measurements.Observations(iMeas).side;

    if Measurements.Observations(iMeas).doneData
        continue
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);

    %% Import Force data

    Measurements.Data(iMeas).isForce = 1;
    tmp = load(fpath);
    fields = fieldnames(tmp);
    Data = tmp.(fields{1});
    lengthScale = 0.001; % mm -> m
    [Forces, sampleRate, COPs, ~, location, Analog] = kraftAusQualisys(Data, lengthScale);
    PowerLineHum = 50;
    nSamples = length(Forces);
    if nSamples < 2^nextpow2(PowerLineHum)/2
        fprintf('\t\tNot enough data (%d samples) -> skipping\n', nSamples);
        Measurements.Data(iMeas).isForce = 0;
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
    if isfield(Data, 'Trajectories') && isfield(Data.Trajectories, 'Labeled')
        labels = Data.Trajectories.Labeled.Labels;
        % get start and stop marker labels
        switch side
            case 'B'
                startLabel = Measurements.Observations(iMeas).Beidbein_start;
                stopLabel = Measurements.Observations(iMeas).Beidbein_stop;
            case {'L', 'R'}
                startLabel = Measurements.Observations(iMeas).Einbein_start;
                stopLabel = Measurements.Observations(iMeas).Einbein_stop;
            otherwise
        end
        % get start and stop marker data
        idx = strcmp(labels, startLabel);
        if any(idx)
            pos_raw = squeeze(Data.Trajectories.Labeled.Data(idx, 1:3, :))/1000;
            startPos = mean(pos_raw, 2);
        end
        idx = strcmp(labels, stopLabel);
        if any(idx)
            pos_raw = squeeze(Data.Trajectories.Labeled.Data(idx, 1:3, :))/1000;
            stopPos = mean(pos_raw, 2);
        end
    end

    %% Store data

    Measurements.Data(iMeas).fileName = fileName;
    Measurements.Data(iMeas).nSamples = nSamples;
    Measurements.Data(iMeas).sampleRate = sampleRate;
    Measurements.Data(iMeas).Time = Time;
    Measurements.Data(iMeas).Force = Force;
    Measurements.Data(iMeas).COP = COP;
    Measurements.Data(iMeas).idxContact = ~idxGaps;
    Measurements.Data(iMeas).startPos = startPos;
    Measurements.Data(iMeas).stopPos = stopPos;    

    % increment number of processed files
    nProc = nProc+1;

    % set flag
    Measurements.Observations(iMeas).doneData = 1;

    % export Measurements structure to base workspace
    fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
    assignin('base', 'Measurements', Measurements);

    fprintf('\t\tFinished in %.3f s\n', toc(ticItem));
end

% save current state
saveState;

fprintf('Finished extracting data from %d files in %.3f s\n\n', nProc, toc(ticAll));

end