function [MeasurementInfo, MeasurementData] = prepareData(dataDir, paramDir, outDir, isLoadFiles)


if nargin < 4
    isLoadFiles = 1;
end

% load table containing subjects info
Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'));

% stages
stages = {'I', 'II', 'III'};

% list content of input folder into cell array
dirInfo = dir(dataDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxContact = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idxContact) = []; % remove folders and hidden files
nMeas = length(fpaths);

MeasurementInfo = table;
MeasurementData(nMeas, 1) = struct;
fprintf('Loading files...\n');
ticLoadingFiles = tic;
for iMeas = 1:nMeas
    measurement = table;
    fpath = fpaths{iMeas};
    [~, fileName, ~] = fileparts(fpath);
    if contains(fileName, 'ungueltig', 'IgnoreCase', true)
        continue
    end

    fprintf('\t-> %s\n', fileName);

    % split file name at underscores
    parts = strsplit(fileName, '_');

    % subject identity and code
    subjectName = parts{1};
    subjectCode = [parts{1}, '_', parts{2}];
    subjectCodes = [Subjects.Code_I, Subjects.Code_II, Subjects.Code_III];
    [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
    if isempty(subjectIdx)
        warning('Subject code %s not found -> skipping', subjectCode);
        continue
    end
    measurement.subjectName = string(subjectName); % store in measurement
    measurement.subjectCode = string(subjectCode);

    % stage
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    measurement.stage = stage; % store in measurement

    % subject properties at time of measurement
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(subjectIdx);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                measurement.(lower(subjectProps{iProp})) = propValue;
                break
            end
            myStage = myStage-1;
        end
    end

    % intervention
    [~, subjectStages] = find(contains(subjectCodes, subjectName));
    switch stageStr
        case 'I'
            intervention = 0;
        case 'II'
            if max(subjectStages) == 3 % subject has been tested at three stages
                intervention = 0;
            else
                intervention = 1;
            end
        case 'III'
            intervention = 1;
    end
    measurement.intervention = intervention;

    % task
    task = parts{3};
    measurement.task = string(task);

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        isMarker = 0;
        side = parts{5};
        trial = parts{6};
    else
        isMarker = 1;
        side = parts{4};
        trial = parts{5};
    end
    measurement.side = string(side);
    measurement.isMarker = isMarker;
    measurement.trial = trial;

    % store file name in measurement
    measurement.fileName = string(fileName);

    % set validity flag
    measurement.isValid = 1;

    if isLoadFiles
        tmp = load(fpath);
        fields = fieldnames(tmp);
        Data = tmp.(fields{1});
        lengthScale = 0.001; % mm -> m
        [Forces, Frequency, COPs, ~, Location, Analog] = kraftAusQualisys(Data, lengthScale);
        PowerLineHum = 50;
        nSamples = length(Forces);
        if nSamples < 2^nextpow2(PowerLineHum)/2
            fprintf('\t\tNot enough data (%d samples) -> skipping\n', nSamples);
            measurement.isValid = 0;
            continue
        end
        Forces = -Forces; % ground reaction forces are inverse of plate forces
        nPlates = size(Forces, 1);
        % nSamples = size(Forces, 3);
        dt = 1/Frequency; % time step size [s]
        LoadThresh = 20; % threshold below which forces are set to zero [N]
        InitForce = 200; % force threshold to cross to start trial [N]
        InitMarg = 0.5; % additional time after crossing InitForce to start trial [s]
        HumPeriod = Frequency/50;
        % CutoffFrequency = 20;
        for iPlate = 1:nPlates
            Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
        end
        [COPs, ~, Forces, ~] = getCOPfromAnalog(Analog, Forces, [], Location, Frequency, [], [], [], COPs);
        % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
        for iPlate = 1:nPlates
            idxContact = (abs(Forces(iPlate, 3, :)) < LoadThresh);
            Forces(iPlate, :, idxContact) = 0;
            COPs(iPlate, :, idxContact) = NaN;
        end
        Force = squeeze(sum(Forces, 1, 'omitnan'));
        absForces = abs(Forces(:, 3, :)); % z-component of force
        weightedCOPs = (absForces .* COPs) ./ sum(absForces, 1);
        COP = squeeze(sum(weightedCOPs, 1, 'omitnan'));

        iStart = find(abs(Force(3,:))>InitForce, 1, 'first');
        iStart = iStart + floor(InitMarg/dt);
        Force = Force(:, iStart:end);
        iStop = find(abs(Force(3,:))>InitForce, 1, 'last');
        iStop = iStop - floor(InitMarg/dt);
        Force = Force(:, 1:iStop);
        COP = COP(:, iStart:iStop+iStart-1);
        nSamples = size(Force, 2);
        Time = (0:nSamples-1)/Frequency;

        % remove non-contact phases
        idxContact = ~any(COP, 1);
        COP(:, idxContact) = NaN;
        idxContact = ~any(isnan(COP), 1);

        % distance to beam (for task "Balance")
        COPx = COP(1, idxContact);
        COPy = COP(2, idxContact);
        coefficients = polyfit(COPx, COPy, 1);
        yBeamFcn = @(x) polyval(coefficients, x);
        beamDistFcn = @(x, y) abs(y - yBeamFcn(x));
        beamDist = beamDistFcn(COP(1, :), COP(2, :));

        MeasurementData(iMeas).fileName = fileName;
        MeasurementData(iMeas).Force = Force;
        MeasurementData(iMeas).COP = COP;
        MeasurementData(iMeas).Time = Time;
        MeasurementData(iMeas).Frequency = Frequency;
        MeasurementData(iMeas).idxContact = idxContact;
        MeasurementData(iMeas).beamDist = beamDist;

        %% Plot data

        % setup figure
        nRows = 4;
        nCols = 1;
        iPlot = 0;
        fig = setupFigure(nCols*800, nRows*200, fileName);

        % plot COP path
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        hold on
        xFit = linspace(min(COPx), max(COPx), nSamples);
        yFit = yBeamFcn(xFit);
        scatter(xFit, yFit, 36, 'red', 'filled');
        scatter(COP(1, :), COP(2, :), 36, 'blue');
        title(sprintf('COP path of "%s"', fileName), 'Interpreter', 'none');
        xlabel('x [m]');
        ylabel('y [m]');        
        axis equal
        legend({'Beam', 'COP'});

        % plot COP components
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        plot(Time', COP');
        xlim([Time(1), Time(end)]);
        title(sprintf('COP components'), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Position [m]');
        legend({'x', 'y', 'z'});

        % plot GRF
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        plot(Time', Force');
        title(sprintf('Ground reaction force'), 'Interpreter', 'none');
        legend({'x', 'y', 'z'});
        xlabel('Time [s]');
        ylabel('Force [N]');
        xlim([Time(1), Time(end)]);

        % plot distance to beam
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        plot(Time', beamDistFcn(COP(1, :), COP(2, :))');
        xlim([Time(1), Time(end)]);
        title(sprintf('Distance to beam'), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Distance [m]');

        % save figure
        ftypes = {'png', 'pdf', 'fig'};
        for iType = 1:length(ftypes)
            ftype = ftypes{iType};
            myOutDir = fullfile(outDir, ftype);
            if ~isfolder(myOutDir)
                mkdir(myOutDir);
            end
            outpath = fullfile(myOutDir, fileName);
            saveFigure(fig, outpath, ftype);
        end
        close(fig);
    end

    %
    % switch measurement.task
    %     case 'Einbein'
    %     case 'Balance'
    %         % path length
    %         pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));
    %         precision = pathLength;
    %         COPx = COP(1, idxContact);
    %         COPy = COP(2, idxContact);
    %         coefficients = polyfit(COPx, COPy, 1);
    %         yBeamFcn = @(x) polyval(coefficients, x);
    %         xFit = linspace(min(COPx), max(COPx), nSamplesContact);
    %         yFit = yBeamFcn(xFit);
    %         beamDistFcn = @(x, y) abs(y - yBeamFcn(x));
    %         beamDist = beamDistFcn(COPx, COPy);
    %     case 'Sprung'
    %     otherwise
    %         continue
    % end
    %
    %

    MeasurementInfo = [MeasurementInfo; measurement]; %#ok<AGROW>
end
sprintf('Finished loading files in %f s.', toc(ticLoadingFiles));

if isLoadFiles
    outpath = fullfile(outDir, 'MeasurementData');
    save(outpath, 'MeasurementData');
    sprintf('Output saved to %s\n', outpath);
end

outpath = fullfile(outDir, 'MeasurementInfo.xlsx');
writetable(MeasurementInfo, outpath, 'WriteMode', 'replacefile');

end


