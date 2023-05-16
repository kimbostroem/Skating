% set input and output folder
inDir = '/Users/kbostroem/sciebo/Skating/Skating_In';
outDir = '/Users/kbostroem/sciebo/Skating/Skating_Out';
paramDir = '/Users/kbostroem/sciebo/Skating/Auswertung';

% load library
addpath(genpath('library'));

% load table containing subjects info
Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'));

% stages
stages = {'I', 'II', 'III'};

% list content of input folder into cell array
dirInfo = dir(inDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxContact = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idxContact) = []; % remove folders and hidden files
nObs = length(fpaths);
% nObs = 2;
Observations = table;
for iObs = 1:nObs
    observation = table;
    fpath = fpaths{iObs};
    [~, fileName, ~] = fileparts(fpath);
    if contains(fileName, 'ungueltig', 'IgnoreCase', true)
        continue
    end   

    % split file name at underscores
    parts = strsplit(fileName, '_');

    % subject identity and code
    subjectName = parts{1};
    subjectCode = [parts{1}, '_', parts{2}];
    subjectCodes = [Subjects.Code_I, Subjects.Code_II, Subjects.Code_III];
    [row, col] = find(strcmp(subjectCodes, subjectCode));
    if isempty(row)
        warning('Subject code %s not found -> skipping', subjectCode);
        continue
    end    
    observation.subjectName = string(subjectName); % store in observation
    observation.subjectCode = string(subjectCode);

    % stage
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    observation.stage = stage; % store in observation    

    % subject properties at time of measurement
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(row);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                observation.(lower(subjectProps{iProp})) = propValue;
                break
            end
            myStage = myStage-1;
        end
    end

    % intervention    
    [row, col] = find(contains(subjectCodes, subjectName));   
    switch stageStr
        case 'I'
            intervention = 0;
        case 'II'
            if max(col) == 3 % subject has been tested at three stages
                intervention = 0;
            else
                intervention = 1;
            end
        case 'III'
            intervention = 1;
    end
    observation.intervention = intervention;

    % task
    task = parts{3};
    observation.task = string(task);

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
    observation.side = string(side);
    observation.isMarker = isMarker;
    observation.trial = trial;

    % store file name in observation
    observation.fileName = string(fileName);
    
    %
    % tmp = load(fpath);
    % fields = fieldnames(tmp);
    % Data = tmp.(fields{1});
    % lengthScale = 0.001; % mm -> m
    % [Forces, Frequency, COPs, ~, Location, Analog] = kraftAusQualisys(Data, lengthScale);
    % Forces = -Forces; % ground reaction forces are inverse of plate forces
    % nPlates = size(Forces, 1);
    % nSamples = size(Forces, 3);
    % dt = 1/Frequency; % time step size [s]
    % LoadThresh = 20; % threshold below which forces are set to zero [N]
    % InitForce = 200; % force threshold to cross to start trial [N]
    % InitMarg = 0.5; % additional time after crossing InitForce to start trial [s]
    % HumPeriod = Frequency/50;
    % CutoffFrequency = 20;
    % Time = (0:nSamples-1)/Frequency;
    % for iPlate = 1:nPlates
    %     Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
    % end
    % [COPs, ~, Forces, Loaded] = getCOPfromAnalog(Analog, Forces, [], Location, Frequency, [], [], [], COPs);
    % % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
    % for iPlate = 1:nPlates
    %     idxContact = (abs(Forces(iPlate, 3, :)) < LoadThresh);
    %     Forces(iPlate, :, idxContact) = 0;
    %     COPs(iPlate, :, idxContact) = NaN;
    % end
    % Force = squeeze(sum(Forces, 1, 'omitnan'));
    % absForces = abs(Forces(:, 3, :)); % z-component of force
    % weightedCOPs = (absForces .* COPs) ./ sum(absForces, 1);
    % COP = squeeze(sum(weightedCOPs, 1, 'omitnan'));
    %
    % iStart = find(abs(Force(3,:))>InitForce, 1, 'first');
    % iStart = iStart + floor(InitMarg/dt);
    % Force = Force(:, iStart:end);
    % iStop = find(abs(Force(3,:))>InitForce, 1, 'last');
    % iStop = iStop - floor(InitMarg/dt);
    % Force = Force(:, 1:iStop);
    % COP = COP(:, iStart:iStop+iStart-1);
    % nSamples = size(Force, 2);
    % Time = (0:nSamples-1)/Frequency;
    %
    %
    % % remove non-contact phases
    % idxContact = ~any(COP, 1);
    % COP(:, idxContact) = NaN;
    % idxContact = ~any(isnan(COP), 1);
    % nSamplesContact = sum(idxContact);
    %
    % switch observation.task
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
    % %% Plot data
    %
    % nRows = 4;
    % nCols = 1;
    % iPlot = 0;
    % fig = setupFigure(nCols*800, nRows*200, fname);
    %
    % % plot COP path
    % iPlot = iPlot+1;
    % subplot(nRows, nCols, iPlot);
    % scatter(COP(1, :), COP(2, :), 'bo');
    % title(sprintf('COP path of "%s"', fname), 'Interpreter', 'none');
    % xlabel('x [m]');
    % ylabel('y [m]');
    % hold on
    % scatter(COP(1, :), yBeamFcn(COP(1, :)), '.r');
    % axis equal
    %
    % % plot COP components
    % iPlot = iPlot+1;
    % subplot(nRows, nCols, iPlot);
    % plot(Time', COP');
    % xlim([Time(1), Time(end)]);
    % title(sprintf('COP components'), 'Interpreter', 'none');
    % xlabel('Time [s]');
    % ylabel('Position [m]');
    % legend({'x', 'y', 'z'});
    %
    % % plot distance to beam
    % iPlot = iPlot+1;
    % subplot(nRows, nCols, iPlot);
    % plot(Time', beamDistFcn(COP(1, :), COP(2, :))');
    % xlim([Time(1), Time(end)]);
    % title(sprintf('Distance to beam'), 'Interpreter', 'none');
    % xlabel('Time [s]');
    % ylabel('Distance [m]');
    %
    % % plot GRF
    % iPlot = iPlot+1;
    % subplot(nRows, nCols, iPlot);
    % plot(Time', Force');
    % title(sprintf('Ground reaction force'), 'Interpreter', 'none');
    % legend({'x', 'y', 'z'});
    % xlabel('Time [s]');
    % ylabel('Force [N]');
    % xlim([Time(1), Time(end)]);
    % filePath = fullfile(outDir, fname);
    % saveFigure(fig, filePath, 'png');
    % close(fig);

    Observations = [Observations; observation]; %#ok<AGROW>
end

outpath = fullfile(outDir, 'Observations.xlsx');
writetable(Observations, outpath, 'WriteMode', 'replacefile');


