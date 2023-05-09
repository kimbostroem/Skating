% set input and output folder
inDir = '/Users/kbostroem/sciebo/Skating/Skating_In';
outDir = '/Users/kbostroem/sciebo/Skating/Skating_Out';

% load library
addpath(genpath('library'));

% load table containing subjects info
Subjects = readtable('Subjects.xlsx');

% list content of input folder into cell array
dirInfo = dir(inDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxContact = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idxContact) = []; % remove folders and hidden files
nObs = length(fpaths);
% nObs = 2;
observations = table;
for iObs = 1:nObs
    observation = table;
    fpath = fpaths{iObs};
    [~, fname, ~] = fileparts(fpath);
    if contains(fname, 'ungueltig', 'IgnoreCase', true)
        continue
    end

    parts = strsplit(fname, '_');

    % get subject name
    switch parts{2}
        case 'II'
            subjectCode = [parts{1}, '_', parts{2}];
            parts = parts(3:end);
        case 'III'
            subjectCode = [parts{1}, '_', parts{2}];
            parts = parts(3:end);
        otherwise
            subjectCode = parts{1};
            parts = parts(2:end);
    end
    subjectCodes = [Subjects.Code, Subjects.Code_post, Subjects.Code_final];
    [row, col] = find(strcmp(subjectCodes, subjectCode));
    if ~isempty(row)
        switch col
            case 1
                subject = subjectCodes{row, 1};
            case 2
                subject = subjectCodes{row, 1};
            case 3
                subject = subjectCodes{row, 1};
        end
    else
        warning('Subject %s not found -> skipping', subjectCode);
    end

    % observation.subject = string(parts{1});
    % observation.task = string(parts{2});
    % observation.side = 'B';
    % observation.time = 1;
    % observation.trial = str2double(parts{end});
    % if length(parts) > 3
    %     if strcmp(parts{3}, 'Post')
    %         observation.time = 2;
    %         observation.side = parts{4};
    %     else
    %         observation.time = 1;
    %         observation.side = string(parts{3});
    %     end
    % end
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

    observations = [observations; observation]; %#ok<AGROW>
end

% subjects = unique(observations.subject);



