doPlot = 1;
inDir = fullfile('..', 'Skating_In');
outDir = fullfile('..', 'Skating_Out');

% list content of input folder into cell array
dirInfo = dir(inDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxContact = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idxContact) = []; % remove folders and hidden files
nObs = length(fpaths);
nObs = 1;
observations = table;
for iObs = 1:nObs
    observation = table;
    fpath = fpaths{iObs};
    [~, fname, ~] = fileparts(fpath);
    if contains(fname, 'ungueltig', 'IgnoreCase', true)
        continue
    end
    parts = strsplit(fname, '_');
    observation.subject = string(parts{1});
    observation.task = string(parts{2});
    observation.side = 'B';
    observation.time = 1;
    observation.trial = str2double(parts{end});
    if length(parts) > 3
        if strcmp(parts{3}, 'Post')
            observation.time = 2;
            observation.side = parts{4};
        else
            observation.time = 1;
            observation.side = string(parts{3});
        end
    end

    tmp = load(fpath);
    fields = fieldnames(tmp);
    Data = tmp.(fields{1});
    lengthScale = 0.001; % mm -> m
    [Forces, Frequency, COPs, ~, Location, Analog] = kraftAusQualisys(Data, lengthScale);
    Forces = -Forces; % ground reaction forces are inverse of plate forces
    nPlates = size(Forces, 1);
    nSamples = size(Forces, 3);
    LoadThresh = 20;
    InitForce = 200;
    HumPeriod = Frequency/50;
    CutoffFrequency = 20;
    Time = (0:nSamples-1)/Frequency;
    for iPlate = 1:nPlates
        Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
    end
    [COPs, ~, Forces, Loaded] = getCOPfromAnalog(Analog, Forces, [], Location, Frequency, [], [], [], COPs);
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

    start = find(abs(Force(3,:))>InitForce, 1, 'first');
    Force = Force(:, start:end);
    stop = find(abs(Force(3,:))<InitForce, 1, 'first');
    Force = Force(:, 1:stop);
    COP = COP(:, start:stop+start-1);
    Time = Time(start:stop+start-1);


    % remove non-contact phases
    idxContact = ~any(COP, 1);
    COP(:, idxContact) = NaN; 
    idxContact = ~any(isnan(COP), 1);
    nSamplesContact = sum(idxContact);

    switch observation.task
        case 'Einbein'
        case 'Balance'
            % path length            
            pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));
            precision = pathLength;
            COPx = COP(1, idxContact);
            COPy = COP(2, idxContact);
            coefficients = polyfit(COPx, COPy, 1);            
            yBeamFcn = @(x) polyval(coefficients, x);
            xFit = linspace(min(COPx), max(COPx), nSamplesContact);
            yFit = yBeamFcn(xFit);
            beamDistFcn = @(x, y) abs(y - yBeamFcn(x));
            beamDist = beamDistFcn(COPx, COPy);
        case 'Sprung'
        otherwise
            continue
    end


    if doPlot
        nRows = 4;
        nCols = 1;
        iPlot = 0;
        fig = setupFigure(nCols*800, nRows*200, fname);

        % plot COP path
        iPlot = iPlot+1;        
        subplot(nRows, nCols, iPlot);
        scatter(COP(1, :), COP(2, :), 'bo');
        title(sprintf('COP path of "%s"', fname), 'Interpreter', 'none');
        xlabel('x [m]');
        ylabel('y [m]');
        hold on
        scatter(COP(1, :), yBeamFcn(COP(1, :)), '.r');
        axis equal
        
        % plot COP components
        iPlot = iPlot+1;        
        subplot(nRows, nCols, iPlot);
        plot(Time', COP');
        xlim([Time(1), Time(end)]);
        title(sprintf('COP components'), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Position [m]');
        legend({'x', 'y', 'z'});
        
        % plot distance to beam
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        plot(Time', beamDistFcn(COP(1, :), COP(2, :))');
        xlim([Time(1), Time(end)]);
        title(sprintf('Distance to beam'), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Distance [m]');

        % plot GRF
        iPlot = iPlot+1;
        subplot(nRows, nCols, iPlot);
        plot(Time', Force');
        title(sprintf('Ground reaction force'), 'Interpreter', 'none');
        legend({'x', 'y', 'z'});
        xlabel('Time [s]');
        ylabel('Force [N]');
        xlim([Time(1), Time(end)]);
        filePath = fullfile(outDir, fname);
        saveFigure(fig, filePath, 'png');
        close(fig);
    end

    observations = [observations; observation]; %#ok<AGROW>
end

subjects = unique(observations.subject);



