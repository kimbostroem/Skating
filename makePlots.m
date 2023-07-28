function makePlots

% load current state
Measurements = loadState();

outDir = evalin('base', 'outDir');

fprintf('Make plots...\n');
ticAll = tic;

% create figure output folders
figTypes = {'pdf', 'png'};
nFigTypes = length(figTypes);
figDirs = cell(nFigTypes, 1);
for iFigType = 1:nFigTypes
    figType = figTypes{iFigType};
    figDir = fullfile(outDir, ['figures_', figType]);
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    figDirs{iFigType} = figDir;
end

nProc = 0; % init number of processed files
nMeas = length(Measurements.MotorData);
for iMeas = 1:nMeas
    ticItem = tic;

    % get variables
    fileName = Measurements.MotorData(iMeas).fileName;

    % if figure already exists -> skip
    outpaths = fullfile(figDirs, strcat([fileName, '.'], figTypes)');
    if all(isfile(outpaths))
        continue
    end

     % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);
    
    % get data
    Time = Measurements.MotorData(iMeas).Time;
    Force = Measurements.MotorData(iMeas).Force;
    COP = Measurements.MotorData(iMeas).COP;
    startPos = Measurements.MotorData(iMeas).startPos;
    stopPos = Measurements.MotorData(iMeas).stopPos;
    task = Measurements.MotorMetrics.Task(iMeas);   

    % setup figure
    nRows = 3;
    nCols = 2;
    iPlot = 0;
    fig = setupFigure(nCols*400, nRows*200, fileName);
    % make subplot go down columns instead of rows
    plotIdx = reshape(1:nRows*nCols, nCols, nRows).';

    %% Plot imported data

    % plot COP path
    iPlot = iPlot+1;
    subplot(nRows, nCols, plotIdx(iPlot));
    scatter(COP(1, :), COP(2, :), 2, 'blue');    
    title(sprintf('COP path'), 'Interpreter', 'none');
    xlabel('x [m]');
    ylabel('y [m]');
    axis equal

    % plot COP components
    iPlot = iPlot+1;
    subplot(nRows, nCols, plotIdx(iPlot));
    plot(Time', COP');
    xlim([Time(1), Time(end)]);
    title(sprintf('COP components'), 'Interpreter', 'none');
    xlabel('Time [s]');
    ylabel('Position [m]');
    legend({'x', 'y'});

    % plot GRF
    iPlot = iPlot+1;
    subplot(nRows, nCols, plotIdx(iPlot));
    hold on
    plot(Time', Force');
    title(sprintf('Ground reaction force'), 'Interpreter', 'none');
    legend({'x', 'y', 'z'});
    xlabel('Time [s]');
    ylabel('Force [N]');
    xlim([Time(1), Time(end)]);

    %% Plot metrics

    switch task

        case 'Balance'

            % get variables
            Deviation = Measurements.MotorData(iMeas).Deviation;
            Jerk = Measurements.MotorData(iMeas).Jerk;
            TargetError = Measurements.MotorMetrics.TargetError(iMeas);
            nSamples = length(Deviation);
            idxContact = Measurements.MotorData(iMeas).idxContact;

            % get targetFcn and targetDistFcn
            coefficients = polyfit(COP(1, idxContact), COP(2, idxContact), 1);
            targetFcn = @(x, y) polyval(coefficients, x);

            % plot COP path with beam
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            hold on
            xFit = linspace(min(COP(1, :)), max(COP(1, :)), nSamples);
            yFit = targetFcn(xFit);
            scatter(xFit, yFit, 2, 'red', '.');
            scatter(COP(1, :), COP(2, :), 2, 'blue');
            title(sprintf('COP path'), 'Interpreter', 'none');
            xlabel('x [m]');
            ylabel('y [m]');
            axis equal
            legend({'Beam', 'COP'});

            % plot distance to beam
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Deviation');
            hold on
            yline(TargetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', TargetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            xlim([Time(1), Time(end)]);
            title(sprintf('Distance to beam'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Jerk');
            xlim([Time(1), Time(end)]);
            title(sprintf('Jerk'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Jerk [m/s^3]');
            legend({'x', 'y', 'z'});

        case 'Einbein'
            % get variables
            Deviation = Measurements.MotorData(iMeas).Deviation;
            Jerk = Measurements.MotorData(iMeas).Jerk;
            TargetError = Measurements.MotorMetrics.TargetError(iMeas);
            MeanCOP = mean(COP, 2);

            % plot COP path with center
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            hold on
            scatter(COP(1, :), COP(2, :), 2, 'blue');
            scatter(MeanCOP(1), MeanCOP(2), 24, 'red', 'x')
            title(sprintf('COP path'), 'Interpreter', 'none');
            xlabel('x [m]');
            ylabel('y [m]');
            axis equal
            legend({'COP', 'center'});

            % plot distance to center
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Deviation');
            hold on
            yline(TargetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', TargetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            xlim([Time(1), Time(end)]);
            title(sprintf('Distance to center'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Jerk');
            xlim([Time(1), Time(end)]);
            title(sprintf('Jerk'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Jerk [m/s^3]');
            legend({'x', 'y', 'z'});

        case 'Sprung'
            % get variables
            Deviation = Measurements.MotorData(iMeas).Deviation;
            Jerk = Measurements.MotorData(iMeas).Jerk;
            TargetError = Measurements.MotorMetrics.TargetError(iMeas);
            jumpStopPos = Measurements.MotorData(iMeas).jumpStopPos;
            [~, targetIdx] = min(Deviation);

            % plot COP path with jump stop position
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            hold on
            scatter(COP(1, :), COP(2, :), 2, 'blue');
            scatter(startPos(1,:)', startPos(2,:)', 5, 'green', 'filled');
            scatter(stopPos(1,:)', stopPos(2,:)', 5, 'red', 'filled');
            scatter(jumpStopPos(1), jumpStopPos(2), 10, 'red', 'x')
            xline(startPos(1), 'g');
            xline(stopPos(1), 'r');
            title(sprintf('COP path'), 'Interpreter', 'none');
            xlabel('x [m]');
            ylabel('y [m]');
            axis equal
            legend({'COP', 'jump start pos', 'jump stop pos', 'landing'});

            % plot distance to jump stop position
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Deviation');
            hold on
            xline(Time(targetIdx), 'red');
            xlim([Time(1), Time(end)]);
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', TargetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            title(sprintf('Distance from landing to jump stop pos'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', Jerk');
            xlim([Time(1), Time(end)]);
            title(sprintf('Jerk'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Jerk [m/s^3]');
            legend({'x', 'y', 'z'});

    end

    % global figure title
    sgtitle(sprintf('%s', fileName), 'Interpreter', 'none');

    % increment number of processed files
    nProc = nProc+1;
    
    % save figure
    for iFigType = 1:nFigTypes
        figType = figTypes{iFigType};
        figDir = figDirs{iFigType};
        outpath = fullfile(figDir, fileName);
        saveFigure(fig, outpath, figType);
    end
    close(fig);

    % report item finish
    fprintf('\t\tFinished in %.3f s\n', toc(ticItem));
end

fprintf('Finished plotting from %d datasets in %.3f s\n\n', nProc, toc(ticAll));

end