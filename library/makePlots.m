function makePlots

% load Measurements structure
load('Measurements.mat'); %#ok<LOAD>

outDir = Measurements.outDir;
nMeas = length(Measurements.Observations);

fprintf('Extract measurement data...\n');
ticAll = tic;
nProc = 0; % init number of processed files
for iMeas = 1:nMeas
    ticItem = tic;

    % get variables
    fileName = Measurements.Data(iMeas).fileName;
    Time = Measurements.Data(iMeas).Time;
    Force = Measurements.Data(iMeas).Force;
    COP = Measurements.Data(iMeas).COP;
    startPos = Measurements.Data(iMeas).startPos;
    stopPos = Measurements.Data(iMeas).stopPos;
    task = Measurements.Observations(iMeas).task;

    if ~(Measurements.Observations(iMeas).doneData && Measurements.Observations(iMeas).doneMetrics)
        continue
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);

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
            targetDist = Measurements.Data(iMeas).targetDist;
            Jerk = Measurements.Data(iMeas).Jerk;
            targetError = Measurements.Observations(iMeas).targetError;
            nSamples = length(targetDist);
            idxContact = Measurements.Data(iMeas).idxContact;

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
            plot(Time', targetDist');
            hold on
            yline(targetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
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
            targetDist = Measurements.Data(iMeas).targetDist;
            Jerk = Measurements.Data(iMeas).Jerk;
            targetError = Measurements.Observations(iMeas).targetError;
            meanCOP = mean(COP, 2);

            % plot COP path with center
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            hold on
            scatter(COP(1, :), COP(2, :), 2, 'blue');
            scatter(meanCOP(1), meanCOP(2), 24, 'red', 'x')
            title(sprintf('COP path'), 'Interpreter', 'none');
            xlabel('x [m]');
            ylabel('y [m]');
            axis equal
            legend({'COP', 'center'});

            % plot distance to center
            iPlot = iPlot+1;
            subplot(nRows, nCols, plotIdx(iPlot));
            plot(Time', targetDist');
            hold on
            yline(targetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
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
            targetDist = Measurements.Data(iMeas).targetDist;
            Jerk = Measurements.Data(iMeas).Jerk;
            targetError = Measurements.Observations(iMeas).targetError;
            jumpStopPos = Measurements.Data(iMeas).jumpStopPos;
            [~, targetIdx] = min(targetDist);

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
            plot(Time', targetDist');
            hold on
            xline(Time(targetIdx), 'red');
            xlim([Time(1), Time(end)]);
            yl = ylim;
            text(Time(end), max(yl), sprintf('targetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
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

    % save figure
    ftypes = {'pdf', 'png'};
    for iType = 1:length(ftypes)
        ftype = ftypes{iType};
        myOutDir = fullfile(outDir, ['figures_', ftype]);
        if ~isfolder(myOutDir)
            mkdir(myOutDir);
        end
        outpath = fullfile(myOutDir, fileName);
        saveFigure(fig, outpath, ftype);
    end
    close(fig);

    fprintf('\t\tFinished in %.3f s\n', toc(ticItem));
end

fprintf('Finished plotting from %d datasets in %.3f s\n\n', nProc, toc(ticAll));

end