function do4_makeMetrics(Measurements)

if nargin < 1
    load('Measurements'); %#ok<LOAD>
end

nMeas = length(Measurements.Observations);
outDir = Measurements.outDir;

fprintf('Create metrics...\n');
ticAll = tic;
for iMeas = 1:nMeas
    fileName = Measurements.Observations(iMeas).fileName;
    
    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);

    subjectWeight = Measurements.Observations(iMeas).weight;
    Time = Measurements.Data(iMeas).Time;
    Force = Measurements.Data(iMeas).Force;
    COP = Measurements.Data(iMeas).COP;
    idxContact = Measurements.Data(iMeas).idxContact;
    sampleRate = Measurements.Data(iMeas).sampleRate;
    dt = 1/sampleRate; % time step size [s]
    nSamples = Measurements.Data(iMeas).nSamples;
    task = Measurements.Observations(iMeas).task;
    startPos = Measurements.Data(iMeas).startPos;
    stopPos = Measurements.Data(iMeas).stopPos;

    switch task

        case {'Balance', 'Einbein'}

            if strcmp(task, 'Balance')
                % distance to beam (for task "Balance")
                coefficients = polyfit(COP(1, idxContact), COP(2, idxContact), 1);
                targetFcn = @(x, y) polyval(coefficients, x);
                targetDistFcn = @(x, y) abs(y - targetFcn(x, y));
            elseif strcmp(task, 'Einbein')
                meanCOP = mean(COP, 2);
                targetDistFcn = @(x, y) vecnorm([x; y] - meanCOP);
            end
            targetDist = targetDistFcn(COP(1, :), COP(2, :));

            % jerk
            dJerk = diff(Force, 1, 2);
            Jerk = [dJerk, dJerk(:, end)] / (dt * subjectWeight);
            Jerk(:, ~idxContact) = 0; % remove jerk around gaps

            % path length
            pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));

            % mean beam distance
            targetError = mean(targetDist, 'omitnan');

            % mean jerk
            meanJerk = mean(vecnorm(Jerk, 2, 1), 'omitnan');
            meanJerkXY = mean(vecnorm(Jerk(1:2, :), 2, 1), 'omitnan');
            Measurements.Data(iMeas).targetDist = targetDist;
            Measurements.Data(iMeas).Jerk = Jerk;
            Measurements.Observations(iMeas).pathLength = pathLength;
            Measurements.Observations(iMeas).targetError = targetError;
            Measurements.Observations(iMeas).meanJerk = meanJerk;
            Measurements.Observations(iMeas).meanJerkXY = meanJerkXY;

        case 'Sprung'
            % distance to jump stop position
            targetDist = abs(Measurements.Data(iMeas).COP(1, :) - stopPos(1));
            [targetError, targetIdx] = min(targetDist);
            jumpStopPos = Measurements.Data(iMeas).COP(:, targetIdx);
            Measurements.Data(iMeas).targetDist = targetDist;
            Measurements.Observations(iMeas).targetError = targetError;
            Measurements.Data(iMeas).jumpStopPos = jumpStopPos;

            % jerk
            dJerk = diff(Force, 1, 2);
            Jerk = [dJerk, dJerk(:, end)] / (dt * subjectWeight);
            Jerk(:, ~idxContact) = 0; % remove jerk around gaps

            % path length
            pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));

            % mean beam distance
            targetError = mean(targetDist, 'omitnan');

            % mean jerk
            meanJerk = mean(vecnorm(Jerk, 2, 1), 'omitnan');
            meanJerkXY = mean(vecnorm(Jerk(1:2, :), 2, 1), 'omitnan');
            Measurements.Data(iMeas).Jerk = Jerk;
            Measurements.Observations(iMeas).pathLength = pathLength;
            Measurements.Observations(iMeas).meanJerk = meanJerk;
            Measurements.Observations(iMeas).meanJerkXY = meanJerkXY;
    end


    %% Store data

    

    %% Plot data

    % setup figure
    nRows = 3;
    nCols = 1;
    iPlot = 0;
    fig = setupFigure(nCols*400, nRows*200, fileName);

    switch task

        case 'Balance'
            % plot COP path with beam
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
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
            subplot(nRows, nCols, iPlot);
            plot(Time', targetDistFcn(COP(1, :), COP(2, :))');
            hold on
            yline(targetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('\ttargetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            xlim([Time(1), Time(end)]);
            title(sprintf('Distance to beam'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
            plot(Time', Jerk');
            xlim([Time(1), Time(end)]);
            title(sprintf('Jerk'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Jerk [m/s^3]');
            legend({'x', 'y', 'z'});

        case 'Einbein'
            % plot COP path with center
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
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
            subplot(nRows, nCols, iPlot);
            plot(Time', targetDistFcn(COP(1, :), COP(2, :))');
            hold on
            yline(targetError, 'r');
            yl = ylim;
            text(Time(end), max(yl), sprintf('\ttargetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            xlim([Time(1), Time(end)]);
            title(sprintf('Distance to center'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');            

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
            plot(Time', Jerk');
            xlim([Time(1), Time(end)]);
            title(sprintf('Jerk'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Jerk [m/s^3]');
            legend({'x', 'y', 'z'});

        case 'Sprung'
            % plot COP path with jump stop position
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
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
            subplot(nRows, nCols, iPlot);
            plot(Time', targetDist');
            hold on
            xline(Time(targetIdx), 'red');
            xlim([Time(1), Time(end)]);
            yl = ylim;
            text(Time(end), max(yl), sprintf('\ttargetError = %.0f mm', targetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            title(sprintf('Distance from landing to jump stop pos'), 'Interpreter', 'none');
            xlabel('Time [s]');
            ylabel('Distance [m]');

            % plot jerk
            iPlot = iPlot+1;
            subplot(nRows, nCols, iPlot);
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
    ftypes = {'png', 'pdf', 'fig'};
    for iType = 1:length(ftypes)
        ftype = ftypes{iType};
        myOutDir = fullfile(outDir, 'metrics', ftype);
        if ~isfolder(myOutDir)
            mkdir(myOutDir);
        end
        outpath = fullfile(myOutDir, fileName);
        saveFigure(fig, outpath, ftype);
    end
    close(fig);
end
sprintf('Finished in %f s\n\n', toc(ticAll));

%% Saving table

fprintf('Saving Observations to table...\n');
MeasurementTable = struct2table(Measurements.Observations);
outpath = fullfile(outDir, 'Observations.xlsx');
writetable(MeasurementTable, outpath, 'WriteMode', 'replacefile');

%% Export Measurements structure to base workspace

assignin('base', 'Measurements', Measurements);

%% Save Measurements structure to MAT file

fprintf('Saving Measurements to MAT file...\n');
save('Measurements', 'Measurements');
fprintf('Moving Measurements.mat to output folder...\n');
movefile('Measurements.mat', outDir);
fprintf('DONE\n\n');

end