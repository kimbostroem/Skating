function fig = makePlot(iMeas, Measurements)

if ~isfield(Measurements, 'MotorMetrics')
    fprintf('\tNo field ''MotorMetrics'' to plot -> abort\n');
    fig = [];
    return
end

ticItem = tic;

% get variables
fileName = Measurements.MotorData(iMeas).fileName;

% report progress
fprintf('\t\t Plot %s (iMeas = %d)\n', fileName, iMeas);

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
        Jerk = Measurements.MotorData(iMeas).Jerk;
        TargetError = Measurements.MotorMetrics.TargetError(iMeas);
        landingPos = Measurements.MotorData(iMeas).landingPos;
        footPos = Measurements.MotorData(iMeas).footPos;

        % plot COP path with jump stop position
        iPlot = iPlot+1;
        subplot(nRows, nCols, plotIdx(iPlot));
        hold on
        scatter(COP(1, :), COP(2, :), 2, 'blue');
        scatter(footPos(1, :), footPos(2, :), 2, 'k');
        scatter(startPos(1), startPos(2), 100, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'g');
        text(startPos(1), startPos(2), 'startPos', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        scatter(stopPos(1), stopPos(2), 100, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'r');
        text(stopPos(1), stopPos(2), 'stopPos', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        scatter(landingPos(1), landingPos(2), 100, 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'k');
        text(landingPos(1), landingPos(2), 'landingPos', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');        
        xline(startPos(1), 'g');
        xline(stopPos(1), 'r');
        title(sprintf('COP path'), 'Interpreter', 'none');
        xlabel('x [m]');
        ylabel('y [m]');
        axis equal
        xlims = xlim;
        ylims = ylim;
        text(xlims(1), mean(ylims), sprintf('targetError = %.0f mm', TargetError * 1000), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
        % legend({'COP', 'Foot pos', 'jump start pos', 'jump stop pos', 'landing'});
        
        % plot foot pos
        iPlot = iPlot+1;
        subplot(nRows, nCols, plotIdx(iPlot));
        plot(Time', footPos');
        hold on
        yline(TargetError, 'r');
        yl = ylim;
        text(Time(end), max(yl), sprintf('targetError = %.0f mm', TargetError * 1000), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
        xlim([Time(1), Time(end)]);
        title(sprintf('Foot position'), 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Position [m]');

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

% report item finish
fprintf('\t\tFinished in %.3f s\n', toc(ticItem));

end