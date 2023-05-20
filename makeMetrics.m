function Measurements = makeMetrics(Measurements)

nMeas = length(Measurements.Info);
outDir = Measurements.outDir;

fprintf('Create metrics...\n');
ticAll = tic;
for iMeas = 1:nMeas
    fileName = Measurements.Info(iMeas).fileName;
    fprintf('\t-> %s\n', fileName);

    subjectWeight = Measurements.Info(iMeas).weight;
    Time = Measurements.Data(iMeas).Time;
    Force = Measurements.Data(iMeas).Force;
    COP = Measurements.Data(iMeas).COP;
    idxContact = Measurements.Data(iMeas).idxContact;
    sampleRate = Measurements.Data(iMeas).sampleRate;
    dt = 1/sampleRate; % time step size [s]
    nSamples = Measurements.Data(iMeas).nSamples;
    task = Measurements.Info(iMeas).task;

    % distance to beam (for task "Balance")
    COPx = COP(1, idxContact);
    COPy = COP(2, idxContact);
    coefficients = polyfit(COPx, COPy, 1);
    yBeamFcn = @(x) polyval(coefficients, x);
    beamDistFcn = @(x, y) abs(y - yBeamFcn(x));
    BeamDist = beamDistFcn(COP(1, :), COP(2, :));

    % jerk
    dJerk = diff(Force, 1, 2);
    Jerk = [dJerk, dJerk(:, end)] / (dt * subjectWeight);

    % path length
    pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));
    
    % mean beam distance
    meanBeamDist = mean(BeamDist, 'omitnan');

    % mean jerk
    meanJerk = mean(vecnorm(Jerk, 2, 1), 'omitnan');
    meanJerkXY = mean(vecnorm(Jerk(1:2, :), 2, 1), 'omitnan');


    %% Store data

    Measurements.Data(iMeas).BeamDist = BeamDist;
    Measurements.Data(iMeas).Jerk = Jerk;  
    Measurements.Info(iMeas).pathLength= pathLength;
    Measurements.Info(iMeas).meanBeamDist = meanBeamDist;
    Measurements.Info(iMeas).meanJerk = meanJerk;
    Measurements.Info(iMeas).meanJerkXY = meanJerkXY;

    %% Plot data

    % setup figure
    nRows = 3;
    nCols = 1;
    iPlot = 0;
    fig = setupFigure(nCols*400, nRows*200, fileName);

    % plot COP path with beam
    iPlot = iPlot+1;
    subplot(nRows, nCols, iPlot);
    hold on
    xFit = linspace(min(COP(1, :)), max(COP(1, :)), nSamples);
    yFit = yBeamFcn(xFit);
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
    plot(Time', beamDistFcn(COP(1, :), COP(2, :))');
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

    % save figure
    ftypes = {'png', 'pdf', 'fig'};
    for iType = 1:length(ftypes)
        ftype = ftypes{iType};
        myOutDir = fullfile(outDir, 'metric', ftype);
        if ~isfolder(myOutDir)
            mkdir(myOutDir);
        end
        outpath = fullfile(myOutDir, fileName);
        saveFigure(fig, outpath, ftype);
    end
    close(fig);
end

%% Saving table
validIdx = [Measurements.Info.isForce];
MeasurementTable = struct2table(Measurements.Info(validIdx));
outpath = fullfile(outDir, 'Measurements.xlsx');
writetable(MeasurementTable, outpath, 'WriteMode', 'replacefile');

sprintf('Finished in %f s\n\n', toc(ticAll));

end