function Measurements = makeData(Measurements)

dataDir = Measurements.dataDir;
outDir = Measurements.outDir;
nMeas = length(Measurements.Info);

fprintf('Extract measurement data...\n');
ticAll = tic;
Measurements.Data(nMeas, 1) = struct;
for iMeas = 1:nMeas
    fileName = Measurements.Info(iMeas).fileName;
    fpath = fullfile(dataDir, fileName);
    fprintf('\t-> %s\n', fileName);

    Measurements.Info(iMeas).isForce = 1;
    tmp = load(fpath);
    fields = fieldnames(tmp);
    Data = tmp.(fields{1});
    lengthScale = 0.001; % mm -> m
    [Forces, sampleRate, COPs, ~, location, Analog] = kraftAusQualisys(Data, lengthScale);
    PowerLineHum = 50;
    nSamples = length(Forces);
    if nSamples < 2^nextpow2(PowerLineHum)/2
        fprintf('\t\tNot enough data (%d samples) -> skipping\n', nSamples);
        Measurements.Info(iMeas).isForce = 0;
        continue
    end
    Forces = -Forces; % ground reaction forces are inverse of plate forces
    nPlates = size(Forces, 1);
    % nSamples = size(Forces, 3);
    dt = 1/sampleRate; % time step size [s]
    LoadThresh = 20; % threshold below which forces are set to zero [N]
    InitForce = 200; % force threshold to cross to start trial [N]
    InitMarg = 0.5; % additional time after crossing InitForce to start trial [s]
    HumPeriod = sampleRate/50;
    % CutoffFrequency = 20;
    for iPlate = 1:nPlates
        Forces(iPlate, :, :) = periodicMedianFilter(squeeze(Forces(iPlate, :, :)), HumPeriod);
    end
    [COPs, ~, Forces, ~] = getCOPfromAnalog(Analog, Forces, [], location, sampleRate, [], [], [], COPs);
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
    Time = (0:nSamples-1)/sampleRate;

    % remove non-contact phases
    idxContact = ~any(COP, 1);
    COP(:, idxContact) = NaN;
    idxContact = ~any(isnan(COP), 1);

    %% Store data
    Measurements.Data(iMeas).fileName = fileName;
    Measurements.Data(iMeas).nSamples = nSamples;
    Measurements.Data(iMeas).sampleRate = sampleRate;
    Measurements.Data(iMeas).Time = Time;
    Measurements.Data(iMeas).Force = Force;
    Measurements.Data(iMeas).COP = COP;
    Measurements.Data(iMeas).idxContact = idxContact;    

    %% Plot data

    % setup figure
    nRows = 3;
    nCols = 1;
    iPlot = 0;
    fig = setupFigure(nCols*400, nRows*200, fileName);

    % plot COP path
    iPlot = iPlot+1;
    subplot(nRows, nCols, iPlot);
    scatter(COP(1, :), COP(2, :), 2, 'blue');
    title(sprintf('COP path'), 'Interpreter', 'none');
    xlabel('x [m]');
    ylabel('y [m]');
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

    % plot GRF
    iPlot = iPlot+1;
    subplot(nRows, nCols, iPlot);
    hold on
    plot(Time', Force');
    title(sprintf('Ground reaction force'), 'Interpreter', 'none');
    legend({'x', 'y', 'z'});
    xlabel('Time [s]');
    ylabel('Force [N]');
    xlim([Time(1), Time(end)]);
    
    % figure title
    sgtitle(sprintf('%s', fileName), 'Interpreter', 'none');

    % save figure
    ftypes = {'png', 'pdf', 'fig'};
    for iType = 1:length(ftypes)
        ftype = ftypes{iType};
        myOutDir = fullfile(outDir, 'data', ftype);
        if ~isfolder(myOutDir)
            mkdir(myOutDir);
        end
        outpath = fullfile(myOutDir, fileName);
        saveFigure(fig, outpath, ftype);
    end
    close(fig);
end
sprintf('Finished in %f s\n\n', toc(ticAll));


end