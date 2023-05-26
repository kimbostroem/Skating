function makeMetrics

% load Measurements structure
load('Measurements.mat'); %#ok<LOAD>

nMeas = length(Measurements.Observations); %#ok<NODEF>
outDir = Measurements.outDir;

fprintf('Create metrics...\n');
ticAll = tic;
nProc = 0; % init number of processed files
for iMeas = 1:nMeas
    ticItem = tic;
    fileName = Measurements.Observations(iMeas).fileName;

    if Measurements.Observations(iMeas).doneMetrics
        continue
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);

    subjectWeight = Measurements.Observations(iMeas).weight;
    Time = Measurements.Data(iMeas).Time;
    Force = Measurements.Data(iMeas).Force;
    COP = Measurements.Data(iMeas).COP;
    idxContact = Measurements.Data(iMeas).idxContact;
    sampleRate = Measurements.Data(iMeas).sampleRate;

    if isempty(sampleRate)
        fprintf('\t\tEmpty dataset - skipping\n');
        continue
    end

    dt = 1/sampleRate; % time step size [s]
    task = Measurements.Observations(iMeas).task;

    switch task

        case {'Balance', 'Einbein'}

            if strcmp(task, 'Balance')
                % distance to beam (for task "Balance")
                coefficients = polyfit(COP(1, idxContact), COP(2, idxContact), 1);
                targetFcn = @(x, y) polyval(coefficients, x);
                targetDistFcn = @(x, y) abs(y - targetFcn(x, y));
            elseif strcmp(task, 'Einbein')
                meanCOP = mean(COP, 2, 'omitnan');
                targetDistFcn = @(x, y) vecnorm([x; y] - meanCOP);
            end
            targetDist = targetDistFcn(COP(1, :), COP(2, :));

            % jerk
            dJerk = diff(Force, 1, 2);
            Jerk = [dJerk, dJerk(:, end)] / (dt * subjectWeight);
            Jerk(:, ~idxContact) = 0; % remove jerk around gaps

            % path length
            pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));

            % target error
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
            Measurements.Observations(iMeas).targetError = targetError;
            Measurements.Observations(iMeas).meanJerk = meanJerk;
            Measurements.Observations(iMeas).meanJerkXY = meanJerkXY;
    end

    % mark invalid datasets
    Measurements.Observations(iMeas).isValid = 1;
    if pathLength == 0 % something is very wrong
        Measurements.Observations(iMeas).isValid = 0;
        Measurements.Observations(iMeas).pathLength = NaN;
        Measurements.Observations(iMeas).targetError = NaN;
        Measurements.Observations(iMeas).meanJerk = NaN;
        Measurements.Observations(iMeas).meanJerkXY = NaN;
    end

    % increment number of processed files
    nProc = nProc+1;

    % set flag
    Measurements.Observations(iMeas).doneMetrics = 1;

    % export Measurements structure to base workspace
    fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
    assignin('base', 'Measurements', Measurements);

    fprintf('\t\t- Saving Observations to table...\n');
    MeasurementTable = struct2table(Measurements.Observations);
    outpath = fullfile(outDir, 'Observations.xlsx');
    writetable(MeasurementTable, outpath, 'WriteMode', 'replacefile');

    fprintf('\t\tFinished in %.3f s\n', toc(ticItem));
end

% save Measurements structure to MAT file
fprintf('\t\t- Saving Measurements structure to MAT file...\n');
save(fullfile(outDir, 'Measurements.mat'), 'Measurements');

fprintf('Finished creating metrics from %d datasets in %.3f s\n\n', nProc, toc(ticAll));

end