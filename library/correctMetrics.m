function correctMetrics

% load current state
loadState;

nMeas = length(Measurements.Observations); %#ok<NODEF>
outDir = '.';

fprintf('Correct metrics...\n');
for iMeas = 1:nMeas
    fileName = Measurements.Observations(iMeas).fileName;

    % if metrics not yet calculated -> skip
    if ~Measurements.Observations(iMeas).doneMetrics
        Measurements.Observations(iMeas).isValid = 0;
        continue
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iMeas, nMeas, iMeas/nMeas*100);

    % get variables
    task = Measurements.Observations(iMeas).task;
    COP = Measurements.Data(iMeas).COP;
    idxContact = Measurements.Data(iMeas).idxContact;
    stopPos = Measurements.Data(iMeas).stopPos;
    pathLength = Measurements.Observations(iMeas).pathLength;

    % mark invalid datasets
    Measurements.Observations(iMeas).isValid = 1;
    if pathLength <= 0 || isnan(pathLength) % something is very wrong
        Measurements.Observations(iMeas).isValid = 0;
        Measurements.Observations(iMeas).pathLength = NaN;
        Measurements.Observations(iMeas).targetError = NaN;
        Measurements.Observations(iMeas).meanJerk = NaN;
        Measurements.Observations(iMeas).meanJerkXY = NaN;
    end

    % get target error
    switch task
        case 'Balance'
            % distance to beam (for task "Balance")
            coefficients = polyfit(COP(1, idxContact), COP(2, idxContact), 1);
            targetFcn = @(x, y) polyval(coefficients, x);
            targetDistFcn = @(x, y) abs(y - targetFcn(x, y));
            targetDist = targetDistFcn(COP(1, :), COP(2, :));
            targetError = mean(targetDist, 'omitnan');
        case 'Einbein'
            meanCOP = mean(COP, 2, 'omitnan');
            targetDistFcn = @(x, y) vecnorm([x; y] - meanCOP);
            targetDist = targetDistFcn(COP(1, :), COP(2, :));
            targetError = mean(targetDist, 'omitnan');
        case 'Sprung'
            targetDist = abs(COP(1, :) - stopPos(1));
            [targetError, ~] = min(targetDist);            
    end
    % store corrected variables
    Measurements.Data(iMeas).targetDist = targetDist;
    Measurements.Observations(iMeas).targetError = targetError;
end

saveState;

end