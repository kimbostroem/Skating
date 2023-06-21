function makeMetrics

paramDir = evalin('base', 'paramDir');

% load current state
loadState;

% load table containing subjects info
Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'));

nFiles = length(Measurements.Data); %#ok<NODEF>
outDir = Measurements.outDir;

fprintf('Create metrics...\n');
ticAll = tic;
iMeas = 1; 
Measurements.Findings = struct([]);
for iFile = 1:nFiles    
    fileName = Measurements.Data(iFile).fileName;

    % skip invalid files
    if isempty(fileName) || contains(fileName, 'ungueltig', 'IgnoreCase', true)
        fprintf('\t-> Skipping invalid file %s\n', fileName);
        continue
    end

    % split file name at underscores
    parts = strsplit(fileName, '_');

    % subject identity and code
    subject = parts{1};
    subjectCode = [parts{1}, '_', parts{2}];
    subjectCodes = [Subjects.Code_I, Subjects.Code_II, Subjects.Code_III];
    [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
    if isempty(subjectIdx)
        warning('Subject code %s not found -> skipping', subjectCode);
        continue
    end
    Measurements.Findings(iMeas).subject = string(subject);
    Measurements.Findings(iMeas).subjectCode = string(subjectCode);

    % stage
    stages = {'I', 'II', 'III'};
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    Measurements.Findings(iMeas).stage = stage;

    % subject properties
    subjectProps = {'ADHS', 'Medication'};
    for iProp = 1:length(subjectProps)
        propName = subjectProps{iProp};
        propValue = Subjects.(propName)(subjectIdx);
        Measurements.Findings(iMeas).(lower(subjectProps{iProp})) = propValue;
    end


    % subject properties with trailing 'I', 'II', or 'III'
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(subjectIdx);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                Measurements.Findings(iMeas).(lower(subjectProps{iProp})) = propValue;
                break
            end
            myStage = myStage-1;
        end
    end

    % intervention
    [~, subjectStages] = find(contains(subjectCodes, subject));
    switch stageStr
        case 'I'
            intervention = 0;
        case 'II'
            if max(subjectStages) == 3 % subject has been tested at three stages
                intervention = 0;
            else
                intervention = 1;
            end
        case 'III'
            intervention = 1;
    end
    Measurements.Findings(iMeas).intervention = intervention;

    % task
    task = parts{3};
    Measurements.Findings(iMeas).task = string(task);

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        side = parts{5};
        trial = parts{6};
    else
        side = parts{4};
        trial = parts{5};
    end
    Measurements.Findings(iMeas).side = string(side);

    % jump position markers
    Measurements.Findings(iMeas).Beidbein_start = string(Subjects.Beidbein_start(subjectIdx));
    Measurements.Findings(iMeas).Beidbein_stop = string(Subjects.Beidbein_stop(subjectIdx));
    Measurements.Findings(iMeas).Einbein_start = string(Subjects.Einbein_start(subjectIdx));
    Measurements.Findings(iMeas).Einbein_stop = string(Subjects.Einbein_stop(subjectIdx));

    % trial number
    Measurements.Findings(iMeas).trial = str2double(trial);

    % store file name in Measurements.Findings(iMeas)
    Measurements.Findings(iMeas).fileName = string(fileName);

    % set flags
    Measurements.Findings(iMeas).donePlots = 0;    

    % report progress
    fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iFile, nFiles, iFile/nFiles*100);

    % get variables
    subjectWeight = Measurements.Findings(iMeas).weight;
    Time = Measurements.Data(iFile).Time;
    Force = Measurements.Data(iFile).Force;
    COP = Measurements.Data(iFile).COP;
    idxContact = Measurements.Data(iFile).idxContact;
    sampleRate = Measurements.Data(iFile).sampleRate;
    stopPos = Measurements.Data(iFile).stopPos;

    if isempty(sampleRate)
        fprintf('\t\tEmpty dataset - skipping\n');
        continue
    end

    dt = 1/sampleRate; % time step size [s]
    task = Measurements.Findings(iMeas).task;

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
            Measurements.Data(iFile).targetDist = targetDist;
            Measurements.Data(iFile).Jerk = Jerk;
            Measurements.Findings(iMeas).pathLength = pathLength;
            Measurements.Findings(iMeas).targetError = targetError;
            Measurements.Findings(iMeas).meanJerk = meanJerk;
            Measurements.Findings(iMeas).meanJerkXY = meanJerkXY;

        case 'Sprung'
            % distance to jump stop position
            targetDist = abs(Measurements.Data(iFile).COP(1, :) - stopPos(1));
            [targetError, targetIdx] = min(targetDist);
            jumpStopPos = Measurements.Data(iFile).COP(:, targetIdx);
            Measurements.Data(iFile).targetDist = targetDist;
            Measurements.Findings(iMeas).targetError = targetError;
            Measurements.Data(iFile).jumpStopPos = jumpStopPos;

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
            Measurements.Data(iFile).Jerk = Jerk;
            Measurements.Findings(iMeas).pathLength = pathLength;
            Measurements.Findings(iMeas).targetError = targetError;
            Measurements.Findings(iMeas).meanJerk = meanJerk;
            Measurements.Findings(iMeas).meanJerkXY = meanJerkXY;
    end

    % mark invalid datasets
    Measurements.Findings(iMeas).isValid = 1;
    if pathLength == 0 % something is very wrong
        Measurements.Findings(iMeas).isValid = 0;
        Measurements.Findings(iMeas).pathLength = NaN;
        Measurements.Findings(iMeas).targetError = NaN;
        Measurements.Findings(iMeas).meanJerk = NaN;
        Measurements.Findings(iMeas).meanJerkXY = NaN;
    end

    % increment number of processed files
    iMeas = iMeas+1;    
end

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

fprintf('\t\t- Saving Findings to table...\n');
MeasurementTable = struct2table(Measurements.Findings);
outpath = fullfile(outDir, 'Findings.xlsx');
writetable(MeasurementTable, outpath, 'WriteMode', 'replacefile');


% save current state
saveState;

fprintf('Finished creating metrics from %d datasets in %.3f s\n\n', iMeas, toc(ticAll));

end