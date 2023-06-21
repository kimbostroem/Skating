function makeMotorMetrics()

outDir = evalin('base', 'outDir');

% load current state
Measurements = loadState();

% load table containing subjects info
Subjects = Measurements.Subjects;

nFiles = length(Measurements.MotorData);

fprintf('Create metrics...\n');
ticAll = tic;
iMeas = 1; 
MotorMetrics = table;
for iFile = 1:nFiles    
    fileName = Measurements.MotorData(iFile).fileName;

    % skip invalid files
    if isempty(fileName) || contains(fileName, 'ungueltig', 'IgnoreCase', true)
        fprintf('\t-> Skipping invalid file %s\n', fileName);
        continue
    end

    % create empty table row
    tableRow = table;

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
    tableRow.Subject = string(subject);
    tableRow.SubjectCode = string(subjectCode);

    % stage
    stages = {'I', 'II', 'III'};
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    tableRow.Stage = stage;

    % subject properties
    subjectProps = {'ADHS', 'Medication'};
    for iProp = 1:length(subjectProps)
        propName = subjectProps{iProp};
        propValue = Subjects.(propName)(subjectIdx);
        tableRow.(subjectProps{iProp}) = propValue;
    end

    % subject properties with trailing 'I', 'II', or 'III'
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(subjectIdx);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                tableRow.(subjectProps{iProp}) = propValue;
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
    tableRow.Intervention = intervention;

    % task
    task = parts{3};
    tableRow.Task = string(task);

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        side = parts{5};
        trial = parts{6};
    else
        side = parts{4};
        trial = parts{5};
    end
    tableRow.Side = string(side);

    % jump position markers
    tableRow.Beidbein_start = string(Subjects.Beidbein_start(subjectIdx));
    tableRow.Beidbein_stop = string(Subjects.Beidbein_stop(subjectIdx));
    tableRow.Einbein_start = string(Subjects.Einbein_start(subjectIdx));
    tableRow.Einbein_stop = string(Subjects.Einbein_stop(subjectIdx));

    % trial number
    tableRow.Trial = str2double(trial);

    % store file name in tableRow
    tableRow.FileName = string(fileName);

    % set flags
    tableRow.DonePlots = 0;    

    % % report progress
    % fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iFile, nFiles, iFile/nFiles*100);

    % get variables
    subjectWeight = tableRow.Weight;
    Time = Measurements.MotorData(iFile).Time;
    Force = Measurements.MotorData(iFile).Force;
    COP = Measurements.MotorData(iFile).COP;
    idxContact = Measurements.MotorData(iFile).idxContact;
    sampleRate = Measurements.MotorData(iFile).sampleRate;
    stopPos = Measurements.MotorData(iFile).stopPos;

    if isempty(sampleRate)
        fprintf('\t\tEmpty dataset - skipping\n');
        continue
    end

    dt = 1/sampleRate; % time step size [s]
    task = tableRow.Task;

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
            Measurements.MotorData(iFile).targetDist = targetDist;
            Measurements.MotorData(iFile).Jerk = Jerk;
            tableRow.PathLength = pathLength;
            tableRow.TargetError = targetError;
            tableRow.MeanJerk = meanJerk;
            tableRow.MeanJerkXY = meanJerkXY;

        case 'Sprung'
            % distance to jump stop position
            targetDist = abs(Measurements.MotorData(iFile).COP(1, :) - stopPos(1));
            [targetError, targetIdx] = min(targetDist);
            jumpStopPos = Measurements.MotorData(iFile).COP(:, targetIdx);
            Measurements.MotorData(iFile).targetDist = targetDist;
            tableRow.TargetError = targetError;
            Measurements.MotorData(iFile).jumpStopPos = jumpStopPos;

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
            Measurements.MotorData(iFile).Jerk = Jerk;
            tableRow.PathLength = pathLength;
            tableRow.TargetError = targetError;
            tableRow.MeanJerk = meanJerk;
            tableRow.MeanJerkXY = meanJerkXY;
    end

    % mark invalid datasets
    tableRow.isValid = 1;
    if pathLength == 0 % something is very wrong
        tableRow.isValid = 0;
        tableRow.PathLength = NaN;
        tableRow.TargetError = NaN;
        tableRow.MeanJerk = NaN;
        tableRow.MeanJerkXY = NaN;
    end

    % append table row to table
    MotorMetrics = [MotorMetrics; tableRow]; %#ok<AGROW>

    % increment number of processed files
    iMeas = iMeas+1;    
end

% append MotorMetrics table to Measurements structure
Measurements.MotorMetrics = MotorMetrics;

fprintf('\t\t- Saving MotorMetrics to table...\n');
outpath = fullfile(outDir, 'MotorMetrics.xlsx');
writetable(Measurements.MotorMetrics, outpath, 'WriteMode', 'replacefile');

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

fprintf('Finished creating metrics from %d datasets in %.3f s\n', iMeas, toc(ticAll));

end