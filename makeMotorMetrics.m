function makeMotorMetrics()

% load current state
Measurements = loadState();

% get outDir
outDir = evalin('base', 'outDir');

% load table containing subjects info
Subjects = Measurements.Subjects;

nFiles = length(Measurements.MotorData);

fprintf('Create metrics...\n');
ticAll = tic;
item = 1; 
MotorMetrics = struct([]);
for iFile = 1:nFiles    
    fileName = Measurements.MotorData(iFile).fileName;

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
    MotorMetrics(item).Subject = string(subject);
    MotorMetrics(item).SubjectCode = string(subjectCode);

    % stage
    stages = {'I', 'II', 'III'};
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    MotorMetrics(item).Stage = stage;

    % subject properties
    subjectProps = {'ADHS', 'Medication'};
    for iProp = 1:length(subjectProps)
        propName = subjectProps{iProp};
        propValue = Subjects.(propName)(subjectIdx);
        MotorMetrics(item).(subjectProps{iProp}) = propValue;
    end

    % subject properties with trailing 'I', 'II', or 'III'
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(subjectIdx);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                MotorMetrics(item).(subjectProps{iProp}) = propValue;
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
    MotorMetrics(item).Intervention = intervention;

    % task
    task = parts{3};
    MotorMetrics(item).Task = string(task);

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        side = parts{5};
        trial = parts{6};
    else
        side = parts{4};
        trial = parts{5};
    end
    MotorMetrics(item).Side = string(side);

    % jump position markers
    MotorMetrics(item).Beidbein_start = string(Subjects.Beidbein_start(subjectIdx));
    MotorMetrics(item).Beidbein_stop = string(Subjects.Beidbein_stop(subjectIdx));
    MotorMetrics(item).Einbein_start = string(Subjects.Einbein_start(subjectIdx));
    MotorMetrics(item).Einbein_stop = string(Subjects.Einbein_stop(subjectIdx));

    % trial number
    MotorMetrics(item).Trial = str2double(trial);

    % store file name in MotorMetrics(item)
    MotorMetrics(item).FileName = string(fileName);

    % set flags
    MotorMetrics(item).DonePlots = 0;    

    % % report progress
    % fprintf('\t-> %s (%d/%d = %.0f%%)\n', fileName, iFile, nFiles, iFile/nFiles*100);

    % get variables
    subjectWeight = MotorMetrics(item).Weight;
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
    % task = MotorMetrics(item).Task;

    switch task

        case {'Balance', 'Einbein'}

            if strcmp(task, 'Balance')
                % Calc distance to beam for task "Balance". The beam is
                % parallel to the x-axis, therefore fit a polynom of 0th
                % order (constant function)
                polyOrder = 0;
                coefficients = polyfit(COP(1, idxContact), COP(2, idxContact), polyOrder);
                targetFcn = @(x, y) polyval(coefficients, x);
                deviationFcn = @(x, y) abs(y - targetFcn(x, y));                
            elseif strcmp(task, 'Einbein')
                meanCOP = mean(COP, 2, 'omitnan');
                deviationFcn = @(x, y) vecnorm([x; y] - meanCOP);
            end
            deviation = deviationFcn(COP(1, :), COP(2, :));
            fluctuation = mean(deviation, 'omitnan');

        case 'Sprung'
            % distance to jump stop position
            deviation = abs(Measurements.MotorData(iFile).COP(1, :) - stopPos(1));
            [fluctuation, targetIdx] = min(deviation);
            jumpStopPos = Measurements.MotorData(iFile).COP(:, targetIdx);
            Measurements.MotorData(iFile).jumpStopPos = jumpStopPos;
    end

    % jerk
    dJerk = diff(Force, 1, 2);
    Jerk = [dJerk, dJerk(:, end)] / (dt * subjectWeight);
    Jerk(:, ~idxContact) = 0; % remove jerk around gaps

    % path length
    pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));

    fluctuationName = 'Fluctuation';

    % mean jerk
    meanJerk = mean(vecnorm(Jerk, 2, 1), 'omitnan');
    meanJerkName = 'Jerk';
    meanJerkXY = mean(vecnorm(Jerk(1:2, :), 2, 1), 'omitnan');
    meanJerkXYName = 'JerkXY';

    % store metrics in structure
    Measurements.MotorData(iFile).Deviation = deviation;
    Measurements.MotorData(iFile).Jerk = Jerk;
    MotorMetrics(item).PathLength = pathLength;
    MotorMetrics(item).(fluctuationName) = fluctuation;
    MotorMetrics(item).(meanJerkName) = meanJerk;
    MotorMetrics(item).(meanJerkXYName) = meanJerkXY;

    % mark invalid datasets
    MotorMetrics(item).isValid = 1;
    if pathLength == 0 % something is very wrong
        MotorMetrics(item).isValid = 0;
        MotorMetrics(item).PathLength = NaN;
        MotorMetrics(item).(fluctuationName) = NaN;
        MotorMetrics(item).(meanJerkName) = NaN;
        MotorMetrics(item).(meanJerkXYName) = NaN;
    end

    % increment number of processed files
    item = item+1;    
end

% append MotorMetrics table to Measurements structure
Measurements.MotorMetrics = struct2table(MotorMetrics);

fprintf('\t\t- Saving MotorMetrics to table...\n');
outpath = fullfile(outDir, 'MotorMetrics.xlsx');
writetable(Measurements.MotorMetrics, outpath, 'WriteMode', 'replacefile');

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

fprintf('Finished creating metrics from %d datasets in %.3f s\n', item, toc(ticAll));

end