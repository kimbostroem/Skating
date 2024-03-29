function makeMotorMetrics()

fprintf('\nMaking Motor Metrics...\n');

% threshold factor to remove outliers (0 = disable, 3 = Matlab default)
outlierThresholdFactor = 4;

% load current state
Measurements = loadState();

% load table containing subjects info
SubjectsTable = Measurements.Subjects;

nFiles = length(Measurements.MotorData);

fprintf('Create metrics...\n');
ticAll = tic;
item = 1; 
MotorMetrics = struct([]);
estimatedBeam = 0;
for iFile = 1:nFiles    
    fname = Measurements.MotorData(iFile).fileName;

    % skip invalid files
    if isempty(fname) || contains(fname, 'ungueltig', 'IgnoreCase', true)
        fprintf('\t-> Skipping invalid file %s\n', fname);
        continue
    end    

    tic

    % split file name at underscores
    parts = strsplit(fname, '_');

    % subject identity and code
    subject = parts{1};
    subjectCode = [parts{1}, '_', parts{2}];
    subjectCodes = [SubjectsTable.Code_I, SubjectsTable.Code_II, SubjectsTable.Code_III];
    [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
    if isempty(subjectIdx)
        continue
    end
    MotorMetrics(item).Subject = string(subject);
    MotorMetrics(item).SubjectCode = string(subjectCode);

    % stage
    stages = {'I', 'II', 'III'};
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    MotorMetrics(item).Stage = stage;

    % subject properties with trailing 'I', 'II', or 'III'
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = SubjectsTable.(propName)(subjectIdx);
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
            if max(subjectStages) == 3 % subject has been tested at three stages, with intervention in 3rd stage
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
    MotorMetrics(item).Beidbein_start = string(SubjectsTable.Beidbein_start(subjectIdx));
    MotorMetrics(item).Beidbein_stop = string(SubjectsTable.Beidbein_stop(subjectIdx));
    MotorMetrics(item).Einbein_start = string(SubjectsTable.Einbein_start(subjectIdx));
    MotorMetrics(item).Einbein_stop = string(SubjectsTable.Einbein_stop(subjectIdx));

    % trial number
    MotorMetrics(item).Trial = str2double(trial);

    % store file name in MotorMetrics(item)
    MotorMetrics(item).FileName = string(fname);

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

        case 'Balance'

                % Calc distance to beam for task "Balance". The beam is
                % parallel to the x-axis, therefore fit a polynom of 0th
                % order (constant function)
                tic
                COPx = COP(1, :);
                COPy = COP(2, :);
                
                % deviation = sum(vecnorm(diff(COPy, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));
                
                % estimate beam function for calculation of distance
                if ~estimatedBeam
                    polyOrder = 0; % 0 = fixed distance from x-axis, 1 = line with slope
                    coefficients = polyfit(COPx(idxContact), COPy(idxContact), polyOrder);
                    beamYFcn = @(x) polyval(coefficients, x);
                    % method 1 using fminsearch (faster)
                    options = optimset('Display', 'off');
                    ppdistanceFcn = @(x, y, x_) vecnorm([x; y] - [x_; beamYFcn(x_)]);
                    estimatedBeam = 1; % set flag that estimation is done
                    fprintf('Succesfully estimated beam distance at %fm\n', coefficients);
                end

                nSamples = size(COPx, 2);
                deviation = nan(1, nSamples);
                exitflag = 1;
                for iSample = 1:nSamples
                    if ~idxContact(iSample)
                        continue
                    end
                    minFcn = @(x_) ppdistanceFcn(COPx(iSample), COPy(iSample), x_);
                    [~, deviation(iSample), exitflag] = fminsearch(minFcn, COPx(iSample), options);
                     if exitflag <= 0
                         break
                     end
                end
                if exitflag <= 0
                    % method 2 using min (slower)
                    BeamX = linspace(min(COPx), max(COPx), nSamples);
                    BeamY = beamYFcn(BeamX);
                    distanceFcn = @(x, y) min(vecnorm([x; y] - [BeamX; BeamY]));
                    deviation = arrayfun(distanceFcn, COPx, COPy);
                end
                targetError = mean(deviation, 'omitnan');

        case 'Einbein'
                meanCOP = mean(COP, 2, 'omitnan');
                deviationFcn = @(x, y) vecnorm([x; y] - meanCOP);
                deviation = deviationFcn(COP(1, :), COP(2, :));
                targetError = mean(deviation, 'omitnan');
            
        case 'Sprung'
            % get jump landing point
            landingPos = getLandingPos(Measurements.MotorData(iFile).COP);
            if isnan(landingPos)
                fprintf('\t\tNo landing position detected -> set targetError to NaN\n');
            end
            % distance to jump stop position
            targetError = abs(landingPos(1) - stopPos(1));
            Measurements.MotorData(iFile).landingPos = landingPos;
    end

    % jerk
    dForce = diff(Force, 1, 2);
    Jerk = [dForce, dForce(:, end)] / (dt * subjectWeight);
    Jerk(:, ~idxContact) = 0; % remove jerk around gaps

    % path length
    pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));

    targetErrorName = 'TargetError';

    % mean jerk
    meanJerk = mean(vecnorm(Jerk, 2, 1), 'omitnan');
    meanJerkName = 'Jerk';
    meanJerkXY = mean(vecnorm(Jerk(1:2, :), 2, 1), 'omitnan');
    meanJerkXYName = 'JerkXY';
    meanJerkZ = mean(vecnorm(Jerk(3, :), 2, 1), 'omitnan');
    meanJerkZName = 'JerkZ';

    % store metrics in structure
    Measurements.MotorData(iFile).Deviation = deviation;
    Measurements.MotorData(iFile).Jerk = Jerk;
    MotorMetrics(item).PathLength = pathLength;
    MotorMetrics(item).(targetErrorName) = targetError;
    MotorMetrics(item).(meanJerkName) = meanJerk;
    MotorMetrics(item).(meanJerkXYName) = meanJerkXY;
    MotorMetrics(item).(meanJerkZName) = meanJerkZ;

    % mark invalid datasets
    MotorMetrics(item).isValid = 1;
    if pathLength == 0 % something is very wrong
        MotorMetrics(item).isValid = 0;
        MotorMetrics(item).PathLength = NaN;
        MotorMetrics(item).(targetErrorName) = NaN;
        MotorMetrics(item).(meanJerkName) = NaN;
        MotorMetrics(item).(meanJerkXYName) = NaN;
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.1f%% in %.3fs)\n', fname, iFile, nFiles, iFile/nFiles*100, toc);

    % increment number of processed files
    item = item+1;    
end

% remove subject code, as it is already encoded in subject ID and stage
MotorMetrics = rmfield(MotorMetrics, 'SubjectCode');

% convert structure array to table
MotorMetrics = struct2table(MotorMetrics);

% remove outliers
if outlierThresholdFactor > 0
    variables = {'PathLength', 'TargetError', 'Jerk', 'JerkXY', 'JerkZ'};
    tasks = unique(MotorMetrics.Task);
    for iVar = 1:length(variables)
        variable = variables{iVar};
        for iTask = 1:length(tasks)
            rows = (Measurements.MotorMetrics.Task == task);
            Data = MotorMetrics.(variable)(rows);
            idx = isoutlier(Data, 'ThresholdFactor', outlierThresholdFactor);
            Data(idx) = NaN;
            MotorMetrics.(variable)(rows) = Data;
            fprintf('Removed %d outliers from variable ''%s'' in task ''%s''\n', sum(idx), variable, task);
        end
    end
end

% append MotorMetrics table to Measurements structure
Measurements.MotorMetrics = MotorMetrics;

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

fprintf('Finished creating metrics from %d datasets in %.3f s\n', item, toc(ticAll));
fprintf('If necessary, save current state using ''saveState''\n');

end