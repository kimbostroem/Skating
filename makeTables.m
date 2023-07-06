function makeTables()

fprintf('\nMaking Tables...\n');


% load current state
Measurements = loadState();

% get outDir
outDir = evalin('base', 'outDir');

%% Motor tables

MotorMetrics = Measurements.MotorMetrics;

depVars = {'PathLength', 'Fluctuation', 'Jerk', 'JerkXY'};
subjects = unique(Measurements.Subjects.Subject);

%% Clean motor table

SourceTable = MotorMetrics;
% remove empty "isValid" rows
rows = isempty(SourceTable.isValid);
SourceTable(rows, :) = [];
% remove rows where "isValid" equals zero
rows = (SourceTable.isValid == 0);
SourceTable(rows, :) = [];
SourceTable = removevars(SourceTable, {'isValid'});
MotorTable_all = SourceTable;
% remove unnecessary variables
rmVars = {'subjectCode', 'date', 'Beidbein_start', 'Beidbein_stop', 'Einbein_start', 'Einbein_stop'};
myRmVars = intersect(MotorTable_all.Properties.VariableNames, rmVars);
MotorTable_all = removevars(MotorTable_all, myRmVars);

%% Motor table in long format (pre and post values in rows)

tasks = unique(MotorTable_all.Task);
stages = unique(MotorTable_all.Stage);

SourceTable = MotorTable_all;
TargetTable = struct([]);
iRow = 1;
newDepVars = {};
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.Subject == string(subject), :);
    for iStage = 1:length(stages)
        stage = stages(iStage);
        stageTable = subjectTable(subjectTable.Stage == stage, :);
        if isempty(stageTable)
            continue
        end
        for iTask = 1:length(tasks)
            task = tasks{iTask};
            taskTable = stageTable(stageTable.Task == task, :);
            if isempty(taskTable)
                continue
            end
            initRow = table2struct(taskTable(1, :));
            variables = setdiff(fieldnames(initRow), [depVars, {'Task'}], 'stable');
            for iVar = 1:length(variables)
                variable = variables{iVar};
                TargetTable(iRow).(variable) = initRow.(variable);
            end
            for iVar = 1:length(depVars)
                depVar = depVars{iVar};
                values = taskTable.(depVar);
                newDepVar = sprintf('%s_%s', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                TargetTable(iRow).(newDepVar) = mean(values, 'omitnan');
                value = std(values, 'omitnan');
                if value == 0
                    value = NaN;
                end
                newDepVar = sprintf('%s_%s_std', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                TargetTable(iRow).(newDepVar) = value;
                newDepVar = sprintf('%s_%s_n', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                TargetTable(iRow).(newDepVar) = length(values);
            end
        end
        % increment row index
        iRow = iRow + 1;
    end
end
% update list of dependent variables
depVars = newDepVars;
% convert structure array to table
TargetTable = struct2table(TargetTable);
% clean up
TargetTable = removevars(TargetTable, {'Trial', 'Side', 'FileName'});
MotorTable_long = TargetTable;

%% Motor table in wide format (pre and post values in columns)

SourceTable = MotorTable_long;
TargetTable = struct([]);
iRow = 1;
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.Subject == string(subject), :);
    if isempty(subjectTable)
        continue
    end
    initRow = table2struct(subjectTable(1, :));
    variables = setdiff(fieldnames(initRow), [depVars, {'Stage'}], 'stable');
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = initRow.(variable);
    end
    for iVar = 1:length(depVars)
        depVar = depVars{iVar};
        % pre
        value = subjectTable.(depVar)(subjectTable.Stage == stages(1));
        if isempty(value)
            value = NaN;
        elseif iscell(value)
            value = cell2mat(value);
        end
        TargetTable(iRow).([depVar, '_pre']) = value;
        % post
        value = subjectTable.(depVar)(subjectTable.Stage == stages(2));
        if isempty(value)
            value = NaN;
        elseif iscell(value)
            value = cell2mat(value);
        end
        TargetTable(iRow).([depVar, '_post']) = value;
        % post - pre
        TargetTable(iRow).([depVar, '_diff']) = (TargetTable(iRow).([depVar, '_post']) - TargetTable(iRow).([depVar, '_pre']));
    end
    % increment row index
    iRow = iRow + 1;
end
% convert structure array to table
TargetTable = struct2table(TargetTable);
MotorTable_wide = TargetTable;


%% Cognition tables

CognitionTable_all = Measurements.CognitionData;

% delete irrelevant columns
variables_orig = CognitionTable_all.Properties.VariableNames;
variables_clean = [
    "Subject"
    "ADHS"
    "Stage"
    "Intervention"
    "Sex"
    "Age_yrs"
    "Height_cm"
    "Weight_kg"
    "Medikation"
    "AD_MW"
    "Hyp_MW"
    "D2_F__SW"
    "D2_BZO_SW"
    "D2_KL_SW"
    "Stroop_FWL_SW"
    "Stroop_FSB_SW"
    "Stroop_INT_SW"
    ];
idx = ~ismember(variables_orig, variables_clean);
CognitionTable_all(:, idx) = [];

%% Combined table long format (pre and post values in rows)

SourceMotorTable = MotorTable_long;
SourceCognitionTable = CognitionTable_all;
TargetTable = struct([]);
iRow = 1;
subjects = unique(SourceMotorTable.Subject);
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectMotorTable = SourceMotorTable(SourceMotorTable.Subject == string(subject), :);
    subjectCognitionTable = SourceCognitionTable(SourceCognitionTable.Subject == string(subject), :);
    for iStage = 1:length(stages)
        stage = stages(iStage);
        stageMotorTable = subjectMotorTable(subjectMotorTable.Stage == stage, :);
        stageCognitionTable = subjectCognitionTable(subjectCognitionTable.Stage == stage, :);
        if isempty(stageMotorTable) || isempty(stageCognitionTable)
            continue
        end
        if size(stageMotorTable, 1) > 1
            warning('Motor table contains for subject %s and stage %d more than 1 entry -> skip additional entries\n', subject, stage);
        end
        if size(stageCognitionTable, 1) > 1
            warning('Cognition table contains for subject %s and stage %d more than 1 entry -> skip additional entries\n', subject, stage);
        end
        initRow = table2struct(stageMotorTable(1, :));
        variables = stageMotorTable.Properties.VariableNames;
        for iVar = 1:length(variables)
            variable = variables{iVar};
            TargetTable(iRow).(variable) = initRow.(variable);
        end
        variables = setdiff(stageCognitionTable.Properties.VariableNames, variables, 'stable');
        depVars = union(depVars, variables, 'stable');
        for iVar = 1:length(variables)
            variable = variables{iVar};
            TargetTable(iRow).(variable) = stageCognitionTable.(variable);
        end
        % increment row index
        iRow = iRow + 1;
    end
end
% convert structure array to table
TargetTable = struct2table(TargetTable);
SkatingTable_long = TargetTable;

%% Combined table wide format (pre and post values in columns)

SourceTable = SkatingTable_long;
TargetTable = struct([]);
nSubjects = length(subjects);
iRow = 1;
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.Subject == string(subject), :);
    if isempty(subjectTable)
        continue
    end
    tic
    initRow = table2struct(subjectTable(1, :));
    variables = setdiff(fieldnames(initRow), [depVars, {'Stage'}], 'stable');
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = initRow.(variable);
    end
    for iVar = 1:length(depVars)
        depVar = depVars{iVar};
        % pre
        value = subjectTable.(depVar)(subjectTable.Stage == stages(1));
        if isempty(value)
            value = NaN;
        elseif iscell(value)
            value = cell2mat(value);
        end
        TargetTable(iRow).([depVar, '_pre']) = value;
        % post
        value = subjectTable.(depVar)(subjectTable.Stage == stages(2));
        if isempty(value)
            value = NaN;
        elseif iscell(value)
            value = cell2mat(value);
        end
        TargetTable(iRow).([depVar, '_post']) = value;
        % post - pre
        TargetTable(iRow).([depVar, '_diff']) = (TargetTable(iRow).([depVar, '_post']) - TargetTable(iRow).([depVar, '_pre']));
    end

    % report progress
    fprintf('\t-> %s (%d/%d = %.1f%% in %.3fs)\n', subject, iSubject, nSubjects, iSubject/nSubjects*100, toc);

    % increment row index
    iRow = iRow + 1;
end
% convert structure array to table
TargetTable = struct2table(TargetTable);
SkatingTable_wide = TargetTable;

%% Append tables to Measurements structure

Measurements.MotorTable_all = MotorTable_all;
Measurements.MotorTable_long = MotorTable_long;
Measurements.MotorTable_wide = MotorTable_wide;
Measurements.CognitionTable_all = CognitionTable_all;
Measurements.SkatingTable_long = SkatingTable_long;
Measurements.SkatingTable_wide = SkatingTable_wide;

%% save Measurements structure

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

%% save table to disk

% write Motor table
fprintf('Saving Motor table...\n');
saveTable(MotorTable_all, 'MotorTable_all', {'xlsx'}, outDir);

% write Cognition table
fprintf('Saving Cognition table...\n');
saveTable(CognitionTable_all, 'CognitionTable_all', {'xlsx'}, outDir);


% write Skating table in long format
fprintf('Saving Skating table in long format...\n');
saveTable(SkatingTable_long, 'SkatingTable_long', {'csv'}, outDir);

% write Skating table in wide format
fprintf('Saving Skating table in wide format...\n');
saveTable(SkatingTable_wide, 'SkatingTable_wide', {'csv'}, outDir);

end
