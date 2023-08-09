function makeTables()

fprintf('\nMaking Tables...\n');


% load current state
Measurements = loadState();

% get outDir
outDir = evalin('base', 'outDir');

%% Motor table

MotorMetrics = Measurements.MotorMetrics;

depVars = {'PathLength', 'TargetError', 'Jerk', 'JerkXY', 'JerkZ'};

% clean table

SourceTable = MotorMetrics;
% remove empty "isValid" rows
rows = isempty(SourceTable.isValid);
SourceTable(rows, :) = [];
% remove rows where "isValid" equals zero
rows = (SourceTable.isValid == 0);
SourceTable(rows, :) = [];
SourceTable = removevars(SourceTable, {'isValid'});
MotorTable = SourceTable;
% remove unnecessary variables
rmVars = {'subjectCode', 'date', 'Beidbein_start', 'Beidbein_stop', 'Einbein_start', 'Einbein_stop'};
myRmVars = intersect(MotorTable.Properties.VariableNames, rmVars);
MotorTable = removevars(MotorTable, myRmVars);

% replace missing height at some stage with mean of heights of all stages
subjects = unique(Measurements.Subjects.Subject, 'stable');
nSubjects = length(subjects);
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    idx = find(MotorTable.Subject == string(subject));
    variable = 'Height';
    if ~ismember(variable, MotorTable.Properties.VariableNames)
        continue
    end
    myCell = MotorTable{idx, variable};
    if iscell(myCell)
        idxNaN = cellfun(@(x) any(isempty(x)), myCell);
    else
        idxNaN = arrayfun(@(x) any(isempty(x)), myCell);
    end
    if any(idxNaN)
        fillValue = mean(cell2mat(myCell), 'omitnan');
        for iIdx = 1:length(idx)
            MotorTable.(variable){idx(iIdx)} = fillValue;
        end
    end
end

% append table to Measurements structure
Measurements.MotorTable = MotorTable;

%% Motor table subject-wise

tasks = unique(MotorTable.Task, 'stable');
stages = unique(MotorTable.Stage, 'stable');
SourceTable = MotorTable;
TargetTable = struct([]);
iRow = 1;
newDepVars = {};
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectSourceTable = SourceTable(SourceTable.Subject == string(subject), :);
    for iStage = 1:length(stages)
        stage = stages(iStage);
        stageTable = subjectSourceTable(subjectSourceTable.Stage == stage, :);
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
            allVariables = setdiff(fieldnames(initRow), [depVars, {'Task'}], 'stable');
            for iVar = 1:length(allVariables)
                variable = allVariables{iVar};
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
MotorTable_subjectMean = TargetTable;

%% Cognition table

CognitionTable_orig = Measurements.CognitionData;

% delete irrelevant columns
cognitionVariables_orig = CognitionTable_orig.Properties.VariableNames;
cognitionVariables_clean = [
    "Subject"
    "Stage"
    "Intervention"
    "AD_MW"
    "Hyp_MW"
    "D2_F__SW"
    "D2_BZO_SW"
    "D2_KL_SW"
    "Stroop_FWL_SW"
    "Stroop_FSB_SW"
    "Stroop_INT_SW"
    ];
cognitionVariables = intersect(cognitionVariables_clean, cognitionVariables_orig, 'stable');
CognitionTable = CognitionTable_orig(:, cognitionVariables);

% rename variables
renameVariables_old = [
    "AD_MW"
    "Hyp_MW"
    "D2_F__SW"
    "D2_BZO_SW"
    "D2_KL_SW"
    "Stroop_FWL_SW"
    "Stroop_FSB_SW"
    "Stroop_INT_SW"
    ];
renameVariables_new = [
    "AttentionDeficit"
    "Hyperactivity"
    "D2_Error"
    "D2_Completed"
    "D2_Concentration"
    "ColorWord"
    "ColorBar"
    "Stroop"
    ];
CognitionTable = renamevars(CognitionTable, renameVariables_old, renameVariables_new);
for iVar = 1:length(renameVariables_old)
    depVars = cellstr(strrep(depVars, renameVariables_old(iVar), renameVariables_new(iVar)));
end

% append table to Measurements structure
Measurements.CognitionTable = CognitionTable;

%% Subject table

subjectVariables_orig = Measurements.Subjects.Properties.VariableNames;
subjectVariables_clean = [
    "Subject"
    "ADHS"
    "Diagnose"
    "Sex"
    "Age_yrs"
    "Medication"
    "School"
    "Grade"
    "Skating_reg"
    "Skating_add"
    "Skating_tot"
    "Skating_tot_bin"
    ];
subjectVariables = intersect(subjectVariables_clean, subjectVariables_orig, 'stable');
SubjectTable = Measurements.Subjects(:, subjectVariables);

%% SkatingTable all motor measurements

TargetTable = struct([]);
nRows = size(MotorTable, 1);
for iRow = 1:nRows
    subject = MotorTable.Subject(iRow);
    % add subject data
    subjectData = SubjectTable(SubjectTable.Subject == subject, :);
    variables = subjectData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = subjectData.(variable);
    end    
    % add cognition data
    stage = MotorTable.Stage(iRow);
    cognitionData = CognitionTable(CognitionTable.Subject == subject & CognitionTable.Stage == stage, :);
    if isempty(cognitionData)
        continue
    end
    variables = cognitionData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = cognitionData.(variable);
    end
    % add motor data
    motorData = MotorTable(iRow, :);
    variables = motorData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = motorData.(variable);
    end
end

% convert structure array to table
SkatingTable = struct2table(TargetTable);

% rename variables
SkatingTable = renamevars(SkatingTable, 'Task', 'MotorTask');

% append table to Measurements structure
Measurements.SkatingTable = SkatingTable;

%% SkatingTable subjectMean

TargetTable = struct([]);
iRow = 1;
subjects = unique(MotorTable_subjectMean.Subject, 'stable');
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectMotorTable = MotorTable_subjectMean(MotorTable_subjectMean.Subject == string(subject), :);
    subjectCognitionTable = CognitionTable(CognitionTable.Subject == string(subject), :);
    subjectSubjectTable = SubjectTable(SubjectTable.Subject == string(subject), :);
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

        % init table row with subject data
        allVariables = {};
        subjectSourceTable = table2struct(subjectSubjectTable);
        newVariables = setdiff(fieldnames(subjectSourceTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            TargetTable(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append motor data
        subjectSourceTable = table2struct(stageMotorTable);
        newVariables = setdiff(fieldnames(subjectSourceTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            TargetTable(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append cognition data
        subjectSourceTable = table2struct(stageCognitionTable);
        newVariables = setdiff(fieldnames(subjectSourceTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            TargetTable(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
            % update list of dependent variables
            depVars = union(depVars, {variable}, 'stable');
        end
        
        % increment row index
        iRow = iRow + 1;
    end
end

% convert structure array to table
SkatingTable_subjectMean = struct2table(TargetTable);

% append table to Measurements structure
Measurements.SkatingTable_subjectMean = SkatingTable_subjectMean;

%% save Measurements structure

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

%% save table to disk

% write Subjects table
fprintf('Saving Subjects table...\n');
saveTable(Measurements.Subjects, 'SubjectsTable', {'xlsx'}, outDir);

% write Motor table
fprintf('Saving Motor table...\n');
saveTable(Measurements.MotorTable, 'MotorTable', {'xlsx'}, outDir);

% write Cognition table
fprintf('Saving Cognition table...\n');
saveTable(Measurements.CognitionTable, 'CognitionTable', {'xlsx'}, outDir);

% write Skating table subjectMean
fprintf('Saving Skating table in long format...\n');
saveTable(SkatingTable_subjectMean, 'SkatingTable_subjectMean', {'csv'}, outDir);

% write Skating table 
fprintf('Saving Skating table in wide format...\n');
saveTable(SkatingTable, 'SkatingTable', {'csv'}, outDir);

fprintf('If necessary, save current state using ''saveState''\n');


end
