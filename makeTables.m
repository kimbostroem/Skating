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
SkatingTable = struct([]);
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
                SkatingTable(iRow).(variable) = initRow.(variable);
            end
            for iVar = 1:length(depVars)
                depVar = depVars{iVar};
                values = taskTable.(depVar);
                newDepVar = sprintf('%s_%s', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                SkatingTable(iRow).(newDepVar) = mean(values, 'omitnan');
                value = std(values, 'omitnan');
                if value == 0
                    value = NaN;
                end
                newDepVar = sprintf('%s_%s_std', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                SkatingTable(iRow).(newDepVar) = value;
                newDepVar = sprintf('%s_%s_n', task, depVar);
                newDepVars = union(newDepVars, {newDepVar}, 'stable');
                SkatingTable(iRow).(newDepVar) = length(values);
            end
        end
        % increment row index
        iRow = iRow + 1;
    end
end
% update list of dependent variables
depVars = newDepVars;
% convert structure array to table
SkatingTable = struct2table(SkatingTable);
% clean up
SkatingTable = removevars(SkatingTable, {'Trial', 'Side', 'FileName'});
MotorTable_subjectMean = SkatingTable;

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
CognitionTableClean = CognitionTable_orig(:, cognitionVariables);

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
CognitionTableClean = renamevars(CognitionTableClean, renameVariables_old, renameVariables_new);
for iVar = 1:length(renameVariables_old)
    depVars = cellstr(strrep(depVars, renameVariables_old(iVar), renameVariables_new(iVar)));
end

%% Subject table

subjectVariables_orig = Measurements.Subjects.Properties.VariableNames;
subjectVariables_clean = [
    "Subject"
    "Group"
    "Skating"
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

%% CognitionTable

CognitionTable = struct([]);
nRows = size(CognitionTableClean, 1);
for iRow = 1:nRows
    subject = CognitionTableClean.Subject(iRow);
    % add subject data
    subjectData = SubjectTable(SubjectTable.Subject == subject, :);
    variables = subjectData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        CognitionTable(iRow).(variable) = subjectData.(variable);
    end    
    % add cognition data
    cognitionData = CognitionTableClean(iRow, :);
    variables = cognitionData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        CognitionTable(iRow).(variable) = cognitionData.(variable);
    end
end

% convert structure array to table
CognitionTable = struct2table(CognitionTable);
% append table to Measurements structure
Measurements.CognitionTable = CognitionTable;

%% SkatingTable all motor measurements

SkatingTable = struct([]);
nRows = size(MotorTable, 1);
for iRow = 1:nRows
    subject = MotorTable.Subject(iRow);
    % add subject data
    subjectData = SubjectTable(SubjectTable.Subject == subject, :);
    variables = subjectData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        SkatingTable(iRow).(variable) = subjectData.(variable);
    end    
    % add cognition data
    stage = MotorTable.Stage(iRow);
    cognitionData = CognitionTableClean(CognitionTableClean.Subject == subject & CognitionTableClean.Stage == stage, :);
    if isempty(cognitionData)
        continue
    end
    variables = cognitionData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        SkatingTable(iRow).(variable) = cognitionData.(variable);
    end
    % add motor data
    motorData = MotorTable(iRow, :);
    variables = motorData.Properties.VariableNames;
    for iVar = 1:length(variables)
        variable = variables{iVar};
        SkatingTable(iRow).(variable) = motorData.(variable);
    end
end

% convert structure array to table
SkatingTable = struct2table(SkatingTable);
% rename variables
SkatingTable = renamevars(SkatingTable, 'Task', 'MotorTask');
% append table to Measurements structure
Measurements.SkatingTable = SkatingTable;

%% SkatingTable subjectMean

SkatingTable_subjectMean = struct([]);
iRow = 1;
subjects = unique(MotorTable_subjectMean.Subject, 'stable');
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectMotorTable = MotorTable_subjectMean(MotorTable_subjectMean.Subject == string(subject), :);
    subjectCognitionTable = CognitionTableClean(CognitionTableClean.Subject == string(subject), :);
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
            SkatingTable_subjectMean(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append motor data
        subjectSourceTable = table2struct(stageMotorTable);
        newVariables = setdiff(fieldnames(subjectSourceTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            SkatingTable_subjectMean(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append cognition data
        subjectSourceTable = table2struct(stageCognitionTable);
        newVariables = setdiff(fieldnames(subjectSourceTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            SkatingTable_subjectMean(iRow).(variable) = subjectSourceTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
            % update list of dependent variables
            depVars = union(depVars, {variable}, 'stable');
        end
        
        % increment row index
        iRow = iRow + 1;
    end
end

% convert structure array to table
SkatingTable_subjectMean = struct2table(SkatingTable_subjectMean);

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
saveTable(Measurements.MotorTable, 'MotorTable', {'csv'}, outDir);

% write Cognition table
fprintf('Saving Cognition table...\n');
saveTable(Measurements.CognitionTable, 'CognitionTable', {'csv'}, outDir);

% write Skating table 
fprintf('Saving Skating table...\n');
saveTable(SkatingTable, 'SkatingTable', {'csv'}, outDir);

% write Skating table subjectMean
fprintf('Saving Skating table averaged per subject...\n');
saveTable(SkatingTable_subjectMean, 'SkatingTable_subjectMean', {'csv'}, outDir);

fprintf('If necessary, save current state using ''saveState''\n');


end
