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

% remove empty "isValid" rows
rows = isempty(MotorMetrics.isValid);
MotorMetrics(rows, :) = [];
% remove rows where "isValid" equals zero
rows = (MotorMetrics.isValid == 0);
MotorMetrics(rows, :) = [];
MotorMetrics = removevars(MotorMetrics, {'isValid'});
MotorTable = MotorMetrics;
% remove unnecessary variables
rmVars = {'subjectCode', 'date', 'Beidbein_start', 'Beidbein_stop', 'Einbein_start', 'Einbein_stop'};
myRmVars = intersect(MotorTable.Properties.VariableNames, rmVars);
MotorTable = removevars(MotorTable, myRmVars);

subjects = unique(Measurements.Subjects.Subject, 'stable');
nSubjects = length(subjects);
for iSubject = 1:nSubjects
    subject = string(subjects{iSubject});

    % replace missing height at some stage with mean of heights of all stages
    idxSubject = find(MotorTable.Subject == subject);
    variable = 'Height';
    if ~ismember(variable, MotorTable.Properties.VariableNames)
        continue
    end
    myCell = MotorTable{idxSubject, variable};
    if iscell(myCell)
        idxNaN = cellfun(@(x) any(isempty(x)), myCell);
    else
        idxNaN = arrayfun(@(x) any(isempty(x)), myCell);
    end
    if any(idxNaN)
        fillValue = mean(cell2mat(myCell), 'omitnan');
        for iIdx = 1:length(idxSubject)
            MotorTable.(variable){idxSubject(iIdx)} = fillValue;
        end
    end

    % Ensure that there are 4 months between the stages
    % MotorTable
    idxStage1 = (MotorTable.Subject == subject & MotorTable.Stage == 1);
    idxStage2 = (MotorTable.Subject == subject & MotorTable.Stage == 2);
    idxStage3 = (MotorTable.Subject == subject & MotorTable.Stage == 3);
    age1 = MotorTable.Age(find(idxStage1, 1));
    MotorTable.Age(idxStage2) = age1 + 4/12;
    if any(idxStage3)
        MotorTable.Age(idxStage3) = age1 + 2*4/12;
    end
    % MotorMetrics
    idxStage1 = (MotorMetrics.Subject == subject & MotorMetrics.Stage == 1);
    idxStage2 = (MotorMetrics.Subject == subject & MotorMetrics.Stage == 2);
    idxStage3 = (MotorMetrics.Subject == subject & MotorMetrics.Stage == 3);
    age1 = MotorMetrics.Age(find(idxStage1, 1));
    MotorMetrics.Age(idxStage2) = age1 + 4/12;
    if any(idxStage3)
        MotorMetrics.Age(idxStage3) = age1 + 2*4/12;
    end
end

% append table to Measurements structure
Measurements.MotorTable = MotorTable;

%% Motor table subject-wise

tasks = unique(MotorTable.Task, 'stable');
stages = unique(MotorTable.Stage, 'stable');
SkatingTable = struct([]);
iRow = 1;
newDepVars = {};
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectMotorTable = MotorTable(MotorTable.Subject == string(subject), :);
    for iStage = 1:length(stages)
        stage = stages(iStage);
        stageTable = subjectMotorTable(subjectMotorTable.Stage == stage, :);
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
    "D2_F_"
    "D2_BZO"
    "D2_KL"
    "Stroop_INT_median"
    ];
cognitionVariables = intersect(cognitionVariables_clean, cognitionVariables_orig, 'stable');
CognitionTableClean = CognitionTable_orig(:, cognitionVariables);

% rename variables
renameVariables_old = [
    "AD_MW"
    "Hyp_MW"
    "D2_F_"
    "D2_BZO"
    "D2_KL"
    "Stroop_INT_median"
    ];
renameVariables_new = [
    "AttentionDeficit"
    "Hyperactivity"
    "D2_Error"
    "D2_Completed"
    "D2_Concentration"
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
    "Age"
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

    % Age
    stage = CognitionTableClean.Stage(iRow);
    idxSubjectStage = find(MotorMetrics.Subject == subject & MotorMetrics.Stage == stage, 1, 'first');
    CognitionTable(iRow).Height = MotorMetrics.Height(idxSubjectStage);
    CognitionTable(iRow).Weight = MotorMetrics.Weight(idxSubjectStage);
    CognitionTable(iRow).Age = MotorMetrics.Age(idxSubjectStage);

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

% rename table values
variables = {'Skating', 'ADHS', 'Medication', 'Intervention'};
for iVar = 1:length(variables)
    variable = variables{iVar};
    CognitionTable.(variable) = string(CognitionTable.(variable));
    CognitionTable.(variable)(CognitionTable.(variable) == "1") = "yes";
    CognitionTable.(variable)(CognitionTable.(variable) == "0") = "no";
end
CognitionTable.Stage = string(CognitionTable.Stage);
CognitionTable.Stage = "t" + CognitionTable.Stage;

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
% rename table values
variables = {'Skating', 'ADHS', 'Medication', 'Intervention'};
for iVar = 1:length(variables)
    variable = variables{iVar};
    SkatingTable.(variable) = string(SkatingTable.(variable));
    SkatingTable.(variable)(SkatingTable.(variable) == "1") = "yes";
    SkatingTable.(variable)(SkatingTable.(variable) == "0") = "no";
end
SkatingTable.Stage = string(SkatingTable.Stage);
SkatingTable.Stage = "t" + SkatingTable.Stage;

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
        subjectMotorTable = table2struct(subjectSubjectTable);
        newVariables = setdiff(fieldnames(subjectMotorTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            SkatingTable_subjectMean(iRow).(variable) = subjectMotorTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append motor data
        subjectMotorTable = table2struct(stageMotorTable);
        newVariables = setdiff(fieldnames(subjectMotorTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            SkatingTable_subjectMean(iRow).(variable) = subjectMotorTable.(variable);
            allVariables = union(allVariables, {variable}, 'stable');
        end

        % append cognition data
        subjectMotorTable = table2struct(stageCognitionTable);
        newVariables = setdiff(fieldnames(subjectMotorTable), allVariables, 'stable');
        for iVar = 1:length(newVariables)
            variable = newVariables{iVar};
            SkatingTable_subjectMean(iRow).(variable) = subjectMotorTable.(variable);
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

%% save tables to disk

% write Subjects table
fprintf('Saving Subjects table...\n');
saveTable(Measurements.Subjects, 'SubjectsTable', {'xlsx'}, outDir);

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
