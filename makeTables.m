function makeTables()

fprintf('\nMaking Tables...\n');


% load current state
Measurements = loadState();

% get outDir
outDir = evalin('base', 'outDir');

%% Motor tables

MotorMetrics = Measurements.MotorMetrics;

depVars = {'PathLength', 'TargetError', 'Jerk', 'JerkXY', 'JerkZ'};
subjects = unique(Measurements.Subjects.Subject, 'stable');

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

%% Motor table

tasks = unique(MotorTable_all.Task, 'stable');
stages = unique(MotorTable_all.Stage, 'stable');
SourceTable = MotorTable_all;
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
MotorTable_long = TargetTable;

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
CognitionTable_all = CognitionTable_orig(:, cognitionVariables);

%% Subject table

subjectVariables_orig = Measurements.Subjects.Properties.VariableNames;
subjectVariables_clean = [
    "Subject"
    "ADHS"
    "ADS"
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
SubjectTable_all = Measurements.Subjects(:, subjectVariables);

%% Combined table long format (pre and post values in rows)

SourceMotorTable = MotorTable_long;
SourceCognitionTable = CognitionTable_all;

TargetTable = struct([]);
iRow = 1;
subjects = unique(SourceMotorTable.Subject, 'stable');
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectMotorTable = SourceMotorTable(SourceMotorTable.Subject == string(subject), :);
    subjectCognitionTable = SourceCognitionTable(SourceCognitionTable.Subject == string(subject), :);
    subjectSubjectTable = SubjectTable_all(SubjectTable_all.Subject == string(subject), :);
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
TargetTable = struct2table(TargetTable);
SkatingTable_long = TargetTable;

%% Fill in missing values that can be derived

% replace missing Diagnose with ADHS
nRows = size(SkatingTable_long, 1);
for iRow = 1:nRows
    variable = 'Diagnose';
    fillValue = 2;
    value = SkatingTable_long{iRow, variable};
    if isnan(value)
        SkatingTable_long.(variable)(iRow) = fillValue;
    end
end

% replace missing height at some stage with mean of heights of all stages
SourceTable = SkatingTable_long;
TargetTable = SkatingTable_long;
subjects = unique(Measurements.Subjects.Subject, 'stable');
nSubjects = length(subjects);
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    idx = find(SourceTable.Subject == string(subject));
    variable = 'Height';
    myCell = SourceTable{idx, variable};
    idxNaN = cellfun(@(x) any(isempty(x)), myCell);
    if any(idxNaN)
        fillValue = mean(cell2mat(myCell), 'omitnan');
        for iIdx = 1:length(idx)
            TargetTable.(variable){idx(iIdx)} = fillValue;
        end
    end
end
SkatingTable_long = TargetTable;

%% Rename variables

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
SkatingTable_long = renamevars(SkatingTable_long, renameVariables_old, renameVariables_new);
for iVar = 1:length(renameVariables_old)
    depVars = cellstr(strrep(depVars, renameVariables_old(iVar), renameVariables_new(iVar)));
end

%% Combined table wide format (pre and post values in columns)

SourceTable = SkatingTable_long;
TargetTable = struct([]);
nSubjects = length(subjects);
iRow = 1;
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    subjectSourceTable = SourceTable(SourceTable.Subject == string(subject), :);
    if isempty(subjectSourceTable)
        continue
    end
    tic
    initRow = table2struct(subjectSourceTable(1, :));
    variables = setdiff(fieldnames(initRow), [depVars, {'Stage'}], 'stable');
    for iVar = 1:length(variables)
        variable = variables{iVar};
        TargetTable(iRow).(variable) = initRow.(variable);
    end
    for iVar = 1:length(depVars)
        depVar = depVars{iVar};
        % pre
        value = subjectSourceTable.(depVar)(subjectSourceTable.Stage == stages(1));
        if isempty(value)
            value = NaN;
        elseif iscell(value)
            value = cell2mat(value);
        end
        TargetTable(iRow).([depVar, '_pre']) = value;
        % post
        value = subjectSourceTable.(depVar)(subjectSourceTable.Stage == stages(2));
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

Measurements.SkatingTable_long = SkatingTable_long;
Measurements.SkatingTable_wide = SkatingTable_wide;

%% Create Skating PB tables (only subjects from the PB group)

fprintf('\t\t- Creating PB tables...\n');
SkatingTable_PB_long = SkatingTable_long(startsWith(SkatingTable_long.Subject, 'PB'), :);
SkatingTable_PB_wide = SkatingTable_long(startsWith(SkatingTable_wide.Subject, 'PB'), :);


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
saveTable(Measurements.MotorMetrics, 'MotorTable', {'xlsx'}, outDir);

% write Cognition table
fprintf('Saving Cognition table...\n');
saveTable(Measurements.CognitionData, 'CognitionTable', {'xlsx'}, outDir);

% write Skating table in long format
fprintf('Saving Skating table in long format...\n');
saveTable(SkatingTable_long, 'SkatingTable_long', {'csv'}, outDir);

% write Skating table in wide format
fprintf('Saving Skating table in wide format...\n');
saveTable(SkatingTable_wide, 'SkatingTable_wide', {'csv'}, outDir);

% write Skating PB table in long format
fprintf('Saving Skating PB table in long format...\n');
saveTable(SkatingTable_PB_long, 'SkatingTable_PB_long', {'csv'}, outDir);

% write Skating PB table in wide format
fprintf('Saving Skating PB table in wide format...\n');
saveTable(SkatingTable_PB_wide, 'SkatingTable_PB_wide', {'csv'}, outDir);

fprintf('If necessary, save current state using ''saveState''\n');


end
