function makeTables()

% load current state
Measurements = loadState();

% get outDir
outDir = evalin('base', 'outDir');

%% Motor tables

MotorMetrics = Measurements.MotorMetrics;

conditions = {'Stage', 'Task'};
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
CleanTable = SourceTable;

%% Long motor table

SourceTable = CleanTable;
LongTable = struct([]);
iRow = 1;
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.Subject == string(subject), :);
    levels = cell(1, length(conditions));
    for iCond = 1:length(conditions)
        condition = conditions{iCond};
        levels{iCond} = unique(subjectTable.(condition));
    end
    condCombs = table2cell(combinations(levels{:}));
    nCondCombs = size(condCombs, 1);
    for iComb = 1:nCondCombs
        combTable = subjectTable;
        % restrict combTable iteratively to match combination of conditions
        for iCond = 1:length(condCombs(iComb, :))
            cond = conditions{iCond};
            value = condCombs{iComb, iCond};
            rows = (combTable.(cond) == value);
            combTable = combTable(rows, :);
        end
        if isempty(combTable)
            continue
        end
        initRow = table2struct(combTable(1, :));
        variables = setdiff(fieldnames(initRow), depVars);
        task = initRow.Task;
        for iVar = 1:length(variables)
            variable = variables{iVar};
            LongTable(iRow).(variable) = initRow.(variable);
        end
        for iVar = 1:length(depVars)
            depVar = depVars{iVar};
            values = combTable.(depVar);
            LongTable(iRow).(sprintf('%s_%s', task, depVar)) = mean(values, 'omitnan');
            value = std(values, 'omitnan');
            if value == 0
                value = NaN;
            end
            LongTable(iRow).(sprintf('%s_%s_std', task, depVar)) = value;
            LongTable(iRow).(sprintf('%s_%s_n', task, depVar)) = length(values);
        end        
        
        % increment row index
        iRow = iRow + 1;
    end
end
% convert structure array to table
LongTable = struct2table(LongTable);
% clean up
LongTable = removevars(LongTable, {'Trial', 'Side', 'FileName'});

%% Wide motor table

prepostVar = 'Stage';
prepostValues = [1 2];
conditions = setdiff(conditions, prepostVar);
newDepVars = depVars;
for iVar = 1:length(depVars)
    depVar = depVars{iVar};
    newDepVars = [newDepVars, {[depVar, '_std'], [depVar, '_n']}]; %#ok<AGROW>
end
depVars = newDepVars;

SourceTable = LongTable;
WideTable = table;
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.Subject == string(subject), :);
    levels = cell(1, length(conditions));
    for iCond = 1:length(conditions)
        condition = conditions{iCond};
        levels{iCond} = unique(subjectTable.(condition));
    end
    condCombs = table2cell(combinations(levels{:}));
    nCondCombs = size(condCombs, 1);
    for iComb = 1:nCondCombs
        combTable = subjectTable;
        for iCond = 1:length(condCombs(iComb, :))
            cond = conditions{iCond};
            value = condCombs{iComb, iCond};
            rows = (combTable.(cond) == value);
            combTable = combTable(rows, :);
        end
        if isempty(combTable)
            continue
        end
        myTable = combTable(combTable.(prepostVar) == prepostValues(1), :);
        if isempty(myTable)
            fprintf('\t- No pre measurement for condition (%s, %s) -> skipping\n', subject, strjoin(string(condCombs(iComb, :)), ', '));
            continue
        end
        myTable = combTable(combTable.(prepostVar) == prepostValues(2), :);
        if isempty(myTable)
            fprintf('\t- No post measurement for condition (%s, %s) -> skipping\n', subject, strjoin(string(condCombs(iComb, :)), ', '));
            continue
        end
        tableRow = myTable;
        for iVar = 1:length(depVars)
            depVar = depVars{iVar};
            value = combTable.(depVar)(combTable.(prepostVar) == prepostValues(1));
            if isempty(value)
                value = NaN;
            end
            tableRow.([depVar, '_pre']) = value;
            value = combTable.(depVar)(combTable.(prepostVar) == prepostValues(2));
            if isempty(value)
                value = NaN;
            end
            tableRow.([depVar, '_post']) = value;
            tableRow = removevars(tableRow, depVar);
        end
        WideTable = [WideTable; tableRow]; %#ok<AGROW>
    end
end
WideTable = removevars(WideTable, {prepostVar});

%% Reduce motor tables

rmVars = {'subjectCode', 'date', 'Beidbein_start', 'Beidbein_stop', 'Einbein_start', 'Einbein_stop', 'doneData', 'doneMetrics', 'donePlots'};
% clean table
myRmVars = intersect(CleanTable.Properties.VariableNames, rmVars);
CleanTable = removevars(CleanTable, myRmVars);
% long mean table
myRmVars = intersect(LongTable.Properties.VariableNames, rmVars);
LongTable = removevars(LongTable, myRmVars);
% wide mean table
myRmVars = intersect(WideTable.Properties.VariableNames, rmVars);
WideTable = removevars(WideTable, myRmVars);


%% Cognitive tables

% % delete irrelevant columns
% variables_orig = CognitionData.Properties.VariableNames;
% variables_clean = [
%     "Probandencode"
%     "ADHS"
%     "Stage"
%     "Intervention"
%     "Sex"
%     "Age_yrs"
%     "Height_cm"
%     "Weight_kg"
%     "Medikation"
%     "AD_MW"
%     "Hyp_MW"
%     "D2_F__SW"
%     "D2_BZO_SW"
%     "D2_KL_SW"
%     "Stroop_FWL_SW"
%     "Stroop_FSB_SW"
%     "Stroop_INT_SW"
%     ];
% idx = ~ismember(variables_orig, variables_clean);
% CognitionData(:, idx) = [];

%% Append tables to Measurements structure

Measurements.CleanTable = CleanTable;
Measurements.LongTable = LongTable;
Measurements.WideTable = WideTable;

%% save Measurements structure

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

%% save tables to disk

% write original table
fprintf('Saving observations to table...\n');
saveTable(MotorMetrics, 'Observations_orig', {'xlsx'}, outDir);

% write clean table
fprintf('Saving clean observations to table...\n');
saveTable(CleanTable, 'Observations_clean', {'csv'}, outDir);

% write mean table
fprintf('Saving averaged observations in long format to table...\n');
saveTable(LongTable, 'Observations_long', {'csv'}, outDir);

% write wide table
fprintf('Saving averaged observations in wide format to table...\n');
saveTable(WideTable, 'Observations_wide', {'csv'}, outDir);

end
