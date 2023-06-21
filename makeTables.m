function makeTables()

init

% load current state
Measurements = loadState();

outDir = strrep(Measurements.outDir, '\', '/');

OrigTable = struct2table(Measurements.Observations);

conditions = {'stage', 'task'};
depVars = {'pathLength', 'targetError', 'meanJerk', 'meanJerkXY'};
subjects = unique([Measurements.Observations.subject]);

%% create clean table

SourceTable = OrigTable;
% remove empty "isValid" rows
rows = ~cellfun(@isempty,SourceTable.isValid);
SourceTable = SourceTable(rows, :);
% remove rows where "isValid" equals zero
rows =  (cell2mat(SourceTable.isValid) == 1);
SourceTable = SourceTable(rows, :);
SourceTable = removevars(SourceTable, {'isValid'});
CleanTable = SourceTable;

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

%% create long mean table

SourceTable = CleanTable;
LongTable = table;
for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = SourceTable(SourceTable.subject == string(subject), :);
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
        tableRow = combTable(1, :);
        for iVar = 1:length(depVars)
            depVar = depVars{iVar};
            values = cell2mat(combTable.(depVar));
            tableRow.(depVar) = mean(values, 'omitnan');
            value = std(values, 'omitnan');
            if value == 0
                value = NaN;
            end
            tableRow.([depVar, '_std']) = value;
            tableRow.([depVar, '_n']) = length(values);
        end
        tableRow = removevars(tableRow, {'trial', 'side', 'fileName'});
        LongTable = [LongTable; tableRow]; %#ok<AGROW>
    end
end

%% create wide mean table

prepostVar = 'stage';
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
    subjectTable = SourceTable(SourceTable.subject == string(subject), :);
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

%% reduce tables

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



%% save tables

% write original table
fprintf('Saving observations to table...\n');
saveTable(OrigTable, 'Observations_orig', {'xlsx'}, outDir);

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
