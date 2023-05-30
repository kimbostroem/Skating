function makeTables

init

% load current state
loadState;

nMeas = length(Measurements.Observations); 
outDir = strrep(Measurements.outDir, '\', '/');

Observations = struct2table(Measurements.Observations);

fields = {'subject', 'stage', 'adhs', 'medication', 'height', 'weight', 'age', 'intervention', 'task', 'side', 'pathLength', 'targetError', 'meanJerk', 'meanJerkXY', 'trial'};

subjectProps = {'subject', 'adhs', 'height', 'weight', 'age'};
conditions = {'stage', 'task', 'side'};
depVars = {'pathLength', 'targetError', 'meanJerk', 'meanJerkXY'};
trialVar = 'trial';

subjects = {Measurements.Observations.subject};

for iSubject = 1:length(subjects)
    subject = subjects{iSubject};
    subjectTable = Observations(Observations.subject == string(subject), :);  
    levels = cell(1, length(conditions));    
    for iCond = 1:length(conditions)
        condition = conditions{iCond};
        levels{iCond} = unique(subjectTable.(condition));
        for iLevel = 1:length(levels{iCond})
            level = levels{iCond}(iLevel);            
            levelsStrs{iCond}(iLevel) = string(condition) + " == " + level;
        end
    end
    condCombs = string(table2cell(combinations(levelsStrs{:})));
    nCondCombs = size(condCombs, 1);
    condStrs = strings(nCondCombs, 1);    
    for iStr = 1:nCondCombs
        myTable = subjectTable;
        for iBla = 1:length(condCombs(iStr, :))
            parts = strsplit(condCombs(iStr, iBla));
            [value, isNum] = str2num(parts{3});
            if isNum
                rows = (myTable.(parts{1}) == value);
            else
                rows = (myTable.(parts{1}) == parts{3});
            end
            myTable = myTable(rows, :);
        end        
        condStrs(iStr) = strjoin(condCombs(iStr, :), " && ");
    end
end


% write original table
fprintf('\t\t- Saving Observations to table...\n');
outpath = fullfile(outDir, 'Observations.xlsx');
writetable(Observations, outpath, 'WriteMode', 'replacefile');








end