function options = kbstatWithOptions(depVars, meanFlags, tasks, distributions, links, depVarUnitss, resultsDir, options)

if nargin < 8
    options = struct;
end

optionsOrig = options;

if nargin < 7
    resultsDir = '';
end
if nargin < 6
    depVarUnitss = '';
end
if nargin < 5
    links = '';
end
if nargin < 4
    distributions = 'normal';
end
if nargin < 3
    tasks = '';
end
if nargin < 2
    meanFlags = 1;
end

if isfield(options, 'constraint')
    constraintOrig = options.constraint;
else
    constraintOrig = '';
end

depVars = cellstr(depVars);
tasks = cellstr(tasks);
distributions = cellstr(distributions);
links = cellstr(links);
depVarUnitss = cellstr(depVarUnitss);

for iVar = 1:length(depVars)
    depVar = depVars{iVar};
    depVarUnits = depVarUnitss{iVar};
    distribution = distributions{iVar};
    link = links{iVar};
    for iTask = 1:length(tasks)
        task = tasks{iTask};
        for iFlag = 1:length(meanFlags)
            meanFlag = meanFlags(iFlag);
            if meanFlag
                options.inFile = '../Skating_Out/All/SkatingTable_subjectMean.csv';
                if ~isempty(task)
                    options.y = sprintf('%s_%s', task, depVar);
                    options.outDir = sprintf('%s/SubjectMean/%s_%s', resultsDir, task, depVar);
                else
                    options.y = depVar;
                    options.outDir = sprintf('%s/SubjectMean/%s', resultsDir, depVar);
                end
                
            else
                options.inFile = '../Skating_Out/All/SkatingTable.csv';
                if ~isempty(task)
                    options.y = depVar;
                    if ~isempty(constraintOrig)
                        options.constraint = sprintf('%s & MotorTask == %s', constraintOrig, task);
                    else
                        options.constraint = sprintf('MotorTask == %s', task);
                    end
                    options.title = sprintf('%s %s', task, depVar);
                    options.outDir = sprintf('%s/NoSubjectMean/%s_%s', resultsDir, task, depVar);
                else
                    options.y = depVar;
                    options.outDir = sprintf('%s/NoSubjectMean/%s', resultsDir, depVar);
                end
                
            end
            options.yUnits = depVarUnits;
            options.distribution = distribution;
            options.link = link;
            kbstat(options);
            options = optionsOrig;
        end
    end
end



end