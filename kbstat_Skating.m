%% Initilization

%       distribution    Distribution used for the GLM fit.
%                       OPTIONAL, default = 'Normal'.
%                       Possible values:
%                       'normal'	        Normal distribution
%                       'logNormal'         Normal distribution on log-values
%                       'binomial'	        Binomial distribution
%                       'poisson'	        Poisson distribution
%                       'gamma'	            Gamma distribution
%                       'inverseGaussian'	Inverse Gaussian distribution
%
%       link            Link function used for the GLM fit.
%                       OPTIONAL, default depends on chosen distribution.
%                       Possible values:
%                       'identity'	    g(mu) = mu.             Default for Normal distribution.
%                       'log'	        g(mu) = log(mu).        Default for Poisson, Gamma, and InverseGaussian.
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution.
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)
%                       Scalar p	    g(mu) = mu.^p           Canonical for Gamma (p = -1) and InverseGaussian (p= -2)

%% Init

fprintf('Initializing...\n');

% clear workspace
clear
close all

% restore default path
restoredefaultpath;

% add library and subfolders to path
addpath(genpath('../kbstat'));
resultsDir = '../Skating_Stats';

%% Global options

options = struct;
options.x = 'Intervention, ADHS';
options.id = 'Subject';
options.within = 'Intervention';
options.interact = 'Intervention, ADHS';
options.posthocMethod = 'ttest';
options.removeOutliers = 'true';
options.isRescale = true;
options.errorBars = 'se';
% options.constraint = 'Medication == 0';
% options.constraint = 'ADHS == 1 & Stage == 2';

%% Analysis of Motor data

depVars = {'Jerk', 'TargetError', 'PathLength'};
depVarUnitss = {'m/s^3', 'm', 'm'};
tasks = {'Balance', 'Sprung', 'Einbein'};
distributions = {'gamma', 'gamma', 'gamma'};
links = {'', '', ''};
meanFlags = [1, 0];

optionsOrig = options;
if isfield(options, 'constraint')
    constraintOrig = options.constraint;
else
    constraintOrig = '';
end
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

%% Analysis of Cognition data

depVars = {'D2_Error', 'Stroop', 'AttentionDeficit', 'Hyperactivity'};
depVarUnitss = {'', '', '', ''};
distributions = {'gamma', 'gamma', 'normal', 'normal'};
links = {'', '', '', ''};

optionsOrig = options;
if isfield(options, 'constraint')
    constraintOrig = options.constraint;
else
    constraintOrig = '';
end
for iVar = 1:length(depVars)
    depVar = depVars{iVar};
    depVarUnits = depVarUnitss{iVar};
    distribution = distributions{iVar};
    link = links{iVar};
    for iTask = 1:length(tasks)
        task = tasks{iTask};
        for iFlag = 1:length(meanFlags)
            options.inFile = '../Skating_Out/All/CognitionTable.csv';
            options.y = depVar;
            options.outDir = sprintf('%s/Cognition/%s', resultsDir, depVar);
            options.yUnits = depVarUnits;
            options.distribution = distribution;
            options.link = link;
            kbstat(options);
            options = optionsOrig;
        end
    end
end
