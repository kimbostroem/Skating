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
options.inFile = '../Skating_Out/All/SkatingTable.csv';
options.id = 'Subject';
options.posthocMethod = 'emm';
options.removeOutliers = 'true';
options.isRescale = true;

metrics = {
    'TargetError'
    };
tasks = {
    'Balance'
    'Sprung'
    'Einbein'
    };
depVarUnitss = {
    'm'
    'm'
    'm'
    };
distributions = {
    'gamma'
    'gamma'
    'gamma'
    };

controlGroups = {'ADHS', 'Skating'};

for iControl = 1:length(controlGroups)
    controlGroup = controlGroups{iControl};

    switch controlGroup
        case 'ADHS'            
            options.x = 'Stage, Skating, Medication';
            options.within = 'Stage';
            options.interact = 'Stage, Skating, Medication';
            options.constraint = 'ADHS == yes & Stage ~= t3';
        case 'Skating'
            options.x = 'Stage, ADHS';
            options.within = 'Stage';
            options.interact = 'Stage, ADHS';
            options.constraint = 'Skating == yes & Stage ~= t3';
    end
    

    %% Analysis of Motor data

    optionsOrig = options;
    constraintOrig = options.constraint;
    iDepVar = 0;
    for iMetric = 1:length(metrics)
        metric = metrics{iMetric};
        for iTask = 1:length(tasks)
            task = tasks{iTask};
            iDepVar = iDepVar+1;

            switch controlGroup
                case 'ADHS'
                    depVar = sprintf('%s_%s_ADHS', task, metric);
                case 'Skating'
                    depVar = sprintf('%s_%s_Skating', task, metric);
            end

            options.y = metric;
            options.constraint = sprintf('%s & MotorTask == %s', constraintOrig, task);
            options.title = sprintf('%s', strrep(depVar, '_', ' '));
            options.outDir = sprintf('%s/Motoric/%s', resultsDir, depVar);

            if length(depVarUnitss) == 1
                options.yUnits = depVarUnitss{1};
            else
                options.yUnits = depVarUnitss{iDepVar};
            end
            if exist('distributions', 'var') && length(distributions) == 1
                options.distribution = distributions{1};
            elseif exist('distributions', 'var')
                options.distribution = distributions{iDepVar};
            end
            if exist('links', 'var') && length(links) == 1
                options.link = links{1};
            elseif exist('links', 'var')
                options.link = links{iDepVar};
            end
            kbstat(options);
            options = optionsOrig;
        end
    end

end

return

%% Analysis of Cognition data

metrics = {
    'D2_Completed'
    'D2_Concentration'
    'Stroop'
    'AttentionDeficit'
    'Hyperactivity'
    };
depVarUnitss = {
    ''
    ''
    ''
    ''
    ''};
distributions = {
    'gamma'
    'gamma'
    'gamma'
    'normal'
    'normal'
    };

optionsOrig = options;
if isfield(options, 'constraint')
    constraintOrig = options.constraint;
else
    constraintOrig = '';
end
for iMetric = 1:length(metrics)
    metric = metrics{iMetric};
    options.y = metric;
    if exist('distributions', 'var') && length(distributions) == 1
        options.distribution = distributions{1};
    elseif exist('distributions', 'var')
        options.distribution = distributions{iMetric};
    end
    if exist('links', 'var') && length(links) == 1
        options.link = links{1};
    elseif exist('links', 'var')
        options.link = links{iMetric};
    end
    if exist('depVarUnitss', 'var') && length(depVarUnitss) == 1
        options.yUnits = depVarUnitss{1};
    elseif exist('depVarUnitss', 'var')
        options.yUnits = lidepVarUnitssnks{iMetric};
    end
    for iTask = 1:length(tasks)
        task = tasks{iTask};
        options.inFile = '../Skating_Out/All/CognitionTable.csv';
        options.outDir = sprintf('%s/Cognition/%s', resultsDir, metric);
        kbstat(options);
        options = optionsOrig;
    end
end
