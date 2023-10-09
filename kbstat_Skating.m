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
options.plotStyle = 'violin';
options.isRescale = 1;
options.separateMulti = 1;
options.transform = 'q50';
options.fitMethod = 'none';
options.posthocMethod = 'utest';
options.outlierMethod = 'auto';
options.removeOutliers = 'prepost';
controlGroups = {'ADHS'};
checkMedication = 0;

%% Analysis of Motor data

options.inFile = '../Skating_Out/All/SkatingTable.csv';
options.id = 'Subject';
options.multiVar = 'MotorTask';

metrics = {
    'TargetError'
    };
depVarUnitss = {
    'm'
    };
distributions = {
    'gamma'
    };

for iControl = 1:length(controlGroups)
    controlGroup = controlGroups{iControl};
    switch controlGroup
        case 'ADHS'
            if checkMedication
                options.x = 'Stage, Skating, Medication'; %#ok<*UNRCH>
                options.interact = 'Stage, Skating, Medication';
            else
                options.x = 'Stage, Skating';
                options.interact = 'Stage, Skating';
            end
            options.within = 'Stage';
            options.constraint = 'ADHS == yes & Stage ~= t3';
        case 'Skating'
            options.x = 'Stage, ADHS';
            options.within = 'Stage';
            options.interact = 'Stage, ADHS';
            options.constraint = 'Skating == yes & Stage ~= t3';
    end

    optionsOrig = options;
    constraintOrig = options.constraint;
    iDepVar = 0;
    for iMetric = 1:length(metrics)
        metric = metrics{iMetric};

        iDepVar = iDepVar+1;
        options.y = metric;

        switch controlGroup
            case 'ADHS'
                depVar = sprintf('%s_ADHS', metric);
                if length(metrics) > 1
                    options.outDir = sprintf('%s/Motoric/%s/ADHS', resultsDir, metric);
                else
                    options.outDir = sprintf('%s/Motoric/ADHS', resultsDir);
                end
            case 'Skating'
                depVar = sprintf('%s_Skating', metric);
                if length(metrics) > 1
                    options.outDir = sprintf('%s/Motoric/%s/Skating', resultsDir, metric);
                else
                    options.outDir = sprintf('%s/Motoric/Skating', resultsDir);
                end
        end

        % options.title = sprintf('%s', strrep(depVar, '_', ' '));

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

%% Analysis of Cognition data

options.inFile = '../Skating_Out/All/CognitionTable.csv';

options.y = {
    'D2_Completed'
    'D2_Concentration'
    'Stroop'
    };
options.yUnits = 'Score';
options.distribution = 'gamma';

for iControl = 1:length(controlGroups)
    controlGroup = controlGroups{iControl};
    switch controlGroup
        case 'ADHS'
            options.outDir = sprintf('%s/Cognition/ADHS', resultsDir);
            if checkMedication
                options.x = 'Stage, Skating, Medication';
                options.interact = 'Stage, Skating, Medication';
            else
                options.x = 'Stage, Skating';
                options.interact = 'Stage, Skating';
            end
            options.within = 'Stage';
            options.constraint = 'ADHS == yes & Stage ~= t3';
        case 'Skating'
            options.outDir = sprintf('%s/Cognition/Skating', resultsDir);
            options.x = 'Stage, ADHS';
            options.within = 'Stage';
            options.interact = 'Stage, ADHS';
            options.constraint = 'Skating == yes & Stage ~= t3';
    end

    optionsOrig = options;
    if isfield(options, 'constraint')
        constraintOrig = options.constraint;
    else
        constraintOrig = '';
    end
    kbstat(options);
    options = optionsOrig;
end

%% Analysis of Symptomatics

options.inFile = '../Skating_Out/All/CognitionTable.csv';

options.y = {
    'AttentionDeficit'
    'Hyperactivity'
    };
options.yUnits = 'Score';
options.distribution = 'normal';

for iControl = 1:length(controlGroups)
    controlGroup = controlGroups{iControl};
    switch controlGroup
        case 'ADHS'
            options.outDir = sprintf('%s/Symptomatics/ADHS', resultsDir);
            if checkMedication
                options.x = 'Stage, Skating, Medication';
                options.interact = 'Stage, Skating, Medication';
            else
                options.x = 'Stage, Skating';
                options.interact = 'Stage, Skating';
            end
            options.within = 'Stage';
            options.constraint = 'ADHS == yes & Stage ~= t3';
        case 'Skating'
            options.outDir = sprintf('%s/Symptomatics/Skating', resultsDir);
            options.x = 'Stage, ADHS';
            options.within = 'Stage';
            options.interact = 'Stage, ADHS';
            options.constraint = 'Skating == yes & Stage ~= t3';
    end

    optionsOrig = options;
    if isfield(options, 'constraint')
        constraintOrig = options.constraint;
    else
        constraintOrig = '';
    end    
    kbstat(options);
    options = optionsOrig;
end


