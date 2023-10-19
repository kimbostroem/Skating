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
%                       'identity'	    g(mu) = mu.             Default for Normal distribution
%                       'log'	        g(mu) = log(mu).        Default for Poisson
%                       'logit'	        g(mu) = log(mu/(1-mu))  Default for Binomial distribution
%                       'loglog'	    g(mu) = log(-log(mu))
%                       'probit'	    g(mu) = norminv(mu)
%                       'comploglog'	g(mu) = log(-log(1-mu))
%                       'reciprocal'	g(mu) = mu.^(-1)        Default for Gamma
%                       Scalar p	    g(mu) = mu.^p           Default for InverseGaussian (p= -2)

%% Init

fprintf('Initializing...\n');

% clear workspace
clear
close all

% restore default path
restoredefaultpath;

% add library and subfolders to path
addpath(genpath('../kbstat'));
resultsDir = '../Statistics';

%% Global options

options = struct;
options.id = 'Subject';
options.showVarNames = 1;
options.levelOrder = 'sorted';
options.markerSize = 4;

%%-------------------------------------------------------------------------
%% Hypothesis 1: ADHS effect

hypoDir = 'H1 ADHS effect for Skating';

options.constraint = 'Skating == yes & Stage ~= t3';
options.x = 'ADHS, Stage';
options.within = '';
options.interact = '';
options.fitMethod = 'REMPL';
options.transform = 'q50';

%% Motor data

options.inFile = '../Data_Out/SkatingTable.csv';
options.multiVar = 'MotorTask';
options.posthocMethod = 'emm';

options.y = {
    'TargetError'
    };
options.yUnits = {
    'm'
    };
options.distribution = {
    'gamma'
    };
options.outDir = sprintf('%s/%s/Motoric', resultsDir, hypoDir);
kbstat(options);

%% Cognition

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'D2_Completed'
    'D2_Concentration'
    'Stroop'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Cognition', resultsDir, hypoDir);
kbstat(options);

%% Symptomatics

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'AttentionDeficit'
    'Hyperactivity'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Symptomatics', resultsDir, hypoDir);
kbstat(options);

%%-------------------------------------------------------------------------
%% Hypothesis 2: Time effect

hypoDir = 'H2 Time effect for ADHS';

options.constraint = 'Stage ~= t3 & ADHS == yes';
options.x = 'Stage, Skating';
options.within = 'Stage';
options.interact = '';
options.fitMethod = 'REMPL';

%% Motor data

options.inFile = '../Data_Out/SkatingTable.csv';
options.multiVar = 'MotorTask';
options.posthocMethod = 'emm';

options.y = {
    'TargetError'
    };
options.yUnits = {
    'm'
    };
options.distribution = {
    'gamma'
    };
options.outDir = sprintf('%s/%s/Motoric', resultsDir, hypoDir);
kbstat(options);

%% Cognition

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'D2_Completed'
    'D2_Concentration'
    'Stroop'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Cognition', resultsDir, hypoDir);
kbstat(options);

%% Symptomatics

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'AttentionDeficit'
    'Hyperactivity'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Symptomatics', resultsDir, hypoDir);
kbstat(options);

%%-------------------------------------------------------------------------
%% Hypothesis 3: Skating effect

hypoDir = 'H3 Skating effect for ADHS at t2';

options.constraint = 'Stage == t2 & ADHS == yes';
options.x = 'Skating, Medication';
options.within = '';
options.interact = '';
options.fitMethod = 'REMPL';

%% Motor data

options.inFile = '../Data_Out/SkatingTable.csv';
options.multiVar = 'MotorTask';
options.posthocMethod = 'emm';

options.y = {
    'TargetError'
    };
options.yUnits = {
    'm'
    };
options.distribution = {
    'gamma'
    };
options.outDir = sprintf('%s/%s/Motoric', resultsDir, hypoDir);
kbstat(options);

%% Cognition

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'D2_Completed'
    'D2_Concentration'
    'Stroop'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Cognition', resultsDir, hypoDir);
kbstat(options);

%% Symptomatics

options.inFile = '../Data_Out/CognitionTable.csv';
options.posthocMethod = 'emm';

options.y = {
    'AttentionDeficit'
    'Hyperactivity'
    };
options.yUnits = {
    'Score'
    };
options.distribution = {
    'normal'
    };
options.outDir = sprintf('%s/%s/Symptomatics', resultsDir, hypoDir);
kbstat(options);
