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
options.fitMethod = 'REMPL';
options.posthocMethod = 'emm';
% options.transform = 'q50';
options.showVarNames = 'names_and_levels';
options.levelOrder = 'sorted';
options.markerSize = 4;

%%-------------------------------------------------------------------------
%% Hypothesis 1: ADHS effect

hypoDir = 'H1 ADHS effect for Skating';

% options.constraint = 'Skating == yes & Stage ~= t3';
% options.x = 'ADHS, Stage';
% options.within = 'Stage';
% options.coVar = 'Age';

options.constraint = 'Skating == yes & Stage == t1';
options.x = 'ADHS';
options.coVar = 'Age';


%% Motor data

options.inFile = '../Data_Out/SkatingTable.csv';
options.multiVar = 'MotorTask';

options.outlierMethod = 'quartiles';
options.preOutlierMethod = 'quartiles';

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

options.outlierMethod = 'none';
options.preOutlierMethod = 'none';

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

options.outlierMethod = 'none';
options.preOutlierMethod = 'none';

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
options.coVar = 'Age, Medication';

%% Motor data

options.inFile = '../Data_Out/SkatingTable.csv';
options.multiVar = 'MotorTask';

options.outlierMethod = 'quartiles';
options.preOutlierMethod = 'quartiles';

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

options.outlierMethod = 'none';
options.preOutlierMethod = 'none';

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

options.outlierMethod = 'none';
options.preOutlierMethod = 'none';

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