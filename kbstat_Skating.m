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

%% Analysis

depVars = {'Jerk', 'TargetError'};
tasks = {'Balance', 'Sprung'};
depVarUnitss = {'m/s^3', 'm'};
distributions = {'gamma', 'gamma'};
links = {'', ''};
meanFlags = [1, 0];
kbstatWithOptions(depVars, meanFlags, tasks, distributions, links, depVarUnitss, resultsDir, options);

depVars = {'D2_Error', 'Stroop', 'AttentionDeficit', 'Hyperactivity'};
tasks = '';
depVarUnitss = {'', '', '', ''};
distributions = {'gamma', 'gamma', 'normal', 'normal'};
links = {'', '', '', ''};
meanFlags = [1, 0];
kbstatWithOptions(depVars, meanFlags, tasks, distributions, links, depVarUnitss, resultsDir, options);
