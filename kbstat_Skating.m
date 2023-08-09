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

fprintf('Initializing...\n');

% clear workspace
clear
close all

% restore default path
restoredefaultpath;

% add library and subfolders to path
addpath(genpath('../kbstat'));
resultsDir = '../Results';


options = struct;
options.inFile = 'SkatingTable_long.csv';
options.x = 'Intervention, ADHS';
options.id = 'Subject';
options.within = 'Intervention';
options.interact = 'Intervention, ADHS';
options.posthocMethod = 'ttest';
options.removeOutliers = 'true';
options.isRescale = true;
options.errorBars = 'se';
% options.constraint = 'Stage < 3';
% options.constraint = 'Medication == 0';

depVar = 'Balance_Jerk';
options.yUnits = 'm/s^3';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'gamma';
kbstat(options);

depVar = 'Balance_TargetError';
options.yUnits = 'm';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'gamma';
kbstat(options);

depVar = 'Sprung_TargetError';
options.yUnits = 'm';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'gamma';
kbstat(options);

depVar = 'D2_Error';
options.yUnits = '';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'gamma';
kbstat(options);

depVar = 'Stroop';
options.yUnits = '';
options.outDir = sprintf('../Results/%s', depVar);
options.y = depVar;
options.distribution = 'gamma';
kbstat(options);

depVar = 'AttentionDeficit';
options.yUnits = '';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'normal';
kbstat(options);

depVar = 'Hyperactivity';
options.yUnits = '';
options.outDir = sprintf('%s/%s', resultsDir, depVar);
options.y = depVar;
options.distribution = 'normal';
kbstat(options);
