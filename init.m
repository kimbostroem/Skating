%% Initilization

fprintf('Initializing...\n');

% clear workspace
clear
% close all open figures
close all

% which dataset to use (subfolder of outDir)
dataSet = 'All';

fprintf('Using dataset ''%s''...\n', dataSet);

% set folders
outDir = sprintf('../Skating_Out/%s', dataSet);
inDir = '../Skating_In';
motorDir = fullfile(inDir, 'Motorisch');
cogDir = fullfile(inDir, 'Kognitiv');
paramDir = fullfile(inDir, 'Parameter');

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));

% suppress warning when table headers do not conform Matlab variable name standard
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

% load state
Measurements = loadState();