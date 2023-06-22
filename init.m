% clear workspace
clear
close all

inDir = '../Skating_In';
motorDir = fullfile(inDir, 'Motorisch');
cogDir = fullfile(inDir, 'Kognitiv');
paramDir = fullfile(inDir, 'Parameter');
outDir = '../Skating_Out';

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));

% suppress warning when table headers do not conform Matlab variable name standard
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

Measurements = loadState();