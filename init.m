%% Initilization

fprintf('Initializing...\n');

% clear workspace
clear
% close all open figures
close all

% which dataset to use (subfolder of outDir)
dataSet = '.'; % PR, PB, PRCO, ProKo, .

fprintf('Using dataset ''%s''...\n', dataSet);

% set folders
allDir = '.';
parentOutDir = '../Data_Out';
outDir = fullfile(parentOutDir, dataSet);
inDir = '../Data_In';
motorDir = fullfile(inDir, 'Motorisch', dataSet);
cogDir = fullfile(inDir, 'Kognitiv');
paramDir = fullfile(inDir, 'Parameter');

% create output folder if necessary
if ~isfolder(outDir)
    mkdir(outDir);
end

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));

% suppress warning when table headers do not conform Matlab variable name standard
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');

% load state
% Measurements = loadState();