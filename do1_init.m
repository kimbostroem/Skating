% restore default path
restoredefaultpath;
% load library
addpath(genpath('library'));
% clear workspace
clear
close all

% init Measurements structure
Measurements = struct;

% set folders
Measurements.dataDir = fullfile('..', 'Skating_In');
Measurements.outDir = fullfile('..', 'Skating_Out');
Measurements.paramDir = fullfile('..', 'Skating_In');

%% Save to MAT file

fprintf('Saving Measurements to MAT file...\n');
save('Measurements', 'Measurements');
fprintf('DONE\n\n');