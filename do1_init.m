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
Measurements.dataDir = '/Users/kbostroem/sciebo/Skating/Skating_In';
Measurements.outDir = '/Users/kbostroem/sciebo/Skating/Skating_Out';
Measurements.paramDir = '/Users/kbostroem/sciebo/Skating/Auswertung';

%% Save to MAT file

fprintf('Saving Measurements to MAT file...\n');
save('Measurements', 'Measurements');
fprintf('DONE\n\n');