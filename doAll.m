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

% get info
Measurements = makeInfo(Measurements);

% get data
Measurements = makeData(Measurements);

% get metrics
Measurements = makeMetrics(Measurements);

% save everything
outpath = fullfile(Measurements.outDir, 'Measurements');
save(outpath, 'Measurements');
