% load library
addpath(genpath('library'));

% set folders
dataDir = '/Users/kbostroem/sciebo/Skating/Skating_In';
outDir = '/Users/kbostroem/sciebo/Skating/Skating_Out';
paramDir = '/Users/kbostroem/sciebo/Skating/Auswertung';
isLoadFiles = 1;

[MeasurementInfo, MeasurementData] = prepareData(dataDir, paramDir, outDir, isLoadFiles);