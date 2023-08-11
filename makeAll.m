init;

diary(fullfile(outDir,'log.txt'));

makeSubjects;
% makeMotorData;
% makeMotorMetrics;
makeCognitionData;
makeTables;
saveState;
% makePlots;

diary off



