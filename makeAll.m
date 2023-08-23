dataSets = {
    'PR'
    'PB'
    'PRCO'
    'ProKo'
    'All'
    };

for iData = 1:length(dataSets)
    clear Measurements;
    close all;
    dataSet = dataSets{iData};

    % if dataSet is 'All' -> join dataSets and quit
    if strcmp(dataSet, allDir)
        diary(fullfile(outDir,'log.txt'));
        fprintf('Joining dataSets into ''%s''...\n', dataSet);
        joinData();
        diary off
        break
    end

    % process dataSet
    
    fprintf('Using dataset ''%s''...\n', dataSet);
    % set folders
    outDir = fullfile(parentOutDir, dataSet);
    motorDir = fullfile(inDir, 'Motorisch', dataSet);
    % create output folder if necessary
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    diary(fullfile(outDir,'log.txt'));

    Measurements = loadState();
    makeSubjects;
    makeMotorData;
    makeMotorMetrics;
    makeCognitionData;
    makeTables;
    saveState;
    % makePlots;

    diary off

end



