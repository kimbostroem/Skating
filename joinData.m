function joinData()

fprintf('\nJoining Data...\n');

allDir = 'All';

% get outDir
outDir = evalin('base', 'outDir');

if ~isfolder(fullfile(outDir, allDir))
    mkdir(fullfile(outDir, allDir));
end

dirInfo = dir(outDir);
outDirs = {dirInfo.name}';
idxExclude = startsWith(outDirs', {'.', '~', allDir}) | ~[dirInfo.isdir];
outDirs(idxExclude) = [];

tableNames = {'MotorTable', 'SkatingTable', 'SkatingTable_subjectMean', 'MotorMetrics'};
Measurements = struct;
for iTable = 1:length(tableNames)
    tableName = tableNames{iTable};
    Measurements.(tableName) = struct([]);
end
Measurements.MotorData = struct([]);

for iDir = 1:length(outDirs)
    myOutDir = fullfile(outDir, outDirs{iDir});

    % load Measurements structure
    fprintf('Loading Measurements structure from %s...\n', myOutDir);
    tic    
    tmp = load(fullfile(myOutDir, 'Measurements.mat'));
    myMeasurements = tmp.Measurements;
    fprintf('DONE in %.3f seconds\n', toc);

    % concatenate tables
    for iTable = 1:length(tableNames)
        tableName = tableNames{iTable};
        Measurements.(tableName) = [Measurements.(tableName); table2struct(myMeasurements.(tableName))];
    end
    Measurements.MotorData = [Measurements.MotorData, myMeasurements.MotorData];
end

% convert structure arrays to tables
for iTable = 1:length(tableNames)
    tableName = tableNames{iTable};
    fprintf('Converting structure array %s to table...\n', tableName);
    Measurements.(tableName) = struct2table(Measurements.(tableName));
end

% write tables to disk and create AllTables
AllTables = struct;
for iTable = 1:length(tableNames)
    tableName = tableNames{iTable};
    fprintf('Saving %s...\n', tableName);
    saveTable(Measurements.(tableName), tableName, {'csv'}, fullfile(outDir, allDir));
    AllTables.(tableName) = Measurements.(tableName);
end

% append SubjectTable
Measurements.Subjects = myMeasurements.Subjects;
AllTables.Subjects = myMeasurements.Subjects;

% append Cognition data
Measurements.CognitionData = myMeasurements.CognitionData;
AllTables.CognitionData = myMeasurements.CognitionData;

% append CognitionTable
Measurements.CognitionTable = myMeasurements.CognitionTable;
AllTables.CognitionTable = myMeasurements.CognitionTable;

% export Measurements structure to base workspace
assignin('base', 'Measurements', Measurements);

% export AllTables structure to base workspace
assignin('base', 'AllTables', AllTables);

fprintf('Saving AllTables structure to MAT file...\n');
tic
% save AllTables structure to MAT file
save(fullfile(outDir, allDir, 'AllTables.mat'), 'AllTables');
fprintf('DONE in %.3f seconds\n', toc);

fprintf('If necessary, save current state using ''saveState''\n');






