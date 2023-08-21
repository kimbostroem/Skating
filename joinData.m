function joinData()

fprintf('\nJoining Data...\n');

allDir = evalin('base', 'allDir');
parentOutDir = evalin('base', 'parentOutDir');
outDir = fullfile(parentOutDir, allDir);
if ~isfolder(outDir)
    mkdir(outDir);
end
dirInfo = dir(parentOutDir);
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
    myOutDir = fullfile(parentOutDir, outDirs{iDir});

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
    saveTable(Measurements.(tableName), tableName, {'csv'}, outDir);
    AllTables.(tableName) = Measurements.(tableName);
end

%% Append tables to Measurements structure

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

%% Save tables to disk

% write Subjects table
fprintf('Saving Subjects table...\n');
saveTable(Measurements.Subjects, 'SubjectsTable', {'xlsx'}, outDir);

% write Motor table
fprintf('Saving Motor table...\n');
saveTable(Measurements.MotorTable, 'MotorTable', {'csv'}, outDir);

% write Cognition table
fprintf('Saving Cognition table...\n');
saveTable(Measurements.CognitionTable, 'CognitionTable', {'csv'}, outDir);

% write Skating table 
fprintf('Saving Skating table...\n');
saveTable(Measurements.SkatingTable, 'SkatingTable', {'csv'}, outDir);

% write Skating table subjectMean
fprintf('Saving Skating table averaged per subject...\n');
saveTable(Measurements.SkatingTable_subjectMean, 'SkatingTable_subjectMean', {'csv'}, outDir);

%% Save AllTables to disk

fprintf('Saving AllTables structure to MAT file...\n');
tic
% save AllTables structure to MAT file
save(fullfile(outDir, 'AllTables.mat'), 'AllTables');
fprintf('DONE in %.3f seconds\n', toc);

%% save Measurements structure to disk

% export outDir to base workspace, so that subsequent saveTable saves into
% the correct folder
assignin('base', 'outDir', outDir);
% save state
saveState






