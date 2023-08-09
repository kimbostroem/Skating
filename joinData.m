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

AllTables = struct;
AllTables.MotorTable = struct([]);
AllTables.CognitionTable = struct([]);
AllTables.SkatingTable = struct([]);
AllTables.SkatingTable_subjectMean = struct([]);

for iDir = 1:length(outDirs)
    myOutDir = fullfile(outDir, outDirs{iDir});

    % load Measurements structure
    fprintf('Loading Measurements structure from %s...\n', myOutDir);
    tic    
    tmp = load(fullfile(myOutDir, 'Measurements.mat'));
    Measurements = tmp.Measurements;
    fprintf('DONE in %.3f seconds\n', toc);

    % concatenate tables
    tableNames = fieldnames(AllTables);
    for iTable = 1:length(tableNames)
        tableName = tableNames{iTable};
        AllTables.(tableName) = [AllTables.(tableName); table2struct(Measurements.(tableName))];
    end
end

% convert structure arrays to tables
tableNames = fieldnames(AllTables);
for iTable = 1:length(tableNames)
    tableName = tableNames{iTable};
    fprintf('Converting structure array %s to table...\n', tableName);
    AllTables.(tableName) = struct2table(AllTables.(tableName));
end

% append Subject table
AllTables.Subjects = Measurements.Subjects;

% write tables to disk
tableNames = fieldnames(AllTables);
for iTable = 1:length(tableNames)
    tableName = tableNames{iTable};
    fprintf('Saving %s...\n', tableName);
    saveTable(AllTables.(tableName), tableName, {'csv'}, fullfile(outDir, allDir));
end

% export AllTables structure to base workspace
assignin('base', 'AllTables', AllTables);

fprintf('Saving AllTables structure to MAT file...\n');
tic
% save AllTables structure to MAT file
save(fullfile(outDir, allDir, 'AllTables.mat'), 'AllTables');
fprintf('DONE in %.3f seconds\n', toc);





