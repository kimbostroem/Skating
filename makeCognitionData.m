function makeCognitionData()

fprintf('\nMaking Cognition Data...\n');


% load current state
Measurements = loadState();

cogDir = evalin('base', 'cogDir');
CognitionTable = readtable(fullfile(cogDir, 'Kognitiv.xlsx'), 'VariableNamesRange', '3:3', 'TextType','string');

% convert subject code (including the stage number) to subject ID
subjectVar = "Probandencode";
subjectCodes = CognitionTable.(subjectVar);
subjectIDs = arrayfun(@(x) regexprep(x, '(\w*)_\w*', '$1'), subjectCodes);
CognitionTable.(subjectVar) = subjectIDs;
CognitionTable = renamevars(CognitionTable, 'Probandencode', 'Subject');

% delete excluded subjects
Subjects = Measurements.Subjects;
excludedSubjects = setdiff(CognitionTable.Subject, Subjects.Subject);
idx = cell2mat(cellfun(@(x) ismember(x, excludedSubjects), CognitionTable.Subject, 'un', 0));
CognitionTable(idx,:) = [];
fprintf('Excluded %d subjects from Cognition table that did not appear in Subjects table\n', length(excludedSubjects));

% append CognitionData to Measurements structure
Measurements.CognitionData = CognitionTable;

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

fprintf('If necessary, save current state using ''saveState''\n');


end