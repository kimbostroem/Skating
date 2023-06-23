function makeCognitionData()

fprintf('\nMaking Cognition Data...\n');


% load current state
Measurements = loadState();

cogDir = evalin('base', 'cogDir');
CognitionData = readtable(fullfile(cogDir, 'Kognitiv.xlsx'), 'VariableNamesRange', '3:3', 'TextType','string');

% convert subject code (including the stage number) to subject ID
subjectVar = "Probandencode";
subjectCodes = CognitionData.(subjectVar);
subjectIDs = arrayfun(@(x) regexprep(x, '(\w*)_\w*', '$1'), subjectCodes);
CognitionData.(subjectVar) = subjectIDs;
CognitionData = renamevars(CognitionData, 'Probandencode', 'Subject');

% delete excluded subjects
Subjects = Measurements.Subjects;
excludedSubjects = setdiff(CognitionData.Subject, Subjects.Subject);
idx = cell2mat(cellfun(@(x) ismember(x, excludedSubjects), CognitionData.Subject, 'un', 0));
CognitionData(idx,:) = [];
fprintf('Excluded %d subjects from Cognition table that did not appear in Subjects table\n', length(excludedSubjects));

% append CognitionData to Measurements structure
Measurements.CognitionData = CognitionData;

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

end