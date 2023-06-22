function makeCognitionData()

% load current state
Measurements = loadState();

% load table containing subjects info
Subjects = Measurements.Subjects;

cogDir = evalin('base', 'cogDir');
CognitionData = readtable(fullfile(cogDir, 'Kognitiv.xlsx'), 'VariableNamesRange', '3:3', 'TextType','string');

% delete excluded subjects
subjectCodes = Subjects.Subject;
subjectVar = "Probandencode";
idx = ~cell2mat(cellfun(@(x) startsWith(x, subjectCodes), CognitionData.(subjectVar), 'un', 0));
CognitionData(idx,:) = [];


% append CognitionData to Measurements structure
Measurements.CognitionData = CognitionData;

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

end