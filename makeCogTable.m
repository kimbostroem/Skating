function CogTable = makeCogTable

inDir = '../Skating_In';

cogFile = 'Kognitiv.xlsx';
cogPath = fullfile(inDir, cogFile);
opts = detectImportOptions(cogPath, 'VariableNamesRange', '3:3');
CogTable = readtable(cogPath, opts);

subjectsFile = 'Subjects.xlsx';
subjectsPath = fullfile(inDir, subjectsFile);
opts = detectImportOptions(subjectsPath, 'VariableNamesRange', '1:1');
Subjects = readtable(subjectsPath, opts);

variables_orig = CogTable.Properties.VariableNames;
variables_clean = [
    "Probandencode"
    "ADHS"
    "Stage"
    "Intervention"
    "Sex"
    "Age_yrs"
    "Height_cm"
    "Weight_kg"
    "Medikation"
    "AD_MW"
    "Hyp_MW"
    "D2_F__SW"
    "D2_BZO_SW"
    "D2_KL_SW"
    "Stroop_FWL_SW"
    "Stroop_FSB_SW"
    "Stroop_INT_SW"
    ];


% delete irrelevant columns
idx = ~ismember(variables_orig, variables_clean);
CogTable(:, idx) = [];

% delete irrelevant rows
subjectVar = 'Probandencode';
idx = cell2mat(cellfun(@(x) contains(x, 'PRCO'), CogTable.(subjectVar), 'un', 0));
CogTable(idx,:) = [];

% delete excluded subjects
excludedSubjectCodes = Subjects.Subject(Subjects.Excluded == 1);
idx = cell2mat(cellfun(@(x) startsWith(x, excludedSubjectCodes), CogTable.(subjectVar), 'un', 0));
CogTable(idx,:) = [];

end