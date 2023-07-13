SubjectsTable = readtable(fullfile(paramDir, 'Subjects.xlsx'), 'TextType','string');
CognitionTable = readtable(fullfile(cogDir, 'Kognitiv.xlsx'), 'VariableNamesRange', '3:3', 'TextType','string');

nSubjects = size(SubjectsTable, 1);
SubjectsTable.Absage = nan(nSubjects, 1);
SubjectsTable.ZusaetzlichSkaten = nan(nSubjects, 1);

for iSubject = 1:nSubjects
    Probandencode = SubjectsTable.Code_I{iSubject};

    idx = find(CognitionTable.Probandencode == string(Probandencode));
    if isempty(idx)
        continue
    end

    varName = "Absage";
    value = CognitionTable.(varName)(idx(1));
    SubjectsTable.(varName)(iSubject) = value;

    varName = "ZusaetzlichSkaten";
    value = CognitionTable.(varName)(idx(1));
    SubjectsTable.(varName)(iSubject) = value;   
end

writetable(SubjectsTable, fullfile(paramDir, 'Subjects_corr.xlsx'));