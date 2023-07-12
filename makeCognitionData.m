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

% calc guided skating units
if ismember('ZusaetzlichSkaten', CognitionData.Properties.VariableNames) && ismember('Absage', CognitionData.Properties.VariableNames)
    addSkate = CognitionData.ZusaetzlichSkaten;
    missSkate = CognitionData.Absage;
    nRows = length(addSkate);
    Skating_reg = nan(nRows, 1);
    Skating_add = nan(nRows, 1);
    for iRow = 1:nRows

        % regular skating units
        switch missSkate(iRow)
            case 0
                Skating_reg(iRow) = 12;
            case 1
                Skating_reg(iRow) = 10;
            case 2
                Skating_reg(iRow) = 7;
            case 3
                Skating_reg(iRow) = 4;
            case 4
                Skating_reg(iRow) = 1;
            otherwise
                Skating_reg(iRow) = 12;
        end

        % additional skating units
        switch addSkate(iRow)
            case 0
                Skating_add(iRow) = 0;
            case 1
                Skating_add(iRow) = 24;
            case 2
                Skating_add(iRow) = 12;
            case 3
                Skating_add(iRow) = 6;
            case 4
                Skating_add(iRow) = 3;
            case 5
                Skating_add(iRow) = 1;
            otherwise
                Skating_add(iRow) = 0;
        end        
    end

    % total skating units
    Skating_tot = sum([Skating_reg, Skating_add], 2, 'omitnan');

    % append variables to table
    CognitionData.Skating_reg = Skating_reg;
    CognitionData.Skating_add = Skating_add;
    CognitionData.Skating_tot = Skating_tot;
end

% append CognitionData to Measurements structure
Measurements.CognitionData = CognitionData;

% export Measurements structure to base workspace
fprintf('\t\t- Exporting Measurements structure to base workspace...\n');
assignin('base', 'Measurements', Measurements);

% save current state
saveState;

end