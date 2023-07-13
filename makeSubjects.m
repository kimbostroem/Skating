function makeSubjects

fprintf('\nMaking Subjects...\n');

% load current state
Measurements = loadState();

paramDir = evalin('base', 'paramDir');
SubjectsTable = readtable(fullfile(paramDir, 'Subjects.xlsx'), 'TextType','string');

% delete excluded subjects
idx = (SubjectsTable.Excluded == 1);
SubjectsTable(idx,:) = [];
% delete "Excluded" and "WhyExcluded" columns
SubjectsTable.Excluded = [];
SubjectsTable.WhyExcluded = [];

% calc guided skating units
if ismember('ZusaetzlichSkaten', SubjectsTable.Properties.VariableNames) && ismember('Absage', SubjectsTable.Properties.VariableNames)
    addSkate = SubjectsTable.ZusaetzlichSkaten;
    missSkate = SubjectsTable.Absage;
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
    SubjectsTable.Skating_reg = Skating_reg;
    SubjectsTable.Skating_add = Skating_add;
    SubjectsTable.Skating_tot = Skating_tot;
end

% delete incomplete subjects
idx = startsWith(SubjectsTable.Subject, "PRCO");
SubjectsTable(idx,:) = [];

% append Subjects table to Measurements structure
Measurements.Subjects = SubjectsTable;

% save current state
saveState;

end