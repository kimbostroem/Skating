function makeSubjects

% load current state
Measurements = loadState();

paramDir = evalin('base', 'paramDir');
Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'), 'TextType','string');

% delete excluded subjects
idx = (Subjects.Excluded == 1);
Subjects(idx,:) = [];
% delete "Excluded" and "WhyExcluded" columns
Subjects.Excluded = [];
Subjects.WhyExcluded = [];

% delete incomplete subjects
idx = startsWith(Subjects.Subject, "PRCO");
Subjects(idx,:) = [];

% append Subjects table to Measurements structure
Measurements.Subjects = Subjects;

% save current state
saveState;

end