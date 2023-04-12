inDir = fullfile('..', 'Skating_In');
outDir = fullfile('..', 'Skating_Out');

% list content of input folder into cell array
dirInfo = dir(inDir);
fpaths = {dirInfo.name}';
fpaths([dirInfo.isdir]) = []; % remove folders
fpaths(startsWith(fpaths, '.')) = []; % remove '.' and '..'
[~, fnames, ~] = fileparts(fpaths);
nObs = length(fnames);

observations = table;
for iObs = 1:nObs
   observation = table;
   fname = fnames{iObs};
   nameTokens = regexp(fname, '(\w*?)_(\w*)_(\d*)', 'tokens');
   condition = nameTokens{1}{2};   
   observation.subject = string(nameTokens{1}{1});
   observation.condition = string(condition);
   observation.trial = str2double(nameTokens{1}{3});
   conditionTokens = regexp(condition, '(\w*?)_(\w*)', 'tokens');
   observation.task = conditionTokens{1}{1};
   remainder = conditionTokens{1}{2};
   observations = [observations; observation]; %#ok<AGROW>
end

subjects = unique(observations.subject);
conditions = unique(observations.condition);
trials = unique(observations.trial);

