function Measurements = do2_makeInfo(Measurements)

if nargin < 1
    load('Measurements'); %#ok<LOAD>
end

% get folders
dataDir = Measurements.dataDir;
paramDir = Measurements.paramDir;
outDir = Measurements.outDir;

% load table containing subjects info
Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'));

% stages
stages = {'I', 'II', 'III'};

% list content of input folder into cell array
dirInfo = dir(dataDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxContact = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idxContact) = []; % remove folders and hidden files
nMeas = length(fpaths);

Measurements.Observations(nMeas, 1) = struct;
fprintf('Extract measurement info...\n');
ticAll = tic;
for iMeas = 1:nMeas
    fpath = fpaths{iMeas};
    [~, fileName, ~] = fileparts(fpath);
    if contains(fileName, 'ungueltig', 'IgnoreCase', true)
        fprintf('\t-> Skipping invalid file %s\n', fileName);
        continue
    end

    % report progress
    fprintf('\t-> %s\n', fileName);

    % split file name at underscores
    parts = strsplit(fileName, '_');

    % subject identity and code
    subjectName = parts{1};
    subjectCode = [parts{1}, '_', parts{2}];
    subjectCodes = [Subjects.Code_I, Subjects.Code_II, Subjects.Code_III];
    [subjectIdx, ~] = find(strcmp(subjectCodes, subjectCode));
    if isempty(subjectIdx)
        warning('Subject code %s not found -> skipping', subjectCode);
        continue
    end
    Measurements.Observations(iMeas).subjectName = string(subjectName); % store in Measurements.Observations(iMeas)
    Measurements.Observations(iMeas).subjectCode = string(subjectCode);

    % stage
    stageStr = parts{2};
    stage = find(strcmp(stages, stageStr), 1, 'first');
    Measurements.Observations(iMeas).stage = stage; % store in Measurements.Observations(iMeas)

    % subject properties at time of Measurements.Observations(iMeas)
    subjectProps = {'Height', 'Weight', 'Age', 'Date'};
    for iProp = 1:length(subjectProps)
        myStage = stage;
        while myStage > 0
            propName = sprintf('%s_%s', subjectProps{iProp}, stages{myStage});
            propValue = Subjects.(propName)(subjectIdx);
            if strcmp(subjectProps{iProp}, 'Date') || ~isnan(propValue)
                Measurements.Observations(iMeas).(lower(subjectProps{iProp})) = propValue;
                break
            end
            myStage = myStage-1;
        end
    end

    % intervention
    [~, subjectStages] = find(contains(subjectCodes, subjectName));
    switch stageStr
        case 'I'
            intervention = 0;
        case 'II'
            if max(subjectStages) == 3 % subject has been tested at three stages
                intervention = 0;
            else
                intervention = 1;
            end
        case 'III'
            intervention = 1;
    end
    Measurements.Observations(iMeas).intervention = intervention;

    % task
    task = parts{3};
    Measurements.Observations(iMeas).task = string(task);

    % side or Kraft
    if strcmp(parts{4}, 'Kraft')
        side = parts{5};
        trial = parts{6};
    else
        side = parts{4};
        trial = parts{5};
    end
    Measurements.Observations(iMeas).side = string(side);

    % jump position markers
    Measurements.Observations(iMeas).Beidbein_start = string(Subjects.Beidbein_start(subjectIdx));
    Measurements.Observations(iMeas).Beidbein_stop = string(Subjects.Beidbein_stop(subjectIdx));
    Measurements.Observations(iMeas).Einbein_start = string(Subjects.Einbein_start(subjectIdx));
    Measurements.Observations(iMeas).Einbein_stop = string(Subjects.Einbein_stop(subjectIdx));

    % trial number
    Measurements.Observations(iMeas).trial = str2double(trial);

    % store file name in Measurements.Observations(iMeas)
    Measurements.Observations(iMeas).fileName = string(fileName);
end
fprintf('Finished in %f s\n', toc(ticAll));

%% Export Measurements structure to base workspace

assignin('base', 'Measurements', Measurements);

%% Save Measurements structure to MAT file

fprintf('Saving Measurements to MAT file...\n');
save('Measurements', 'Measurements');
fprintf('Copying Measurements.mat to output folder...\n');
copyfile('Measurements.mat', outDir);
fprintf('DONE\n\n');

end


