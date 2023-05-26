% clear workspace
clear
close all

inDir = fullfile('..', 'Skating_In');
outDir = fullfile('..', 'Skating_Out');
paramDir = fullfile('..', 'Skating_In');

if ~isfolder(inDir)
    inDir = uigetdir('..', 'Select input folder');
end
if ~isfolder(outDir)
    outDir = uigetdir('..', 'Select input folder');
end
if ~isfolder(paramDir)
    paramDir = uigetdir('..', 'Select parameter folder');
end

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));
% add output folder to path
addpath(outDir);

response = 'n';
if isfile(fullfile(outDir, 'Measurements.mat'))
    response = input('Measurements.mat already exists in output folder.\nDo you want to load it and continue from where you left? ([y]/n) ', 's');
    if strcmp(response, 'n') % user wants to restart from scratch
        fprintf('\nAlright, let''s ditch the file and restart from scratch!\n\n');
    else
        fprintf('\nOK, let''s load the file and continue!\n\n');
    end
end

if strcmp(response, 'n') % user wants to restart from scratch
    fprintf('Creating new Measurements structure...\n');

    % init Measurements structure
    Measurements = struct;

    % delete content of output folder
    fprintf('Deleting content of output folder...\n');
    dirInfo = dir(outDir);
    fdirs = {dirInfo.folder}';
    fnames = {dirInfo.name}';
    idxFolders = [dirInfo.isdir];
    idxNoDelete = matches(fnames, {'.', '..'});
    fdirs(idxNoDelete) = [];
    fnames(idxNoDelete) = [];
    idxFolders(idxNoDelete) = [];
    fpaths = strcat(fdirs, filesep, fnames);
    for iPath = 1:length(fpaths)
        fpath = fpaths{iPath};
        if isfile(fpath)
            delete(fpath);
        else
            rmdir(fpath, 's');
        end
    end
    Measurements.inDir = inDir;
    Measurements.outDir = outDir;
    Measurements.paramDir = paramDir;

    % load table containing subjects info
    Subjects = readtable(fullfile(paramDir, 'Subjects.xlsx'));

    % stages
    stages = {'I', 'II', 'III'};

    % list content of input folder into cell array
    dirInfo = dir(fullfile(inDir, '*.mat'));
    fdirs = {dirInfo.folder}';
    fnames = {dirInfo.name}';
    idxExclude = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
    fnames(idxExclude) = [];
    fdirs(idxExclude) = [];
    fpaths = strcat(fdirs, filesep, fnames);
    nMeas = length(fpaths);

    fprintf('Extract measurement info...\n');
    ticAll = tic;
    for iMeas = 1:nMeas
        fpath = fpaths{iMeas};
        [~, fileName, ~] = fileparts(fpath);
        if contains(fileName, 'ungueltig', 'IgnoreCase', true)
            fprintf('\t-> Skipping invalid file %s\n', fileName);
            continue
        end

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

        % set flags
        Measurements.Observations(iMeas).doneData = 0;
        Measurements.Observations(iMeas).doneMetrics = 0;
        Measurements.Observations(iMeas).donePlots = 0;
    end

    % export Measurements structure to base workspace
    fprintf('Exporting Measurements structure to base workspace...\n');
    assignin('base', 'Measurements', Measurements);

    % save Measurements structure to MAT file
    fprintf('Saving Measurements to MAT file...\n');
    save('Measurements', 'Measurements');

    % move Measurements.mat to output folder
    fprintf('Moving Measurements.mat to output folder...\n');
    movefile('Measurements.mat', outDir);

    fprintf('Finished all in %f s\n', toc(ticAll));

end

clearvars -except Measurements

makeData;
makeMetrics;
makePlots;



