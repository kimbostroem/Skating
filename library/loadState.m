function Measurements = loadState()

if ~evalin('base','exist(''outDir'')')
    evalin('base', 'init');
end

outDir = evalin('base', 'outDir');

if evalin('base', 'exist(''Measurements'', ''var'')')
    fprintf('Importing Measurements structure from base workspace...\n');
    Measurements = evalin('base', 'Measurements');
elseif isfile(fullfile(outDir, 'Measurements.mat'))
    fprintf('Loading Measurements structure from disk...\n');
    tic
    % load Measurements structure
    tmp = load(fullfile(outDir, 'Measurements.mat'));
    Measurements = tmp.Measurements;
    fprintf('DONE in %.3f seconds\n', toc);
else
    fprintf('Creating new Measurements structure...\n');
    % init Measurements structure
    Measurements = struct;
    % delete content of output folder
    fprintf('Deleting content of output folder...\n');
    dirInfo = dir(outDir);
    fdirs = {dirInfo.folder}';
    fnames = {dirInfo.name}';
    idxNoDelete = matches(fnames, {'.', '..'});
    fdirs(idxNoDelete) = [];
    fnames(idxNoDelete) = [];
    fpaths = strcat(fdirs, filesep, fnames);
    for iPath = 1:length(fpaths)
        fpath = fpaths{iPath};
        if isfile(fpath)
            delete(fpath);
        else
            rmdir(fpath, 's');
        end
    end
end

end


