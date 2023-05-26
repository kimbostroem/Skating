function updateState
% update state with new MAT-Files in input folder

% load current state
loadState;

% stored file names
oldFileNames = {Measurements.Data.fileName};

% list content of input folder into cell array
dirInfo = dir(fullfile(inDir, '*.mat'));
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idxExclude = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fnames(idxExclude) = [];
fdirs(idxExclude) = [];

[newFileNames, newIdx] = setdiff(fnames, oldFileNames, 'sorted');

fpaths = strcat(fdirs, filesep, fnames);
nMeas = length(fpaths);

% save current state
saveState;

end