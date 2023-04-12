fprintf('Reset paths...\n');

restoredefaultpath; % restore default path

% add paths to full model
oldPaths = strsplit(path, pathsep);
addpath(genpath('library'));
paths = strsplit(path, pathsep);
for i = 1:(length(paths)-length(oldPaths))
    myPath = strrep(paths{i},'\','\\'); % escape backslashes
    msg = sprintf('Added path %s\n', myPath);
    fprintf(msg);
end
% add paths for client app
oldPaths = strsplit(path, pathsep);
addpath(genpath(fullfile('MyonardoPro'))); % include MyonardoPro app plus library
paths = strsplit(path, pathsep);
for i = 1:(length(paths)-length(oldPaths))
    myPath = strrep(paths{i},'\','\\'); % escape backslashes
    msg = sprintf('Added path %s\n', myPath);
    fprintf(msg);
end

clear oldPaths paths i myPath msg RESTOREDEFAULTPATH_EXECUTED