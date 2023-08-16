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
end

end


