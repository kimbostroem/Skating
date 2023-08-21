function makePlots

% load current state
Measurements = loadState();

outDir = evalin('base', 'outDir');

fprintf('Make plots...\n');
ticAll = tic;

% create figure output folders
figTypes = {'pdf', 'png'};
nFigTypes = length(figTypes);
figDirs = cell(nFigTypes, 1);
for iFigType = 1:nFigTypes
    figType = figTypes{iFigType};
    figDir = fullfile(outDir, ['figures_', figType]);
    if ~isfolder(figDir)
        mkdir(figDir);
    end
    figDirs{iFigType} = figDir;
end

nProc = 0; % init number of processed files
nMeas = length(Measurements.MotorData);
for iMeas = 1:nMeas

    fig = makePlot(iMeas, Measurements);

    % increment number of processed files
    nProc = nProc+1;
    
    % save figure
    for iFigType = 1:nFigTypes
        figType = figTypes{iFigType};
        figDir = figDirs{iFigType};
        outpath = fullfile(figDir, fileName);
        saveFigure(fig, outpath, figType);
    end
    close(fig);
end

fprintf('Finished plotting from %d datasets in %.3f s\n\n', nProc, toc(ticAll));

end