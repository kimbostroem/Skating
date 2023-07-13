fprintf('\t\t- Saving Measurements structure to MAT file...\n');
outDir = evalin('base', 'outDir');
% export Measurements structure to base workspace
assignin('base', 'Measurements', Measurements);
% save Measurements structure to MAT file
save(fullfile(outDir, 'Measurements.mat'), 'Measurements');
