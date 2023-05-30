% save Measurements structure to MAT file
fprintf('\t\t- Saving Measurements structure to MAT file...\n');
outDir = Measurements.outDir;
save(fullfile(outDir, 'Measurements.mat'), 'Measurements');
