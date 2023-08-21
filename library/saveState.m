fprintf('Saving Measurements structure to MAT file...\n');
tic
outDir = evalin('base', 'outDir');
% export Measurements structure to base workspace
assignin('base', 'Measurements', Measurements);
% save Measurements structure to MAT file
save(fullfile(outDir, 'Measurements.mat'), 'Measurements', '-v7.3');
fprintf('DONE in %.3f seconds\n', toc);
