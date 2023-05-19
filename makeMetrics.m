function Measurements = makeMetrics(Measurements)

dataDir = Measurements.dataDir;
outDir = Measurements.outDir;
nMeas = length(Measurements.Info);

fprintf('Loading files...\n');
ticAll = tic;
Measurements.Data(nMeas, 1) = struct;
for iMeas = 1:2
    fileName = Measurements.Info(iMeas).fileName;
    fpath = fullfile(dataDir, fileName);
    fprintf('\t-> %s\n', fileName);

    

end

end