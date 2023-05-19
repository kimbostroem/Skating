% load library
addpath(genpath('library'));
% clear workspace
clear
close all

Measurements = struct;

% set folders
Measurements.dataDir = '/Users/kbostroem/sciebo/Skating/Skating_In';
Measurements.outDir = '/Users/kbostroem/sciebo/Skating/Skating_Out';
Measurements.paramDir = '/Users/kbostroem/sciebo/Skating/Auswertung';

% get info
Measurements = makeInfo(Measurements);

% get data
Measurements = makeData(Measurements);

% save everything
outpath = fullfile(Measurements.outDir, 'Measurements');
save(outpath, 'Measurements');


%
% switch Measurements.Info(iMeas).task
%     case 'Einbein'
%     case 'Balance'
%         % path length
%         pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2, 'omitnan')/sum(diff(Time));
%         precision = pathLength;
%         COPx = COP(1, idxContact);
%         COPy = COP(2, idxContact);
%         coefficients = polyfit(COPx, COPy, 1);
%         yBeamFcn = @(x) polyval(coefficients, x);
%         xFit = linspace(min(COPx), max(COPx), nSamplesContact);
%         yFit = yBeamFcn(xFit);
%         beamDistFcn = @(x, y) abs(y - yBeamFcn(x));
%         beamDist = beamDistFcn(COPx, COPy);
%     case 'Sprung'
%     otherwise
%         continue
% end
%
%