% function makePlots

% load current state
loadState;

outDir = Measurements.outDir;
nMeas = length(Measurements.Observations);

fprintf('Make plots...\n');
ticAll = tic;

%% Linear model evaluation
%% Load raw data

outDir = 'Results';
methods = {'Normal', 'Log'};
% methods = {'Normal'};
nMethods = length(methods);
% Observations = readtable('Observations.xlsx');
Observations = struct2table(Measurements.Observations);
% make categorical variables
categories = {'subject'};
for c = 1:length(categories)
    Observations.(categories{c}) = categorical(Observations.(categories{c}));
end
% responses = {'targetError', 'pathLength', 'meanJerk', 'meanJerkXY'}';
responses = {'targetError', 'pathLength'}';
nResponses = length(responses);
formulas = cell(nResponses, nMethods);
models = cell(nResponses, nMethods);
results = cell(nResponses, nMethods);
ps = nan(nResponses, nMethods);
Fs = nan(nResponses, nMethods);
DF1s = nan(nResponses, nMethods);
DF2s = nan(nResponses, nMethods);
etaPSquareValues = nan(nResponses, nMethods);
effectSizes = cell(nResponses, nMethods);
resultTables = [];
pNormals = nan(nResponses, nMethods);
hNormals = nan(nResponses, nMethods);
if ~isfolder(outDir)
    mkdir(outDir);
end

% restrict task
Observations = Observations(Observations.task == "Balance", :);

%% Data Plot

idx1 = matches(Observations.Properties.VariableNames, responses);
Before = cell2mat(table2array(Observations(Observations.stage == 1 & Observations.intervention == 0, idx1)));
After = cell2mat(table2array(Observations(Observations.stage == 2 & Observations.intervention == 1, idx1)));
% Responses  

figW = 800;
figH = 400;
fig = figure('Position', [0, 0, figW, figH]);

% mean values and std
% members in column (2nd) dimension, groups in row (1st) dimension
nMembers = 2;
nGroups = nResponses;
X = 1:nResponses;
Y = [nanmean(Before); nanmean(After)]';
Yerr = [nanstd(Before); nanstd(After)]';
hBar = bar(X, Y);
errX = nan(nGroups, nMembers);
errY = nan(nGroups, nMembers);
for iMember = 1:nMembers
    errX(:, iMember) = bsxfun(@plus, hBar(1).XData, [hBar(iMember).XOffset]);
    errY(:, iMember) = hBar(iMember).YData;
    for iGroup = 1:nGroups
        text(errX(iGroup, iMember), errY(iGroup, iMember)+Yerr(iGroup, iMember), sprintf('%.2f', Y(iGroup, iMember)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
hold on
errorbar(errX, errY, Yerr, '.r');
hold off
title('Data plot');
legend({'Before', 'After'}, 'Location', 'best');
xticklabels(responses);
ylabel('(multiple units)');

% save figure
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'PaperUnits', 'points');
set(fig, 'PaperSize', [figW figH]);
set(fig, 'renderer', 'painters');
print(fig, fullfile(outDir, 'DataPlot.pdf'), '-dpdf', sprintf('-r%.0f', 300));