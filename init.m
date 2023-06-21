% clear workspace
clear
close all

inDir = '../Skating_In';
motorDir = fullfile(inDir, 'Motorisch');
paramDir = fullfile(inDir, 'Parameter');
outDir = '../Skating_Out';

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));
% addpath(inDir);
% addpath(motorDir);
% addpath(paramDir);
% addpath(outDir);