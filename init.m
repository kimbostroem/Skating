% clear workspace
clear
close all

inDir = '../Skating_In';
outDir = '../Skating_Out';
paramDir = '../Skating_In';

if ~isfolder(inDir)
    inDir = uigetdir('..', 'Select input folder');
end
if ~isfolder(outDir)
    outDir = uigetdir('..', 'Select input folder');
end
if ~isfolder(paramDir)
    paramDir = uigetdir('..', 'Select parameter folder');
end

% restore default path
restoredefaultpath;
% add library and subfolders to path
addpath(genpath('library'));
% add output folder to path
addpath(outDir);