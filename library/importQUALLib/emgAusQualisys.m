function [EMG, Freq, Labels] = emgAusQualisys(Data, Label)
%% Get the data structure (Dat) from mat exported from QTM (Qualisys track manager)
% - Singleton dimensions are removed from Pos (so if just one 
%   coordinate of one marker is requested, an array is returned)
% 
% SYNTAX
% [EMG, Freq, Labels] = emgAusQualisys(Data, Label)
% 
% INPUT
%   Dat    (struct) from QTM-exported .mat file
%   Label  (cell)   label names (if empty, then data for all labels are returned)
% 
% OUTPUT
%   EMG    (double [nLabels x nSampEMG])
%   Freq   (double)   Measurement frequency
%   Labels (cell)     list of all labels
% 
% EXAMPLES
%   [ ~,Freq,Labels  ] = emgAusQualisys( Dat ) % return only measurement frequency and list of labels
%   Pos = emgAusQualisys(Dat, {'Lab1',Lab2'})  % return EMG of 2 Markers (2 x nSamp)
%   Pos = emgAusQualisys(Dat, '')              % return EMG of all Markers (nMark x nSamp)
% 
%
% See also: datenAusQualisys
% 
% (c) 2017 by Movement Science, WWU Muenster
% Author: Marc de Lussanet, WWU Muenster
% Version 171020 (MdL) Case-Insensitive Label-selection; based on datenAusQualisys()
% Version 190415 (MdL) Enable Several Boards (e.g., Kistler and Noraxon)
% Version 230421 (MdL) Allow Cometa EMG (in addition to Noraxon); revised header; some formatting

%% Init
Repair = 1; % flag for repair of duplicate labels
Board  = 0;
EMG    = [];
Freq   = [];
Labels = [];

% check if any EMG data are present in the data
if ~isfield(Data,'Analog')
    return;
end
for i=1:length(Data.Analog)
    if matches(Data.Analog(i).BoardName, {'Noraxon','Cometa'})==true, Board=i; end
end
if ~Board
    %warning('No EMG Data found');
    return;
end

% Get the data and parameters
Freq   = Data.Analog(Board).Frequency;
Labels = Data.Analog(Board).Labels;
% if just one argument then return Labels and freq
if nargin ==1
    return;
end
% if the function is called without a label list, then return everything
if isempty(Label)
    Label = Labels;
end

NLabels  = size(Label,2);
NSamples = size(Data.Analog(Board).Data,2);
EMG      = NaN([NLabels NSamples]);

% copy by label, to asign the correct labels to each other
for i = 1 : NLabels
    %Index  = contains(Labels, Label{i}, 'Ignorecase',true);
    Index  = strcmpi(Labels, Label{i});
    if sum(Index)>0
        if sum(Index)>1
            warning('datenAusQualisys: multiple definition of Label "%s"',Label{i});
            % if a label exists more than once, take the first occurrence
            if Repair==1
                x = find(Index,1);
                Index(:) = 0;
                Index(x) = 1;
            end
        end
        EMG(i,:) = Data.Analog(Board).Data(Index,:);
    end
end
end
