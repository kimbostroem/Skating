function [ EMG,Freq,Labels ] = emgAusQualisys( Data, Label )
%% Marc de Lussanet, WWU Muenster
%% Version 1 (20.10.2017)
%     Case-Insensitive Label-selection
%     based on datenAusQualisys()
%% Version 2 (15.4.2019)
%     Enable Several Boards (e.g., Kistler and Noraxon)

%% this function gets the data structure (Dat) from mat exported from QTM (Qualisys track manager)
%% Output:
%   EMG    (optional) matrix of [nLabels , nSampEMG]
%   Freq   (optional) Measurement frequency
%   Labels (optional) cell list of all labels
%% Input:
%   Dat    Structure from QTM-exported .mat file
%   Label  (optional) cellstring of label names (if empty, then data of all labels are returned)
%% - Singleton dimensions are removed from Pos (so if just one coordinate of one marker is requested, an array is returned)
%% Examples:
%   [ ~,Freq,Labels  ] = emgAusQualisys( Dat ) % return only measurement frequency and list of labels
%   Pos = emgAusQualisys(Dat, {'Lab1',Lab2'})  % return EMG of 2 Markers (2 x nSamp)
%   Pos = emgAusQualisys(Dat, '')              % return EMG of all Markers (nMark x nSamp)

Repair = 1; % flag for repair of duplicate labels
Board  = 0;
EMG    = [];
Freq   = [];
Labels = [];
if ~isfield(Data,'Analog')
    return;
end
for i=1:length(Data.Analog)
    if strcmp(Data.Analog(i).BoardName, 'Noraxon')==1, Board=i; end
end
if ~Board
    %warning('No EMG Data found');
    return;
end

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

Nlab   = size(Label,2);
Nsam   = size(Data.Analog(Board).Data,2);
EMG    = NaN([Nlab Nsam]);

% copy by label, to asign the correct labels to each other
for i = 1 : Nlab
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
