function Stream = extractEMG(Stream,Data,Info)
%%
%
% Version 210831 outsourced from importQUAL
% Version 210902 - FIXED: added missing parameters
%                - checked behavior
% Version 211027 Implemented MNData structure

[ EMGRaw,EMGFreq,EMGLabels ] = emgAusQualisys(Data,'');
if isempty(EMGRaw)
    return;
end
    
DoPlot    = 0;

DoRemoveEMGHum = Info.params.DoRemoveEMGHum;
HumFreq        = Info.params.HumFreq;
EMGhighpass    = Info.params.EMGhighpass;
EMGlowpass     = Info.params.EMGlowpass;
EMGfiltorder   = Info.params.EMGfiltorder;
NSintQTM       = Info.NSintQTM;
QTMFreqSRint   = Info.QTMFreqSRint;

if any(any(isnan(EMGRaw))) % cope with nan values  % MdL 210111
    fprintf('EMG data contains NANs; Fix: these will be filled with mean\n');
    for Ch=1:size(EMGRaw,1)
        EMGRaw(Ch,isnan(EMGRaw(Ch,:))) = mean(EMGRaw(Ch,:),'omitnan');
    end
end

Temp            = EMGRaw; %#ok<NASGU>
HumPeriod       = EMGFreq / HumFreq;          % sample length of one period power line hum
EMGFiltered     = EMGRaw;                     % init
if DoRemoveEMGHum % remove power line hum. This is advisable in wired EMG, only if hum is present
    % periodic median subtraction filter
    EMGFiltered = periodicMedianFilter(EMGFiltered, HumPeriod, [],0,HumFreq);
end
% remove movement artifacts, rectify and smooth
EMGFiltered = preprocessEMG(EMGFiltered,EMGFreq,EMGhighpass,EMGlowpass,EMGfiltorder);
[nChannels,~] = size(EMGRaw);             % number of EMG channels
% resample to internal frequency
EMGRaw      = resample(EMGRaw,     QTMFreqSRint,EMGFreq,'dimension',2);
EMGFiltered = resample(EMGFiltered,QTMFreqSRint,EMGFreq,'dimension',2);
if length(EMGRaw)>NSintQTM
    EMGRaw      = EMGRaw(     :,1:NSintQTM);
    EMGFiltered = EMGFiltered(:,1:NSintQTM);
end

% check for rounding errors
if size(EMGRaw,2)==NSintQTM
    % fine
elseif size(EMGRaw,2)==NSintQTM+1 % too long: remove last sample
    EMGRaw      = EMGRaw(:,1:end-1);
    EMGFiltered = EMGFiltered(:,1:end-1); % too short: duplicate last
elseif size(EMGRaw,2)==NSintQTM-1
    EMGRaw      = [EMGRaw      EMGRaw(:,end)];
    EMGFiltered = [EMGFiltered EMGFiltered(:,end)];
else
    fprintf('Error in downsampling of EMG ns=%d should be %d\n',size(EMGRaw,2),NSintQTM);
end

% compose the emg substructure for the stream
for i=1:nChannels
    field = matlab.lang.makeValidName(EMGLabels{i});
    Stream.EMG(i).id   = i;
    Stream.EMG(i).name = field;
    Stream.EMG(i).dynamic.raw.units = 'mV';
    Stream.EMG(i).dynamic.raw.data = EMGRaw(i,:)';
    Stream.EMG(i).dynamic.EMG.units = 'mV';
    Stream.EMG(i).dynamic.EMG.data = EMGFiltered(i,:)';
end

if DoPlot
    ch = 1; %#ok<UNRCH> 
    Traw = (0:length(Temp)-1)/EMGFreq;
    TimeInternal = (0:length(EMGRaw)-1)/EMGFreq;
    figure;hold on;
    plot(Traw,Temp(ch,:))
    plot(TimeInternal,EMGRaw(ch,:))
    plot(TimeInternal,EMGFiltered(ch,:),'linewidth',2)
    xlabel('time (s)');ylabel('Amplitude');title(sprintf('function extractEMG() : EMG plot for channel %d',ch))
    legend('Raw','Resampled','smooth rectified');
end

end
