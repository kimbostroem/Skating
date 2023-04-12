doPlot = 1;
inDir = fullfile('..', 'Skating_In');
outDir = fullfile('..', 'Skating_Out');

% list content of input folder into cell array
dirInfo = dir(inDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
fpaths = strcat(fdirs, filesep, fnames);
fpaths([dirInfo.isdir]) = []; % remove folders
fpaths(startsWith(fpaths, {'.', '~'})) = []; % remove hidden files
nObs = length(fpaths);

observations = table;
nObs = 3;
for iObs = 1:nObs
   observation = table;
   fpath = fpaths{iObs};
   [~, fname, ~] = fileparts(fpath);
   if contains(fname, 'ungueltig', 'IgnoreCase', true)
       continue
   end
   parts = strsplit(fname, '_');   
   observation.subject = string(parts{1});   
   observation.task = string(parts{2});   
   observation.side = 'B';
   observation.time = 1;
   observation.trial = str2double(parts{end});
   if length(parts) > 3
       if strcmp(parts{3}, 'Post')
           observation.time = 2;
           observation.side = parts{4};
       else
           observation.time = 1;
           observation.side = string(parts{3});
       end
   end

   tmp = load(fpath);
   fields = fieldnames(tmp);
   Data = tmp.(fields{1});
   lengthScale = 0.001; % mm -> m
   % [Forces, ~, COP0, ~, Location, Analog] = kraftAusQualisys(Data, lengthScale);   
   [Forces, Frequency, COPs, ~, Location, Analog] = kraftAusQualisys(Data, lengthScale);
   COPs = COPs + squeeze(Location(:,4,:));
   % [COPs,~,Forces,Loaded] = getCOPfromAnalog(Analog,Forces,[],Location,Frequency,[],[],[],COPs);
   % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
   LoadThresh0 = 200;
   HumPeriod = Frequency/50;
   CutoffFrequency = 20;
   nPlates = size(Forces, 1);
   nDims = size(Forces, 2);
   nSamples = size(Forces,3);
   Loaded = false(size(Forces,1),1);
   Force = zeros(3, nSamples);
   COP = zeros(3, nSamples);
   Time = (0:nSamples-1)/Frequency;
   for plate = 1:nPlates
       if min(Forces(plate,3,:),[],'all') < -LoadThresh0
       % if Loaded(plate)
           % Loaded(plate) = true;
           Force_filt   = squeeze(Forces(plate,:,:));
           COP_filt   = squeeze(COPs(plate,:,:));
           % 50 Hz Netzbrumm entfernen
           for d = 1:nDims
               Force_filt(d,:) = periodicMedianFilter(Force_filt(d,:), HumPeriod);
               COP_filt(d,:) = periodicMedianFilter(COP_filt(d,:), HumPeriod);
           end
           % Glaetten
           Force_filt = nanfilth(CutoffFrequency,Frequency,2,Force_filt,'l',2);
           Force    = Force + Force_filt;
           COP_filt = nanfilth(CutoffFrequency,Frequency,2,COP_filt,'l',2);
           COP    = COP + COP_filt;
       end
   end
   Force = -Force; % force -> reaction force
   % COP = [-COP(1,:); -COP(2,:); COP(3,:)]; % rotate 180° on XY-plane to improve readability (negative values -> positive)

   switch observation.task
       case 'Einbein'
       case 'Balance'
       case 'Sprung'
       otherwise
           continue
   end
   % path length
   dt = median(diff(Time));
   T = Time(end) - Time(1);
   pathLength = sum(vecnorm(diff(COP, 1, 2), 2, 1), 2)/T;

   precision = sqrt(sum(var(COP,[],2)));


   if doPlot
       fig = setupFigure(800, 600, fname);
       subplot(2,1,1);       
       plot(Time',COP');
       title(sprintf('%s COP', fname), 'Interpreter', 'none');
       legend({'x', 'y', 'z'});
       xlabel('Time [s]');
       ylabel('Position [m]');
       xlim([Time(1), Time(end)]);
       subplot(2,1,2);       
       plot(Time',Force');
       title(sprintf('%s GRF', fname), 'Interpreter', 'none');
       legend({'x', 'y', 'z'});
       xlabel('Time [s]');
       ylabel('Force [N]');
       xlim([Time(1), Time(end)]);
       filePath = fullfile(outDir, fname);
       saveFigure(fig,filePath,'png');
       close(fig);
   end

   observations = [observations; observation]; %#ok<AGROW>
end

subjects = unique(observations.subject);

function [Filtered,Filter,FgHandle] = periodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum)
%% Remove power line interference ("mains hum") with harmonics from signal.
% This subtraction filter analyses the shape of the power line interference ("Hum noise") and
% subtracts it from the time series RawData.
% Method: moving estimate of Hum noise by periodic median of the high-pass filtered signal.
%
% This method is highly robust to any harmonic contributions, frequency and amplitude flutuations
% and transients in the signal. It performs really good at signal onset and offset.
%
% The optimal window is about 1 sec, but if the data have a strong frequency component at 1 Hz (e.g.
% tuck jumps), the filter can fail dramatically. A check and solution (halving the window) is
% implemented to deal with this case.
%
%% Syntax:
% Filtered = periodicMedianFilter(RawData, HumPeriod);
% [Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window);
% [Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window, Plot);
% [Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum);
%
%% mandatory parameters:
% RawData    (n x nsamples double) series of data with power line hum
% Humperiod  (double) the (whole) number of samples over which the hum repeats (e.g. 1000/50=20 if sample freq=1000)
%
%% optional parameters
% Window       [double] the number of periods on the sliding window (default 50; 0=entire range)
% Plot         [double] default 0. 1=A comparison plot and FFT; 2=also plot the filter
% PowerLineHum [double] the frequency of the hum (default 50 Hz)
%
%% output
% Filtered  (n x nsamples double) clean data
% Filter    (n x nsamples double) shape of the subtraction filter
% FgHandle  (handle) figure handle
%
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 11 (24.8.2022) Integrate a test for the goodness with fix in case the filter produces an
% error (see header). This used to be "periodicMedianFilterWithCatch"
% Version 12 (21.10.2022) fixed: crash if RawData contains nan values
% Version 13 (11.11.2022) remove gaps again after filtering

%% handle optional parameters
narginchk(2,5);
if nargin<3 || isempty(Window),       Window       = 50;    end
if nargin<4 || isempty(Plot),         Plot         = false; end
if nargin<5 || isempty(PowerLineHum), PowerLineHum = 50;    end

%% Test whether the current window works for the file
% The periodic median filter has almost always excellent performance, except when the signal has a
% main frequency with a period in the range of the window (i.e. typically 1 sec). This can be the
% case e.g., for tuck jumps.
% In that case (i.e. if a warning is produced), try reduce the window. If that does not help, no
% filtering is applied.
[Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum);
if any(Warn,'all')
    [Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window/2, Plot, PowerLineHum);
    if any(Warn,'all')
        warning('Reducing the window did not help: no hum filtering is applied')
        Filtered = RawData;
    end
end
end

function [Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum)
%% Subfunction that does the actual computation

%% handle optional parameters
narginchk(2,5);
if nargin<3 || isempty(Window),       Window       = 50;    end
if nargin<4 || isempty(Plot),         Plot         = false; end
if nargin<5 || isempty(PowerLineHum), PowerLineHum = 50;    end

%% Constants
HPCutFreq  = 20; % high pass frequency for constructing the subtraction filter
HPOrder    = 4 ; % (effective order is doubled by filtfilt)
MessFreq   = PowerLineHum * HumPeriod;
Sz         = size(RawData);
PlCh       = 1;  % channel that is plotted

%% Init
Signal  = RawData;
Filtered= RawData;
Filter  = [];
FgHandle= [];
Warn    = '';

%% error handling
if length(Sz)>2 || Sz(2)<Sz(1)
    error('Data must be shaped as Data(channels : timeseries)');
end
if length(RawData) < 2^nextpow2(PowerLineHum)/2
    disp('data too short for filtering'); return;
end
if any(all(isnan(RawData),2))
    % The data contains only NaN values in at least one dimension: there is nothing to be filtered
    return;
end
Gaps = isnan(RawData);
if any(isnan(RawData),'all')
    % Data contains NANs; these are filled with mean and removed again after filtering
    for Ch=1:Sz(1)
        RawData(Ch,Gaps(Ch,:)) = mean(RawData(Ch,:),'omitnan');
    end
    Signal  = RawData;
end

%% prepare the signal: remove offset and drift
% high-pass filter removes drifts and variations. 2nd Order (4th order effectively) and cutoff
% not too close to hum frequency, to prevent deformations.
Signal  = filth(HPCutFreq,MessFreq,HPOrder,Signal,'h');

%% if window longer than data, then simply take the entire signal (added 190902)
if HumPeriod*Window > length(Signal)
    Window=0;
end

%% create the filter. Either from entire signal or moving window
DoUseMovingWindow = Window ~= 0;
if DoUseMovingWindow
    for Ch=1:Sz(1)
        %% create a moving window of Window periods for which to compute a gradually changing filter
        % make a repeating array of the signal
        FiltMat  = reshape(repmat(Signal(Ch,:),[1,Window]),1,length(Signal)*Window);
        % make a matrix of this in which each column starts one period later in the signal
        RepeatLen= HumPeriod*Window;
        FiltMat  = reshape([FiltMat FiltMat(1:RepeatLen)], [length(Signal) + HumPeriod,Window]);
        % the filter is the median
        Median = median( FiltMat,2)';
        % the end of the Filter contains increasingly more of the start of the signal: omit this
        Len      = size(FiltMat,1) - RepeatLen;
        Median(Len+1 : end) = [];
        % Repeat the first and the last Periods and add those to the filter
        Vor      = repmat(Median(1:HumPeriod)        ,1,ceil((Window-1)/2));
        Nach     = repmat(Median(end-HumPeriod+1:end),1,2*Window);
        if Ch==1,      Filter=zeros(Sz(1),length(Vor)+length(Median)+length(Nach));      end
        Filter(Ch,:) = [Vor Median Nach];
        Filter(Ch,:) = Filter(Ch,:)-mean(Filter(Ch,:)); %zero-offset
    end
    % subtract
    Filtered = Filtered-Filter(:, 1:length(RawData));
    FilterRep= squeeze(Filter(PlCh, 1:length(RawData)));
    NoiseAmplitude = max(abs(Filter(:, 1:length(RawData))),[],2);
else
    for Ch=1:Sz(1)
        if Ch==1,      Filter=zeros(Sz(1),HumPeriod);      end
        %% chop the signal in pieces of one period and take the median as filter
        NPer           = floor(length(Signal)/HumPeriod);
        Repeat         = reshape(Signal(Ch,1:HumPeriod*NPer),HumPeriod,NPer);
        Filter(Ch,:)   = median(Repeat,2)';
        Filter(Ch,:)   = Filter(Ch,:)-mean(Filter(Ch,:));%zero-offset
        % repeat the filter for all periods
        FilterRep      = reshape(squeeze(repmat(Filter(Ch,:),1,NPer+1)),1,HumPeriod*(NPer+1));
        % ... and subtract it from the signal
        Filtered(Ch,:) = RawData(Ch,:) - FilterRep(1:length(Signal));
    end
    NoiseAmplitude = max(abs(Filter),[],2);
end
% remove the gaps again
Filtered(Gaps) = nan;
% report warnings
Error = abs(RawData-Filtered);
Warn = Error>2*NoiseAmplitude;
if any(Warn,'all')
    warning('The filter is not well suited for quickly repeated jumps (mean error = %f). Try using a shorter window',mean(Error,'all'));
end

%% figures if desired
if Plot
    RawCh = RawData(PlCh,:); FiltCh = Filtered(PlCh,:); %(this will plot the first channel only!)
    %% Frequency spectrum
    NFFT     = 2^nextpow2(length(FiltCh)); % Next power of 2 from length of y
    Freqs    = MessFreq/2*linspace(0,1,NFFT/2+1);
    Time     = (0:length(FiltCh)-1) / MessFreq;
    PowerVor = fft(RawCh -mean(RawCh),NFFT)/length(RawCh);
    PowerNach= fft(FiltCh-mean(RawCh),NFFT)/length(FiltCh);
    PowerVor = 2*abs(PowerVor( 1:NFFT/2+1));
    PowerNach= 2*abs(PowerNach(1:NFFT/2+1));

    FgHandle=figure('name','periodicMedianFilter');
    subplot(2,1,1);hold on; plot(Time,RawCh); plot(Time,FiltCh);
    title('Dataseries before (blue) and after PMS filter');   xlabel('time (s)');
    subplot(2,1,2);hold on; title('FFT spectrum before (blue) and after PMS filter')
    xlabel('frequency (s^{-1})'); ylabel('power');
    plot(Freqs,PowerVor);
    plot(Freqs,PowerNach);
    axis([0 inf 0 1.05*max(PowerVor(Freqs>PowerLineHum*0.9 & Freqs<PowerLineHum+1))]);
    if Plot==2
        figure;hold on; plot(RawData); plot(FilterRep);   title('dataseries / filter')
    end
end
end

%% =================================================================================================
%% =================================================================================================

function [Filtered]=filth(CutFreq,MessFreq,Order,Data,f,Dim)

%% Version 2 : 24.10.2017 Marc de Lussanet, WWU Muenster
%     high pass of low pass filter (f= 'high'  'stop' of 'low')
%     der effektive Order wird verdoppelt durch filtfilt
%% Version 3 : 15.4.2019
%     conventional order of Matrix (data in second dimension)

%% Checks
if nargin < 6
    Dim = 2; % default dimension of dataseries
end
Dims=size(Data);
if length(Dims)>3, error('Maximal dimensionality of the Data is 3.'); end
Dims(Dim)=[];

%% Filter type
if f=='h'
    [B,A]=butter(Order,(2*CutFreq)/MessFreq,'high');
elseif f=='s'
    [B,A]=butter(Order,(2*CutFreq)/MessFreq,'stop');
else % if f=='l' %% lowpass
    [B,A]=butter(Order,(2*CutFreq)/MessFreq);
end

%% Dimsionality issues
n = Dims(1);
if length(Dims)==2, m=Dims(2); else, m=0; end

%% loop for filtering
Filtered = Data;
for i=1:n
    if ~m
        if Dim==2
            Filtered(i,:) = filtfilt(B,A,Filtered(i,:));
        else
            Filtered(:,i) = filtfilt(B,A,Filtered(:,i));
        end
    else
        for j=1:m
            switch Dim
                case 1
                    Filtered(i,:,j) = filtfilt(B,A,Filtered(i,:,j));
                case 2
                    Filtered(:,i,j) = filtfilt(B,A,Filtered(:,i,j));
                otherwise
                    Filtered(i,j,:) = filtfilt(B,A,Filtered(i,j,:));
            end
        end
    end
end
end

%% =================================================================================================

function [Filtered]=nanfilth(CutFreq,MessFreq,Order,Data,f,Dim)
%% High pass oder low pass filter (f= 'high'  'stop' of 'low')
% "zero phase lag" (filtfilt): der effektive Order wird verdoppelt
% conventional order of Matrix (data in second dimension)
%
% (c) 2017, Movement Science, WWU Muenster
% Author: Marc de Lussanet
% Version 4 : 11.12.2020: filtered signal is returned with the same nan values as original

%% Checks
if nargin < 6
    Dim = 2; % default dimension of dataseries
end
Dims=size(Data);
if length(Dims)>3, error('Maximal dimensionality of the Data is 3.'); end
Dims(Dim)=[];

%% Filter type
if f=='h'
    [B,A]=butter(Order,(2*CutFreq)/MessFreq,'high');
elseif f=='s'
    [B,A]=butter(Order,(2*CutFreq)/MessFreq,'stop');
else % if f=='l' %% lowpass
    [B,A]=butter(Order,(2*CutFreq)/MessFreq);
end

%% Dimsionality issues
n = Dims(1);
if length(Dims)==2, m=Dims(2); else, m=0; end

%% loop for naninterpolation and filtering
Nans     = isnan(Data);
Filtered = Data;
for i=1:n
    if ~m
        if Dim==2
            Filtered(i,:) = naninterp(   Filtered(i,:));
            Filtered(i,:) = filtfilt(B,A,Filtered(i,:));
        else
            Filtered(:,i) = naninterp(   Filtered(:,i));
            Filtered(:,i) = filtfilt(B,A,Filtered(:,i));
        end
    else
        for j=1:m
            switch Dim
                case 1
                    Filtered(i,:,j) = naninterp(   Filtered(i,:,j));
                    Filtered(i,:,j) = filtfilt(B,A,Filtered(i,:,j));
                case 2
                    Filtered(:,i,j) = naninterp(   Filtered(:,i,j));
                    Filtered(:,i,j) = filtfilt(B,A,Filtered(:,i,j));
                otherwise
                    Filtered(i,j,:) = naninterp(   Filtered(i,j,:));
                    Filtered(i,j,:) = filtfilt(B,A,Filtered(i,j,:));
            end
        end
    end
end
Filtered(Nans) = nan;
end

%% =================================================================================================

function X = naninterp(X)
% Interpolate over NaNs, extrapolate with mean value.
% See INTERP1 for more info
M=mean(X,'omitnan');
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'linear',M);
end

