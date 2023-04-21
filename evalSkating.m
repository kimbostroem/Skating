doPlot = 1;
inDir = fullfile('..', 'Skating_In');
outDir = fullfile('..', 'Skating_Out');

% list content of input folder into cell array
dirInfo = dir(inDir);
fdirs = {dirInfo.folder}';
fnames = {dirInfo.name}';
idx = startsWith(fnames', {'.', '~'}) | [dirInfo.isdir];
fpaths = strcat(fdirs, filesep, fnames);
fpaths(idx) = []; % remove folders and hidden files
nObs = length(fpaths);
nObs = 3;
observations = table;
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
   % COPs = COPs + squeeze(Location(:, 4, :));
   [COPs, ~, Forces, Loaded] = getCOPfromAnalog(Analog, Forces, [], Location, Frequency, [], [], [], COPs);
   % ignore data if Force in Z -direction is smaller than -LoadThresh Newton
   LoadThresh0 = 200;
   HumPeriod = Frequency/50;
   CutoffFrequency = 20;
   nPlates = size(Forces, 1);
   nDims = size(Forces, 2);
   nSamples = size(Forces, 3);
   Force = zeros(3, nSamples);
   COP = zeros(3, nSamples);
   Time = (0:nSamples-1)/Frequency;
   COP = squeeze(sum(COPs, 1, 'omitnan'));
   Force = squeeze(sum(Forces, 1, 'omitnan'));
   % for iPlate = 1:nPlates
   %     if min(Forces(iPlate, 3, :), [], 'all') < -LoadThresh0
   %     % if Loaded(plate)
   %         % Loaded(plate) = true;
   %         % Force_filt   = squeeze(Forces(iPlate, :, :));
   %         % COP_filt   = squeeze(COPs(iPlate, :, :));           
   %         % % 50 Hz Netzbrumm entfernen
   %         % for iDim = 1:nDims
   %         %     Force_filt(iDim, :) = periodicMedianFilter(Force_filt(iDim, :), HumPeriod);
   %         %     COP_filt(iDim, :) = periodicMedianFilter(COP_filt(iDim, :), HumPeriod);
   %         % end
   %         % Glaetten
   %         % Force_filt = nanfilth(CutoffFrequency, Frequency, 2, Force_filt, 'l', 2);
   %         Force    = Force + Force_filt;
   %         % COP_filt = nanfilth(CutoffFrequency, Frequency, 2, COP_filt, 'l', 2);
   %         COP    = COP + COP_filt;
   %     end
   % end
   Force = -Force; % force -> reaction force
   % COP = [-COP(1, :); -COP(2, :); COP(3, :)]; % rotate 180° on XY-plane to improve readability (negative values -> positive)

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

   precision = sqrt(sum(var(COP, [], 2)));


   if doPlot
       fig = setupFigure(800, 600, fname);
       subplot(2, 1, 1);       
       plot(Time', COP');
       title(sprintf('%s COP', fname), 'Interpreter', 'none');
       legend({'x', 'y', 'z'});
       xlabel('Time [s]');
       ylabel('Position [m]');
       xlim([Time(1), Time(end)]);
       subplot(2, 1, 2);       
       plot(Time', Force');
       title(sprintf('%s GRF', fname), 'Interpreter', 'none');
       legend({'x', 'y', 'z'});
       xlabel('Time [s]');
       ylabel('Force [N]');
       xlim([Time(1), Time(end)]);
       filePath = fullfile(outDir, fname);
       saveFigure(fig, filePath, 'png');
       close(fig);
   end

   observations = [observations; observation]; %#ok<AGROW>
end

subjects = unique(observations.subject);



