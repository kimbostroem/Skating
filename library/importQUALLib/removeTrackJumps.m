function [PosClean, VelClean, Q, LargeJumps] = removeTrackJumps(Pos,Freq,Unit,MaxJump,SaveFigAs)
% Kinematic tracking of passive markers by multiple cameras leads to "jumps"
%   in the trajectories and spikes in the velocity.
% This function locates these spikes in the velocity signal and removes them.
% The threshold for jumps is detected automatically
%
% Minimal usage:
% PosClean = removeTrackJumps(Pos,Freq);
%
% Pos       [3D,NSamples] dataseries
% Freq      Sampling frequency
% ThrFact   (default 3) factor for determining the threshold
% MaxJump   Optional threshold to detect possible marker swaps and mis-identifications
% DoFigures Optional create explanatory figures
% DataName  Optional file name (for figs; MaxJump)
%
% PosClean   with steps removed
% VelClean   according velocity
% Q          Quality estimates : NSpikes, Error, CumLen
%            NSpikes = Number of spikes detected; 
%            Error   = Cumulative difference between Pos and PosClean (1D)
%            CumLen  = total length of all detected spikes
% LargeJumps Sample numbers of possible  marker flips (see MaxJump)

% Marc de Lussanet & Charlotte Le Mouel, Movement Science, WWU Muenster,  5.11.2020
%
% version  2 (12.11.2020) : handle jumps at beginning and end of time series
% version  3 (13.11.2020) : better threshold
% version  4 (23.11.2020) : clever HPcut;
% version  5 (25.11.2020) : cut each component individually
% version  6 (01.12.2020) : MaxJump; iterative Threshold
% version  7 (02.02.2021) : Check sign of spike (CosAA)
% version  8 (11.02.2021) : some small improvements (MaxJump)
% version  9 (22.02.2021) : Detect and interpolate clusters of jumps and
%                          fine-tuning of params based on Quality control (Q) 
% version 10 (28.02.2021) : Replaced VelHp with DotAA, in 3D
% version 11 (29.06.2021) : save jump figure
% version 12 (09.07.2021) : LargeJumps nargout was wrong


%% Constants:
narginchk(2,5);
if nargin < 3 || isempty(Unit),    Unit = 1; end  % Default 1=[m] for scaling the threshold
if nargin < 4 || isempty(MaxJump), MaxJump = 0; end  % if >0 warn at large steps [mm]
if nargin < 5,                     SaveFigAs = []; end  % do not save the large jump figure

% cutoff freq for low-pass filter for correcting possible drift of clean signal
LPcut    = 0.5; % .5 cutoff for smoothing the error signal (see below)

%% Initialize
Vel      = diff(Pos,1,2) * Freq;
Ct       = 1:length(Vel); % sample counter

%% Angle between subequent accelerations and length
Acc      = diff([Vel(:,1) Vel Vel(:,end)],1,2) * Freq;
DotAA    = dot(Acc(:,2:end),Acc(:,1:end-1),1); % dot product from sample to sample
JumpLen  = [real(sqrt(-DotAA * Freq^-4)) 0]; % Charlottes equation

%% the CLM Criterion
DotAA3D  = Acc(:,2:end) .* Acc(:,1:end-1);  % product from sample to sample
%QuotAA3D= Acc(:,2:end) ./ Acc(:,1:end-1);  % Quotient from sample to sample
Cutoffmm = 0.5/1000;
Cutoff   = -(Cutoffmm*Unit * Freq^2)^2;
% detect the jumps
%IsJump3D= (DotAA3D < Cutoff) & (QuotAA3D > -5) & (QuotAA3D < -1/5);
IsJump3D = (DotAA3D < Cutoff);
IsJump   = sum(IsJump3D)>0; % spikes in any direction


%% find clusters of jumps ("crappy parts") in each direction
Vsmth      = Vel;              % initialize
VelFill    = Vel;              % initialize
MovMnJumps = movmean(IsJump3D,10,2); % find denser groups of spikes by moving average
IsCluster  = zeros(size(MovMnJumps));
for Ax=1:3 % Loop the x,y,z coords
	Clusters = find(MovMnJumps(Ax,:)>0.4); % define clusters by threshold
	% for each sample above threshold, search Start and End of cluster
	St = zeros(size(Clusters)); % start 
	En = zeros(size(Clusters)); % end
	Next=0;
	for i=1:length(Clusters) 
		if Next>Clusters(i), continue; end % skip duplicates to speed up the search
		St(i) = find([0   MovMnJumps(Ax,1:Clusters(i))]==0,1,'last');
		En(i) = find([MovMnJumps(Ax,Clusters(i):end) 0]==0,1)+Clusters(i)-2;
		Next = En(i);
	end
	% handle clusters at begining and end of the time series
	if ~isempty(St) && St(1)  >En(1),   St=[1 St];                  end %#ok<*AGROW>
	if ~isempty(St) && St(end)>En(end), En=[En length(MovMnJumps)]; end
	En(En==0)=[]; % omit zeros
	St(St==0)=[]; % omit zeros
	% now create a boolean array that is based on the selected Starts and Ends 
	IsCluster(Ax,St)= 1;
	IsCluster(Ax,En)=-1;
	IsCluster(Ax,:) = cumsum(IsCluster(Ax,:));
	%% interpolate the clusters of spikes
	% interpolate the largest spikes in the velocity signal ...
	IsJumpAx = IsJump3D(Ax,:);
	Vsmth(Ax,IsJumpAx) = interp1(Ct(~IsJumpAx),Vsmth(Ax,~IsJumpAx),Ct(IsJumpAx)); %,'pchip',0);
	% ... and smooth the velocity
	Vsmth(Ax,:) = savitzkyGolayFilt(Vsmth(Ax,:),5,0,51,[],2);
	VelFill(Ax,IsCluster(Ax,:)==1) = Vsmth(Ax,IsCluster(Ax,:)==1);
end

%% Replace the velocity with the filled one
Vel = VelFill;

%% Analyse the jumps in each direction.
% If the jump is oblique, it is smaller in the affected coordinates
% allocate
VelClean  = Vel; % initialize: fill the jumps
PosClean  = Pos; % initialize: cleaned position signal
for Ax=1:3 % Loop the x,y,z coords
	% Select jumps that were not yet filled
	IsJumps = IsJump3D(Ax,:);
	% skip clusters because these were filled already
	IsJumps(IsCluster(Ax,:)>0)=0;
	% interpolate the transients
	VelClean(Ax,IsJumps) = interp1(Ct(~IsJumps),Vel(Ax,~IsJumps),Ct(IsJumps));
end
% Get the jump-free position by integrating
PosClean = cumsum([PosClean(:,1) VelClean/Freq],2);

%% Estimate the smooth drift error (the jumps may not add up to zero) 
Err      = filth(LPcut,Freq,2,   PosClean-Pos   ,'l',2);
PosClean = PosClean-Err; %  correct for the drift
% Recalculate the vel from the error-corrected pos (very minor changes)
VelClean = diff(PosClean,1,2) * Freq;

%% Quality control
% Count the no of spikes
Q.NSpikes = sum(IsJump);
Q.CumLen  = sum(JumpLen(IsJump));
Q.Error   = sum(sqrt(sum((PosClean-Pos).^2))); % the final error with respect to the raw signal
% fprintf('CumLen Error NSpikes %f %f %d\n',CumLen,Error*1000/lenght(Pos),Q.NSpikes);
% figure;hold on;plot(JumpLen);plot(IsJump/1000)

%% Detect possible marker swaps and mis-identifications
if MaxJump
	LargeJumps = find(JumpLen>MaxJump); % Samples in which the jump size is larger than MaxJump 
	LargeJumps(diff(LargeJumps)==1) = []; % Up and down flanks are found: remove the latter:
	if ~isempty(LargeJumps)
      if nargout<4 % warning in case LargeJumps is not taken as output argument
         warning('Possible marker swap: %d very large jumps detected (max %f)', ...
            length(LargeJumps),max(JumpLen));
      end
      if ~isempty(SaveFigAs)
         Fg=figure;hold on;
         plot(Pos');plot(LargeJumps,Pos(1,LargeJumps),'o');
         xlabel('samples');ylabel('position (m)');
         [~,Fname] = fileparts(SaveFigAs);
         title(sprintf('large jumps detected: %s',strrep(Fname,'_','\_')));
         %saveCurrentFigure(Fg,SaveFigAs,'');
         close(Fg);
      end
	end
end

end


%% ==========================================================================================


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

%% loop for naninterpolation and filtering
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
