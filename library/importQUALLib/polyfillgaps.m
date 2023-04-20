function [DataFilled,Gaps,EstErrors,Residuals] = polyfillgaps(DataWGaps,Freq,Tail,Degree,MaxGaplen,Extrapol,Unit,Range,Dim)
%% This function fills all NaN gaps of an array or matrix, using polyfit
%
% minimal usage:
% DataFilled = polyfillgaps(DataWithGaps);
% [DataFilled,Gaps,EstErrors] = polyfillgaps(DataWithGaps,Freq);
%
% Input:
%     Gaps :   - There can be any number of gaps; Gaps can be at the ends of the array. All are filled
%              - NGaps   = size(Gaps,1);
%              - Gapsize = Gaps(:,2)-Gaps(:,1)+1;
%     EstErrors    - estimated errors, the estimated mean error of the fill. Note, that the 
%                    median error is only about 1/10th of the mean and the 90%
%                    of the errors are smaller than the mean.
%     Residuals    -
%     DataWithGaps - [Nsamp, Cols] see Dim 
%     Freq         = [1/s] Measurment frequency (for EstErrors) 
%     Tail         - is the max length of the flanks left and right to the gap (Default 10) 
%     Degree       - is the Polynomial degree (Default 5)
%     MaxGaplen    - [s] default 0=no limit (if isempty(Freq), then [samples])
%     Extrapol     - [bool] (Default true)
%     Unit         - (Default m). Enter 0.001 if the unit is mm
%     Range        - [s] (optional) range to be filled
%     Dim          - (optional) Dimension along which to fill the gaps (default 1)
% 
% Outputs: 
%     DataFilled - Data with filled gaps
%     Gaps       - matrix with startpoints and endpoints of detected gaps
%     EstErrors  - Estimated Error of the fillings
%     Residuals  - Residuals of the fitted tails of the polynom
% 
% Local functions: savitzkyGolayFilt, savitzkyGolay
%
% See also: gapfilth
% 
% (c) 2017 by Movement Science, WWU Muenster
% Marc de Lussanet
% Version 18: 230119 (MdL) use min() for DegreeLongGap1, 2

%% Defaults
narginchk(1,9);
if nargin<3 || isempty(Tail) || isempty(Freq) || isempty(Tail)
	if nargin==1 || isempty(Freq)
		Freq = 1; % if Freq not given, then MaxGapLen is in samples
		Tail = 30;
	else
		Tail  = max(30,round(.15*Freq)); % 30 30 45 60 75 for 100 200 ... 500 Hz (60 ms)
		%Tail = max(10,round(.06*Freq)); % 10 12 18 24 30 for 100 200 ... 500 Hz (60 ms) 
	end
end
if nargin<4 || isempty(Degree),    Degree     = 8;    end
if nargin<5 || isempty(MaxGaplen), MaxGaplen  = 0;    end % 0 fill all
if nargin<6 || isempty(Extrapol),  Extrapol   = true; end
if nargin<7 || isempty(Unit),      Unit       = 1;    end
if nargin<8 || isempty(Range),     Range      = 0;    end % 0=fill complete range
if nargin<9 || isempty(Dim),       Dim        = 1;    end % Default: Column data

%% Parameters
% DegreeExtrap = Degree-1; % extrapolation with a lower degree (max 1)  
% Terminal  = 3*Tail; % longer flank for extrapolation
DegreeExtrap   = 0; % extrapolate with the same value ("padding")
Terminal       = 1; % "padding" on the basis of one sample, to prevent jumps.
DegreeLongGap1 = min(Degree,2); % lower degree for long gaps
DegreeLongGap2 = min(Degree,3); % lower degree for long gaps
if Freq==1 % frequency is not set
    LongGap    = 5*Tail;
else
    LongGap    = Freq; % Freq = 1 sec
end
PlotGaps       = 0; % 0: no
                    % 1: plot each gap into the current (existing) figure
                    % 2: a figure for each gap
                    % 3: a figure for each gap and save each figure

%% initialize
Gaps       = [];
Residuals  = [];
EstErrors  = [];

% Handle dimensions                                                 
if length(size(DataWGaps))>2
    warning('Too high matrix dimension: The length of size(Data)=%d > 2. (DATA NOT FILLED)',length(size(DataWGaps))); 
    return;
end
% Check if the data are transposed (default is transpose the data)
if Dim==1
    Data = DataWGaps';
	Transposed = true;
else
    Data = DataWGaps;
	Transposed = false;
end
DataFilled = Data;

[nDim,Ns] = size(Data);
if Ns<2*Tail+2
    warning('Data (Ns=%d, nDim=%d) too short for filling',Ns,nDim);
    if Ns<4 && nDim>4
        error('... the data have a highly unplausible dimensionality!')
    end
    DataFilled = DataWGaps;
    return;
end

Time    = 1:Ns;
% test if all dims have the same gaps
IsNan = isnan(Data);
Tst = sum(~IsNan,1,'omitnan');
Tst(Tst==0) = [];
if any(Tst~=nDim), fprintf('Not all dimensions of Data have the same gaps. Using the longest gaps\n'); end

%% Detect an order gaps
% find the beginnings of each gap
IsGap = any(IsNan,1);
St = [-1 find(IsGap)];   Df = [1 diff(St)];
St(Df==1)=[];	
% find the endings of the gaps
En = [find(IsGap) Ns+2]; Df = [diff(En) 1];
En(Df==1) = [];

% remove leading and trailing gaps
if ~Extrapol
	if St(1)==1
		St(1)=[];
		En(1)=[];
	end
	if En(end)==Ns
		St(end)=[];
		En(end)=[];
	end
end

% check the maximum gap length (s)
GapLen = (En-St+1);
if MaxGaplen^2 < 0.001
    Gaps = []; % fill all gaps
elseif MaxGaplen > 0
	En(GapLen/Freq > MaxGaplen) = [];
	St(GapLen/Freq > MaxGaplen) = [];
else
    return;
end

% check the range to be filled
if Range(1) && length(Range)>1
	for g=1:length(St)
		if St(g) > Range(2)*Freq || En(g) < Range(1)*Freq
			St(g) = 0;
			En(g) = 0;
		end
	end
end
St(St==0) = [];
En(En==0) = [];

% Analyse gaps
Gaps       = [St; En]';
NGaps      = length(St); 
Residuals  = zeros(1,NGaps);

%% only fit if gaps are present and data are not empty or only nans %171212
if NGaps>0 && ~isempty(Data) && ~all(isnan(Data),'all')
	%% Start and End of tails when the tails reach across nan entries
	% remove nans for fitting
	DataContinuous = Data; DataContinuous(:,isnan(Data(1,:))) = [];
	TimeContinuous = Time; TimeContinuous(  isnan(Data(1,:))) = [];
	% Start and End of Tails (STl, ETl) skipping the nan entries
	StEdge  = St-1;
	EnEdge  = En+1;
	STlCont = find(ismember([0 TimeContinuous Time(end)+1],StEdge)) - Tail-1;
	ETlCont = find(ismember([0 TimeContinuous Time(end)+1],EnEdge)) + Tail-1;
	StEdge(StEdge<=0) = 1;
	EnEdge(EnEdge>=Ns)= Ns;
	% Alternative flank for extrapolation
   if DegreeExtrap == 0
      % just one sample if degree == 0 (i.e. fill with the same value)
      if St(1)  ==1;  ETlCont(1)   = ETlCont(1)  -Tail;  end
      if En(end)==Ns; STlCont(end) = STlCont(end)+Tail;  end
   else
      % increase tail
      if St(1)  ==1;  ETlCont(1)   = ETlCont(1)  +Terminal;  end
      if En(end)==Ns; STlCont(end) = STlCont(end)-Terminal;  end
   end
	% truncate tails to 1:Ns
	STlCont(STlCont<1)  = 1;
	ETlCont(ETlCont>length(TimeContinuous)) = length(TimeContinuous);
	STlCont(STlCont>length(TimeContinuous)) = length(TimeContinuous)-1;
	STl = TimeContinuous(STlCont);
	ETl = TimeContinuous(ETlCont);
	
	%% loop through the gaps, if present
	for i=1: NGaps 
		
		Testl   = DataContinuous(:,STlCont(i):ETlCont(i));
		Timel   = TimeContinuous(  STlCont(i):ETlCont(i));
		TailFit = zeros(nDim,length(Timel));
		Interp  = zeros(nDim,Time(En(i))-Time(St(i))+1);
		
		%% loop through the dimensions
		for dd=1:nDim

			%% Deg -the degree of the polynomial
			if St(i)==1 || En(i)==Ns,     Deg=DegreeExtrap;
			elseif GapLen(i)>LongGap,     Deg=DegreeLongGap1;
			elseif GapLen(i)>LongGap/2,   Deg=min(DegreeLongGap2,Degree);
			else,                         Deg=Degree;
			end
			% If few points are fitted, then reduce Deg
			if Deg >= length(Timel),      Deg = length(Timel)-1;   end      %191010
			% centered and normalized fit, for better result (mu)
			[Par,~,mu]    = polyfit(Timel,Testl(dd,:),Deg);
			Interp(dd,:)  = polyval(Par,((Time(St(i)):Time(En(i)))-mu(1))/mu(2), mu);
			% residual for points in the fitted tails-fitted curve
			TailFit(dd,:) = polyval(Par,(Timel-mu(1))/mu(2), mu);
		end
		
		Residuals(i) = mean(rms(TailFit-Testl));
		%figure; plot((DataFilled)');plot((TailFit)')

		%% fill the fit into the curve: 
		DataFilled(:,St(i):En(i))=Interp;
		
		%% plot each gap
		if PlotGaps
			for dd=1:nDim
				%create figure for each gap plot filled data for each gap,  
				%otherwise plot into current figure
				if PlotGaps>=2
					str = sprintf('Gap %2d\n',i);
					H=figure; hold on; plot(Time,DataFilled(dd,:)); title(str); 
				end  %#ok<*UNRCH>
				% plot interpolated points 'o', and tails '+'
				TmIntp = St(i):En(i);
				plot(TmIntp,Interp(dd,:),'o'); plot(Timel,Testl,'+'); 
				% plot fitted curve
				x1 = linspace(Timel(1),Timel(end));	y1 = polyval(Par,(x1-mu(1))/mu(2), mu); 
				plot(x1,y1);
				% plot fitted tails
				P1   = []; P2   = []; %#ok<*NASGU>
				if -(En(i)+1 - ETl(i)) > Deg+1
					P1   = polyfit(Time(En(i)+1:ETl(i)),Data(dd,En(i)+1:ETl(i)),Degree);
				end
				if -(STl(i) - St(i)-1) > Deg+1
					P2   = polyfit(Time(STl(i):St(i)-1),Data(dd,STl(i):St(i)-1),Degree);
				end
				if ~isempty(P1), y1 = polyval(P1,x1); plot(x1,y1); end
				if ~isempty(P2), y1 = polyval(P2,x1); plot(x1,y1); end
				if PlotGaps==3, Name=sprintf('%d_gap_%02d',H.Number,i); savefig(Name);	end
			end
		end
	end
	
	if nargout>2 && nargin>1 && Freq>1
		%% estimate the error of the gapfilling. For this a regression model was developed
		% The simplest regression model weighs the standard deviation of the
		% velocity one second before and after the gap, and the length of the
		% gap in seconds. 
		% (interestingly, the product of the squareroots [interaction] correlates
		% positively with the error but both individual factors negatively)
		GapVel = -savitzkyGolayFilt(DataWGaps,3,1,21,[],2)*Freq*Unit;
		StDev  = zeros(1,length(St));
		for i=1:length(St)
			S1=max(1,St(i)-Freq); S2=min(Ns,En(i)+Freq);
			StDev(i) = sqrt( ... 
				std(GapVel(1,S1:S2),'omitnan')^2 + ...
				std(GapVel(2,S1:S2),'omitnan')^2 + ...
				std(GapVel(3,S1:S2),'omitnan')^2 );
		end
		p1 = [0.0060	-0.020	-0.033	0.12];
		GapLn       = (En-St+1)/Freq;
		x1          = sqrt(GapLn);
		x2          = sqrt(StDev);
		EstErrors   = p1(1) + p1(2)*x1 + p1(3)*x2 + p1(4)*x1.*x2;
		%% These are two more complex regression models,
		% which, however, do not seem to give better results
		VelPoly    = diff(DataFilled,1,2)*Freq;
		VelPolyTan = sqrt(sum(VelPoly.^2));
		VelTan     = sqrt(sum(GapVel.^2));
		MeanVels   = (VelTan(StEdge)+VelTan(EnEdge))/2;
		PathPoly   = cumsum([0 VelPolyTan])/Freq;
		p2         = [0.0017	-0.006	-0.007	-17.9	0.033	65];
		x3         = Residuals * Unit;
		ErrEstComplex = p2(1) + p2(2)*x1 + p2(3)*x2 + p2(4)*x3 + p2(5)*x1.*x2 + p2(6)*x1.*x3;
		p3         = [0.00092	-0.0023	0.0052	-42	-0.059	76.49	0.15	26	-35];
		PathPolyFilled= PathPoly(En) -  PathPoly(St);
		EstGapDetours = MeanVels.*(En-St+2) / Freq;
		x4         = abs(EstGapDetours-PathPolyFilled);
		ErrEstComplex2 = p3(1) + p3(2)*x1 + p3(3)*x2 + p3(4)*x3 + p3(5)*x4 + ...
			p3(6)*x1.*x3 + p3(7)*x1.*x4 + p3(8)*x2.*x3 + p3(9)*x3.*x4;
		EstErrors = ErrEstComplex2;
		
		% extrapolation is always suspect
		if St(1)==1
			EstErrors(1) = 1; 
		end
		if En(end)==Ns
			EstErrors(end) = 1; 
		end
		EstErrors(isnan(EstErrors)) = 1;
	end
end
if Transposed
	DataFilled = DataFilled';
end
end %function


%% ========================================================================


function y=savitzkyGolayFilt(x,N,DN,F,W,DIM)
%savitzkyGolayFilt Savitzky-Golay Filtering.
%   savitzkyGolayFilt(X,N,DN,F) 
%     X : filters the signal X using a Savitzky-Golay(polynomial) filter.  
%     N : The polynomial order, N, must be less than the
%     F : frame size, F, and F must be odd.  
%     DN specifies the differentiation order (DN=0 is smoothing). 
%     For a DN higher than zero, you'll have to scale the output by 1/T^DN to acquire the DNth 
%     smoothed derivative of input X, where T is the sampling interval. 
%     The length of the input X must be >= F.  If X is a matrix, the filtering is done on the columns of X.
%
%   Note that if the polynomial order N equals F-1, no smoothing
%   will occur.
%
%   savitzkyGolayFilt(X,N,DN,F,W) specifies a weighting vector W with
%   length F containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   savitzkyGolayFilt(X,N,DN,F,[],DIM) or savitzkyGolayFilt(X,N,DN,F,W,DIM)
%   operates along the dimension DIM.
%
%   See also savitzkyGolay, FILTER, sgolayfilt

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2009/08/11 15:47:54 $

narginchk(4,6);

% Check if the input arguments are valid
if round(F) ~= F, error(generatemsgid('MustBeInteger'),'Frame length must be an integer.'), end
if rem(F,2) ~= 1, error(generatemsgid('SignalErr'),'Frame length must be odd.'), end
if round(N) ~= N, error(generatemsgid('MustBeInteger'),'Polynomial order must be an integer.'), end
if N > F-1, error(generatemsgid('InvalidRange'),'The Polynomial order must be less than the frame length.'), end
if DN > N, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

if nargin < 5 || isempty(W)
	% No weighting matrix, make W an identity
	W = ones(F,1);
else
	% Check for right length of W
	if length(W) ~= F, error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
	% Check to see if all elements are positive
	if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
end

if nargin < 6, DIM = []; end

% Compute the projection matrix B
pp = fix(-F./2):fix(F./2);
B = savitzkyGolay(pp,N,DN,pp,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(generatemsgid('InvalidDimensions'),'Dimension specified exceeds the dimensions of X.')
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(generatemsgid('InvalidDimensions'),'The length of the input must be >= frame length.'), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on (note, this is different than in sgolayfilt,
% they had an optimization leaving out some transposes that is only valid
% for DN==0)
y(1:(F+1)/2-1,:) = fliplr(B(:,(F-1)/2+2:end)).'*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B(:,(F-1)./2+1),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = fliplr(B(:,1:(F-1)/2)).'*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
end


%% ========================================================================


function [fc, df] = savitzkyGolay(x,n,dn,x0,W,flag)
% Function:
%       Savitzky-Golay Smoothing and Differentiation Filter
%       The Savitzky-Golay smoothing/differentiation filter (i.e., the
%       polynomial smoothing/differentiation filter, or  the least-squares
%       smoothing/differentiation filters) optimally fit a set of data
%       points to polynomials of different degrees.
%       See for details in Matlab Documents (help sgolay). The sgolay
%       function in Matlab can deal with only symmetrical and uniformly
%       spaced data of even number.
%       This function presented here is a general implement of the sgolay
%       function in Matlab. The Savitzky-Golay filter coefficients for even
%       number, nonsymmetrical and nonuniformly spaced data can be
%       obtained. And the filter coefficients for the initial point or the
%       end point can be obtained too. In addition, either numerical
%       results or symbolical results can be obtained. Lastly, this
%       function is faster than MATLAB's sgolay.
%
% Usage:
%       [fc,df] = savitzkyGolay(x,n,dn,x0,flag)
%   input:
%       x    = the original data point, e.g., -5:5
%       n    = polynomial order
%       dn   = differentation order (0=smoothing),  default=0
%       x0   = estimation point, can be a vector    default=0
%       W    = weight vector, can be empty
%              must have same length as x0          default=identity
%       flag = numerical(0) or symbolical(1),       default=0
%
%   output:
%       fc   = filter coefficients obtained (B output of sgolay).
%       df   = differentiation filters (G output of sgolay).
% Notes:
% 1.    x can be arbitrary, e.g., odd number or even number, symmetrical or
%       nonsymmetrical, uniformly spaced or nonuniformly spaced, etc.
% 2.    x0 can be arbitrary, e.g., the initial point, the end point, etc.
% 3.    Either numerical results or symbolical results can be obtained.
% Example:
%       sgsdf([-3:3],2,0,0,[],0)
%       sgsdf([-3:3],2,0,0,[],1)
%       sgsdf([-3:3],2,0,-3,[],1)
%       sgsdf([-3:3],2,1,2,[],1)
%       sgsdf([-2:3],2,1,1/2,[],1)
%       sgsdf([-5:2:5],2,1,0,[],1)
%       sgsdf([-1:1 2:2:8],2,0,0,[],1)
% Author:
%       Diederick C. Niehorster <dcniehorster@hku.hk> 2011-02-05
%       Department of Psychology, The University of Hong Kong
%
%       Originally based on
%       http://www.mathworks.in/matlabcentral/fileexchange/4038-savitzky-golay-smoothing-and-differentiation-filter
%       Allthough I have replaced almost all the code (partially based on
%       the comments on the FEX submission), increasing its compatibility
%       with MATLABs sgolay (now supports a weight matrix), its numerical
%       stability and it speed. Now, the help is pretty much all that
%       remains.
%       Jianwen Luo <luojw@bme.tsinghua.edu.cn, luojw@ieee.org> 2003-10-05
%       Department of Biomedical Engineering, Department of Electrical Engineering
%       Tsinghua University, Beijing 100084, P. R. China
% Reference
%[1]A. Savitzky and M. J. E. Golay, "Smoothing and Differentiation of Data
%   by Simplified Least Squares Procedures," Analytical Chemistry, vol. 36,
%   pp. 1627-1639, 1964.
%[2]J. Steinier, Y. Termonia, and J. Deltour, "Comments on Smoothing and
%   Differentiation of Data by Simplified Least Square Procedures,"
%   Analytical Chemistry, vol. 44, pp. 1906-1909, 1972.
%[3]H. H. Madden, "Comments on Savitzky-Golay Convolution Method for
%   Least-Squares Fit Smoothing and Differentiation of Digital Data,"
%   Analytical Chemistry, vol. 50, pp. 1383-1386, 1978.
%[4]R. A. Leach, C. A. Carter, and J. M. Harris, "Least-Squares Polynomial
%   Filters for Initial Point and Slope Estimation," Analytical Chemistry,
%   vol. 56, pp. 2304-2307, 1984.
%[5]P. A. Baedecker, "Comments on Least-Square Polynomial Filters for
%   Initial Point and Slope Estimation," Analytical Chemistry, vol. 57, pp.
%   1477-1479, 1985.
%[6]P. A. Gorry, "General Least-Squares Smoothing and Differentiation by
%   the Convolution (Savitzky-Golay) Method," Analytical Chemistry, vol.
%   62, pp. 570-573, 1990.
%[7]Luo J W, Ying K, He P, Bai J. Properties of Savitzky-Golay Digital
%   Differentiators, Digital Signal Processing, 2005, 15(2): 122-136.
%
%See also:
%       sgolay, savitzkyGolayFilt

% Check if the input arguments are valid and apply defaults
narginchk(2,6);

if round(n) ~= n, error(generatemsgid('MustBeInteger'),'Polynomial order (n) must be an integer.'), end
if round(dn) ~= dn, error(generatemsgid('MustBeInteger'),'Differentiation order (dn) must be an integer.'), end
if n > length(x)-1, error(generatemsgid('InvalidRange'),'The Polynomial Order must be less than the frame length.'), end
if dn > n, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

% set defaults if needed
if nargin<6
	flag=false;
end
if nargin < 5 || isempty(W)
	% No weighting matrix, make W an identity
	W = eye(length(x0));
else
	% Check W is real.
	if ~isreal(W), error(generatemsgid('NotReal'),'The weight vector must be real.'),end
	% Check for right length of W
	if length(W) ~= length(x0), error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
	% Check to see if all elements are positive
	if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
	% Diagonalize the vector to form the weighting matrix
	W = diag(W);
end
if nargin<4
	x0=0;
end
if nargin<3
	dn=0;
end

% % prepare for symbolic output
% if flag
% 	x=sym(x);
% 	x0=sym(x0);
% end

Nx  = length(x);
x=x(:);
Nx0 = length(x0);
x0=x0(:);

if flag
	A=ones(length(x),1);
	for k=1:n
		A=[A x.^k]; %#ok<*AGROW>
	end
	df = inv(A'*A)*A';                          %#ok<MINV> % backslash operator doesn't work as expected with symbolic inputs, but the "slowness and inaccuracy" of this method doesn't matter when doing the symbolic version
else
	df = cumprod([ones(Nx,1) x*ones(1,n)],2) \ eye(Nx);
end
df = df.';

hx = [(zeros(Nx0,dn)) ones(Nx0,1)*prod(1:dn)];  % order=0:dn-1,& dn,respectively
for k=1:n-dn                                    % order=dn+1:n=dn+k
	hx = [hx x0.^k*prod(dn+k:-1:k+1)];
end

% filter coeffs
fc = df*hx'*W;
end