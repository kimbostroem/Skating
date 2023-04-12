function [Filtered]=gapfilth(CutFreq,MessFreq,Order,Data,FilterType,Dim,EdgeMethod)
%% Butterworth filter without phase shift, highpass, lowpass or bandstop, that deals with gaps 
%% containing of NaNs or Zeros.
% There are two methods: 
%   - piecewise only smoothes each piece between the NaN-gaps. 
%   - interp1 fills the gaps, smoothes and deletes the gaps again.
% Gaps are defined from Zeros and NaNs:
%   - Zeros: to avoid coincidental gaps, ALL dimentions must be zero at the sample of the gap.
%   - NaNs : the filtfilt function only works if there aren't any NaN values, so samples in which 
%            there is any NaN are considered as gap.
% 
% Syntax:
%     Filtered = gapfilth(CutFreq,MessFreq,Order,Data);
%     Filtered = gapfilth(CutFreq,MessFreq,Order,Data,'',[],EdgeMethod);
% 
% Inputs:
%     CutFreq    - (double) cutoff frequency
%     MessFreq   - (double) measurement frequency
%     Order      - (double) filter order
%     Data       - (double matrix or 3-tensor) with Dim defining the direction of the Samples
%     FilterType - (char) Default 'l' = low-pass; or 's' (band stop) or 'h' (high pass)
%     Dim        - (double) Default 2; dimension on which to apply the filter
%     EdgeMethod - (char*) Default 'piecewise' or 'interp1'
% Outputs:
%     Filtered   - filtered data
% 
% Other m-files required: getTensorLayoutType
% Local functions: gapinterp
% 
% See also: savGol, diffSameLength, filth, filterData
% 
% (c) 2016 by Movement Science, WWU Muenster
% Marc de Lussanet
% Version 7, 230123 (MdL) initialize output parameter "Filtered"

%% Init
narginchk(4,7);
if nargin < 5 || isempty(FilterType),          FilterType = 'l'; end 
if nargin < 6 || isempty(Dim),        Dim = []; end
if nargin < 7 || isempty(EdgeMethod), EdgeMethod = 'piecewise'; end

%% do nothing if the data contain only nan values
Filtered = Data;
Isnan = isnan(Data);
if all(Isnan,'all')
    return;
end

%% Filter type
if CutFreq == 0 || CutFreq == Inf % do not filter
    A = 1;
    B = 1;
elseif FilterType == 'h'
    [B,A] = butter(Order,(2*CutFreq)/MessFreq,'high');
elseif FilterType == 's'
    [B,A] = butter(Order,(2*CutFreq)/MessFreq,'stop');
else % if f=='l' (Default)
    [B,A] = butter(Order,(2*CutFreq)/MessFreq,'low');
end

%% Dimensionality issues
% get the size and the layout type of the data (which one of five possibilities)
[Layout,Sizes,Dim] = getTensorLayoutType(Data,Dim);
% the dimentionalities that are not the time series
SizesNonDim = [Sizes 1]; % elongate, to deal with dimensionality of 2
SizesNonDim(Dim)=[]; % remove the dimensionalty of the time series 

%% Detect gaps 
Iszero = Data==0; % find all zeros
% get the indices of the non- timeseries dimensionalities
NonDim=1:length(Sizes); NonDim(Dim) = [];
% A gap is when ANY nan or when ALL zero (see explanation above)
IsGap = squeeze(all(Iszero,NonDim) | any(Isnan,NonDim)); 
IsGap = IsGap(:)';
NSamples = length(IsGap); 

%% loop for naninterpolation and filtering
MinPieceLen = 2*(Order+1)+1; % filtfilt cannot handle very short pieces 
if NSamples<MinPieceLen
    return;
end
Filtered = Data;
switch EdgeMethod
    case 'piecewise'
        %% find pieces that are separated by nans or zeros
        % find the beginnings of each segment
        St = [-1 find(~IsGap)];      Df = [1 diff(St)]; St(Df==1) = [];
        % find the endings of the segment
        En = [find(~IsGap) NSamples+2]; Df = [diff(En) 1]; En(Df==1) = [];
        % do not filter very short pieces: remove those
        PieceLengths = En-St+1;
        St(PieceLengths<=MinPieceLen) = [];
        En(PieceLengths<=MinPieceLen) = []; 
        %% loop pieces and dimensions
        for Piece = 1:length(St) % pieces
            Range = St(Piece):En(Piece); % the range of the current piece
            % loop the dimensions
            for i=1:SizesNonDim(1)
                for j=1:SizesNonDim(2)
                    % the five possible Layouts
                    switch Layout
                        case 1,  Filtered(Range,i) = filtfilt(B,A,Filtered(Range,i));
                        case 2,  Filtered(i,Range) = filtfilt(B,A,Filtered(i,Range));
                        case 11, Filtered(Range,i,j) = filtfilt(B,A,squeeze(Filtered(Range,i,j)));
                        case 12, Filtered(i,Range,j) = filtfilt(B,A,squeeze(Filtered(i,Range,j)));
                        case 13, Filtered(i,j,Range) = filtfilt(B,A,squeeze(Filtered(i,j,Range)));
                        otherwise
                        error('piecewise: Layout %d is not defined',Layout)
                    end
                end
            end
        end
    case 'interp1'
        %% interpolate the gaps consisting of nans or zeros and smooth the entire series
        for i=1:SizesNonDim(1)
            for j=1:SizesNonDim(2)
                switch Layout
                    case 1
                        Series = gapinterp(Filtered(:,i),IsGap);
                        Filtered(:,i)  = filtfilt(B,A,Series);
                    case 2
                        Series = gapinterp(Filtered(i,:),IsGap);
                        Filtered(i,:)  = filtfilt(B,A,Series);
                    case 11
                        Series = squeeze(Filtered(:,i,j));
                        Series = gapinterp(Series,IsGap);
                        Filtered(:,i,j)  = filtfilt(B,A,Series);
                    case 12
                        Series = squeeze(Filtered(i,:,j));
                        Series = gapinterp(Series,IsGap);
                        Filtered(i,:,j)  = filtfilt(B,A,Series);
                    case 13
                        Series = squeeze(Filtered(i,j,:));
                        Series = gapinterp(Series,IsGap);
                        Filtered(i,j,:)  = filtfilt(B,A,Series);
                    otherwise
                        error('interp1: Layout %d is not defined',Layout)
                end
            end
        end
        % remove the filled pieces again
        Filtered(Isnan) = NaN;
        Filtered(Iszero) = 0;
    otherwise
        error('Method "%s" is not defined',EdgeMethod)
end
end

%% =================================================================================================

function X = gapinterp(X,IsGap)
% Interpolate over NaNs, do not extrapolate.
X(IsGap) = interp1(find(~IsGap), X(~IsGap), find(IsGap));
% pad the beginning and end with constant value
if isnan(X(1)),   Sm=find(~isnan(X),1);        X(1:Sm) = X(Sm); end
if isnan(X(end)), Sm=find(~isnan(X),1,'last'); X(Sm:end) = X(Sm); end
end
