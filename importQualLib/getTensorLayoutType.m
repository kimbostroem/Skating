function [LayoutType,Sizes,Dim] = getTensorLayoutType(Data,Dim)
%% Get the type of layout of matrix or tensor data timeseries
% The problem of time series is that the level at which samples are defined, versus the direction 
% in which to find the time series is relevant for various operations such as filtering, 
% adding a sample at beginning or end, etc. There are five cases:
% [NSxn] (Dim==1, Layout==1)
% [nxNS] (Dim==2, Layout==2)
% [NSxnxm] (Dim==1, Layout==11)
% [nxNSxm] (Dim==2, Layout==12)
% [nxmxNS] (Dim==3, Layout==13)
% 
% Inputs:
%     Data      - data, dimensionality (length of size) less than 3
%     Dim       - [1 2 3] dimension representing the time series (default 1)
% Outputs:
%     LayoutType - the tyope of layout of the martrix or 3-tensor of data
%     Sizes      - size(Data)
%     Dim        - see above
% 
% Usage:
%     [Layout,Sizes,Dim] = getTensorLayoutType(Data,Dim);
%
% See also: diffSameLength, gapfilth, savGol
% 
% (c) 2022 by Predimo GmbH
% Marc de Lussanet
% Version: 220301 

%% init
% get the size of Data
Sizes = size(Data);
if length(Sizes)>3, error('Maximal dimensionality of the Data is 3.'); end

% Argument checks: if Dim is empty, set the default. 
narginchk(1,2);
if nargin<2 || isempty(Dim)
    % if Data is an array, it can be [1xNS] or [NSx1]. In the first case, we are sure that Dim=2
    IsOnexNSampArray = length(Sizes)==2 && Sizes(1)==1 && Sizes(2)>1;
    if IsOnexNSampArray
        Dim = 2; % An array can be 1xN, which is equivalent to Nx1: do this implicitly  
    else
        Dim = 1; % default dimension of dataseries
    end
end

%% assign one of five cases to LayoutType
if length(Sizes)==3 % dimensionality is 3
    LayoutType = 10+Dim; % Dim is the dimensionality along which the time series is provided
else % dimensionality is 2
    LayoutType = Dim; % Dim is the dimensionality along which the time series is provided
end

%% error check:
DefinedTypes = [1, 2, 11, 12, 13];% These are the allowed values for LayoutType
if ~ismember(LayoutType,DefinedTypes)
    error('Layout %d is not defined',LayoutType)
end
end
