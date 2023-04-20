function [DataPose1,DataPose2] = selectStablePoseWindow(Data,MaxNSamp,Dim,Data2)
%% select a time window for the stable (t-)pose.
% Since gaps at the beginning are not filled, these must be omitted first. 
% If there are only nans in the series, return empty
% the first ten samples are ignored to ensure that stable pose has been reached
%
% SYNTAX
%   DataPose = selectStablePoseWindow(Data,MaxNSamp,Dim);
%
% INPUT
%     Data     (double) Data Nsamp x D or D x NSamp 
%     MaxNSamp (double) optional number of samples to select
%     Dim      (double) optional Dimension
%     Data2    (double) optional Data Nsamp x D or D x NSamp 
%
% OUTPUT
%     DataPose1 (double) Data Nsamp x D or D x NSamp 
%     DataPose2 (double) Data Nsamp x D or D x NSamp 
%
% Examples
%   DataPose = selectStablePoseWindow(Data,MaxNSamp,Dim);
%   DataPose = selectStablePoseWindow(Data,[],Dim);    % only remove leading and trailing nan values
%   DataPose = selectStablePoseWindow(Data,MaxNSamp);  % select the Dim automatically
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 221026 (MdL) first version
% version 221107 (MdL) correct error in Dim

%% Init
narginchk(1,4);
if nargin<2 || isempty(MaxNSamp)
    MaxNSamp = [];
end
% handle dimension: assume number of samples is longer than dimension
if nargin<3 || isempty(Dim)
    Sz = size(Data);
    [~,Dim] = min(Sz);
    % however, if the data is just one sample long, it will be 1x3 or 3x1, so the shortest 
    if max(Sz)==3 && min(Sz)==1
        [~,Dim] = max(Sz);
    end
end
if nargin<4
    Data2 = Data;
end

DataPose1 = [];
DataPose2 = [];

%% find and omit nan values (may be present at beginning and end)
Omit = all(isnan(Data),Dim) | all(isnan(Data2),Dim); % the values that are NaN
Omit(find(~Omit,10)) = true; % skip the first ten samples
if all(Omit) 
    return;
end

%% select (not more than) the first MaxNSamp
DataPose1 = Data;
DataPose2 = Data2;
if ~isempty(MaxNSamp)
    % use not more than MaxNSamp
    if Dim == 1
        DataPose1(:,Omit) = [];
        DataPose1(:,MaxNSamp+1:end) = [];
        DataPose2(:,Omit) = [];
        DataPose2(:,MaxNSamp+1:end) = [];
    else
        DataPose1(Omit,:) = [];
        DataPose1(MaxNSamp+1:end,:) = [];
        DataPose2(Omit,:) = [];
        DataPose2(MaxNSamp+1:end,:) = [];
    end
end
end
