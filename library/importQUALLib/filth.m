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