function [SimInMN,Info] = getTime(Data,SimInMN,Info,UserPref)
%% Determine the timing variables (start, stop, freq) also in case of other import files.
% We have start and end time and frequency for Qualisys and each subsystem.
% QTMTime0 : Moreover for the Qualisys there is the startsample, in case the file was
% not exported entirely.
% There is the computational "internal" frequency
% For possible further input files (e.g. BVH from XSENS) there is a start and stop time and frequency
%
% SYNTAX
% [SimInMN,Info] = getTime(Data,SimInMN,Info,UserPref)
%
% INPUT
%     Data     (QTM struct) 
%     SimInMN  (MNData struct) output stream according MNData
%     Info     (struct) settings and constants
%     UserPref (struct) user preferences
%           
% OUTPUT
%     SimInMN      (MNData struct) output stream according MNData
%     Info         (struct) settings and constants
%
% See also: importQUAL
% 
% (c) 2022 by Predimo GmbH
% Website: http://www.predimo.com
% Author: Marc de Lussanet
% version 230203 (MdL & LK) cleanup of Info and UserPref; removed DoUseQTMStartTime

% Init
QTMFreq      = Data.FrameRate;
QTMTime      = ((1:Data.Frames)-1)/QTMFreq;
QTMStartTime = (Data.StartFrame-1)/QTMFreq;
TimeEnd      = QTMTime(end);

% If isCombine, the time base (params.time) comes from a different file (e.g., XSENS-MVNX).
if UserPref.isCombine
    QTMFreqSRint = getSampleRate(UserPref.time);
    QTMTime0     = QTMStartTime;
    TimeInternal = UserPref.time;
else
    QTMFreqSRint = UserPref.internalSR;
    QTMTime0     = 0;
    TimeInternal = 0: 1/QTMFreqSRint : TimeEnd;
end
NSinternal   = length(TimeInternal); % number of samples
QTMTimeSRint = 0: 1/QTMFreqSRint : TimeEnd;
NSintQTM     = length(QTMTimeSRint); % number of samples
% fill the Info structure
Info.TimeInternal          = TimeInternal(:); % force to column array
Info.QTMTime               = QTMTime;
Info.NSintQTM              = NSintQTM;
Info.QTMFreqSRint          = QTMFreqSRint;
Info.QTMTimeSRint          = QTMTimeSRint;
Info.QTMTime0              = QTMTime0;
Info.NSinternal            = NSinternal;

if ~isempty(SimInMN)
    % Time in the original QTM measurement
    SimInMN.signals.time.units = 's';
    SimInMN.signals.time.data  = TimeInternal(:); % force to column array
    SimInMN.meta.startTimeOrig = QTMStartTime;
    SimInMN.meta.stopTimeOrig  = QTMStartTime+TimeEnd;
end
end
