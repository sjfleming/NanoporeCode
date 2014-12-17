function [data,si] = importdata(file,startTime,endTime,quietFlag)

% Stephen Fleming
% July 1, 2013

% Purpose: Imports data from an Axopatch .abf file and puts it in a useful
% format for Matlab.

% Inputs:
% file = the file name, in full
% startTime = the experiment clock time where you want to start importing
% endTime = the experiment clock time where you want to stop importing

% Outputs:
% data(:,1) is time
% data(:,2) is the measured current
% si = the data sampling interval when data was taken, in microseconds

if quietFlag == 1
    [d,si] = abfloadQuiet(file,'start',startTime,'stop',endTime); % use abfloadQuiet.m to import the data
else
    [d,si] = abfload(file,'start',startTime,'stop',endTime); % use abfload.m to import the data
end

% t = linspace(startTime,endTime,length(d))'; % construct a time vector
ts = si*10^(-6);
t = [startTime:ts:startTime+(length(d)-1)*ts]'; % construct a time vector
data = [t, d];

end