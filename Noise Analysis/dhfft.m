function [freq,p] = dhfft(data, ts, ys, seg, w, numpts)
% [freq,p] = dhfft(data, ts, ys, seg, w, numpts)
%
% dhfft takes an input data set data
% with a sampling rate ts (in sec)
% and vertical data scale ys
% averages a number of segments seg (use 1 for the entire trace)
% and smooths by an average number of points w (use 1 for no filtering)
% for an (optional) # of points numpts
% and returns a power spectrum 'p'
% and frequency scale 'freq'
%
% more memory-efficient 5/11/09

data = data*ys;
if nargin==5
    n = floor(length(data)/seg);
else
    n = floor(numpts/seg);
end


[p, freq] = pwelch(data, n, [], n, 1/ts);      % 50% overlap
%[p, freq] = pwelch(data, blackman(n), [], n, 1/ts);      % 50% overlap with blackman window
clear data


% adjustment of range required to avoid DC frequency explosion

if (w > 1)
    p=filtfilt(ones(1,w)/w,1,fliplr(p));
    p=fliplr(p(w+2:length(p)));
    freq=freq(2:(length(freq)-w));
%    p=filtfilt(ones(1,w)/w,1,p);
%    p=p(w+2:length(p));
%    freq=freq(w+2:length(freq));
else
    p=p(w+2:length(p));
    freq=freq(w+2:length(freq));
end