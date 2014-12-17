function [subsampled,f] = levelFilter(data,fs)

% Stephen Fleming
% 6/19/14

% levelFilter(data) takes in a vector data, sampled at fs, and applies a
% 10kHz Bessel filter, downsamples using median to a new sampling frequency
% of 10kHz, and then median filters extensively until convergence using a
% huge 10ms window.  This filter is rather slow.

%% Bessel filter at 1kHz

wn = 1000/(fs/2);
[b,a] = butter(4, wn, 'low');
f = filtfilt(b,a,data);

%% Downsample to 2kHz

n = round(fs/2000);
f = medfilt1(f,n,1e5);
%disp(size(f))
f = accumarray(1+floor((1:numel(f))/25.0)',f',[],@mean)';
subsampled = f;


%% Median filter

i = 1;
delta = 10;
while delta>0.0001
    f2 = medfilt1(f,51,1e5); % window is ~2ms: longer than one level, shorter than two
    delta = sum(abs(f2-f));
    %display([num2str(i) ' iterations. Convergence metric is ' num2str(delta)])
    i = i+1;
    f = f2;
end

figure()
plot(f)

end
