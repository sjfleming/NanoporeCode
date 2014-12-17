function [fAvg, powerAvg] = noiseplot(data,varargin)

% Stephen Fleming
% 07/15/2013

% Purpose: plot a figure that shows noise in a given signal as the
% amplitude at a given frequency.

% Syntax: noiseplot(data,legend,color)

% Inputs: ['data' is required, the rest are optional]
% data is an array of timeseries current data.
% Column 1 is time (seconds), column 2 is current (nA).
% varargin is standard syntax for other input options, which are color of
% plot and a legend entry

% Other variables:
% color is an array containing the 3 number hex color code for the plot color
% leg is the string containing the legend entry for the plot

% For MATLAB help document see "Fast Fourier Transform"

%% Defining variables
color = [0 0 0]; % if no color, set to black
if length(varargin) == 1 % one input only is a legend
    leg = varargin{1};
elseif length(varargin) == 2
    leg = varargin{1};
    color = varargin{2};
end
    
si = data(2,1) - data(1,1); % sampling interval in seconds
fs = 1/si; % samples per unit time

%% Creating the noise plot

% Calculate the Fourier transform
m = length(data(:,2));           % Window length
n = m;                              % Transform length
y = fft(data(:,2),n);            % DFT
f = (0:n-1)*(fs/n);                 % Frequency range
power = y.*conj(y)/n;               % Power of the DFT

% Cut off the first point, which seems anomalous
power(1) = [];

% Bin the frequency amplitudes and average them so it looks good on a
% log-log scale
bins = 100;
fAvg = [5,logspace(1,log10((n-1)*(fs/n)),bins)]; % bin centers
dividers = fAvg./2; % places about half way between bin centers
tempStart = 1;
for j = 1:bins+1
    tempEnd = find(f>dividers(j),1,'first'); % 
    powerAvg(j) = mean(power(tempStart:tempEnd-1)); % average around bin centers
    tempStart = tempEnd;
end

%% Plot the average power as a function of binned frequency
figure(6)
loglog(fAvg,powerAvg,'-','Color',color)
hold on

%% Finish up the figure
xlabel('Frequency (Hz)')
ylabel('Amplitude (nA^2/Hz)')
title('{\bf Noise}')
% Add to the legend if there is an entry specified
if length(varargin) >= 1
    [legend_h,object_h,plot_h,text_strings] = legend();
    text_strings{length(text_strings)+1} = leg;
    legend(text_strings)
end

end