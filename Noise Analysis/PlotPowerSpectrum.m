% Stephen Fleming
% Nov 14, 2013

% This script imports a specified data file and calls the function
% 'dhfft.m' (written by Dave Hoogerheide) to make a plot of the noise.

% Calls 'dhfft.m' and 'importdata.m', which calls 'abfload.m' and/or
% 'abfloadQuiet.m'

% Uses .abf files.  I made these programs for data taken in gap-free mode.

%% Definitions

clear all
%location = '/home/fleming/Documents/Physics/Golovchenko/Data/Axopatch/20130725/'; % file folder
location = '/Users/Stephen/Documents/Stephen/Research/Data/Biopore/20140602/';
name = '2014_06_02_0002.abf';

%% Load and calculate noise

% Import data
file = [location name]; % putting whole file name together
startTime = 12;
endTime = 24;
[data,si] = importdata(file,startTime,endTime,1);
ts = si*10^(-6); % sampling interval in seconds

% Calculate the noise power spectrum
w = 3; % number of points to average over
[freq, power] = dhfft(data(:,2),ts,1,1,w);

% Bin and average
bins = 500;
n = size(data,1);
fs = 1/ts;
fAvg = logspace(-1,log10((n-1)*(fs/n)),bins); % bin centers
firstCenter = logspace(log10(fAvg(1)),log10(fAvg(2)),3);
lastCenter = logspace(log10(fAvg(bins-1)),log10(fAvg(bins)),3);
dividers = [logspace(log10(firstCenter(2)),log10(lastCenter(2)),bins-1) fAvg(bins)]; % places half way between bin centers
tempStart = 1;
for j = 1:bins
    tempEnd = find(freq>dividers(j),1,'first'); % 
    powerAvg(j) = mean(power(tempStart:tempEnd-1)); % average around bin centers
    tempStart = tempEnd;
end

%% Plot the noise
figure(1)
loglog(freq,power,'Color',[0.8 0.8 0.8])
hold on
loglog(fAvg,powerAvg,'-b','LineWidth',1)
h = gca;
set(h,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',18)
ylabel('Power (nA^2/Hz)','FontSize',18)
title('Noise Power Spectrum +100mV','FontSize',18)
axis tight
ylim([1e-15 1e-4])
box on

%% Plot the average noise in a cumulative figure
% figure(2)
% loglog(fAvg,powerAvg,'-r','LineWidth',1)
% hold on
% h = gca;
% set(h,'FontSize',14)
% xlabel('Frequency (Hz)','FontSize',16)
% ylabel('Power (nA^2/Hz)','FontSize',16)
% title('Noise Power Spectra: Patch Model Cell','FontSize',16)
% axis auto
% box on
% % legend('Directly connected','Connected through pulse gen','Connected through pulse gen with compensation')

